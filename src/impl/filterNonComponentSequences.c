/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>

#include "sonLib.h"
#include "cactus.h"
#include "substitutions.h"
#include "adjacencyClassification.h"
#include "scaffoldPaths.h"
#include "referenceCommon.h"
#include "cactusMafs.h"

stSortedSet *connectedSequences;
int32_t minCoordinate = INT32_MAX;
int32_t maxCoordinate = 0;

static void getConnectedSequences(Block *block, FILE *fileHandle) {
    Block_InstanceIterator *instanceIterator1 = block_getInstanceIterator(block);
    Segment *segment;
    bool b = 0;
    while ((segment = block_getNext(instanceIterator1)) != NULL) {
        Event *event = segment_getEvent(segment);
        if (strcmp(event_getHeader(event), referenceEventString) == 0) {
            if (segment_getSequence(segment) != NULL && sequence_getLength(segment_getSequence(segment)) >= 1000000) {
                b = 1;
            }
        }
        if (strcmp(event_getHeader(event), otherReferenceEventString) == 0) {
            assert(segment_getSequence(segment) != NULL);
            segment = segment_getStrand(segment) ? segment : segment_getReverse(segment);
            assert(segment_getStart(segment) <= segment_getStart(segment_getReverse(segment)));
            if (segment_getStart(segment) - sequence_getStart(segment_getSequence(segment)) < minCoordinate) {
                minCoordinate = segment_getStart(segment) - sequence_getStart(segment_getSequence(segment));
            }
            if (segment_getStart(segment) + segment_getLength(segment) - sequence_getStart(segment_getSequence(segment)) > maxCoordinate) {
                maxCoordinate = segment_getStart(segment) + segment_getLength(segment) - sequence_getStart(segment_getSequence(segment));
            }
        }
    }
    block_destructInstanceIterator(instanceIterator1);
    if (b) {
        instanceIterator1 = block_getInstanceIterator(block);
        while ((segment = block_getNext(instanceIterator1)) != NULL) {
            Event *event = segment_getEvent(segment);
            if (strcmp(event_getHeader(event), referenceEventString) != 0) {
                if (segment_getSequence(segment) != NULL) {
                    stSortedSet_insert(connectedSequences, (void *) sequence_getHeader(segment_getSequence(segment)));
                }
            }
        }
        block_destructInstanceIterator(instanceIterator1);
    }
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "filterNonComponentSequences");
    assert(referenceEventString != NULL);
    assert(otherReferenceEventString != NULL);

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    connectedSequences = stSortedSet_construct3((int(*)(const void *, const void *)) strcmp, NULL);
    getMAFs(flower, fileHandle, getConnectedSequences);
    fprintf(fileHandle, "<disconnectedSequences minOtherReferenceCoordinate=\"%i\" maxOtherReferenceCoordinate=\"%i\" referenceEventString=\"%s\" otherReferenceEventString=\"%s\">\n", minCoordinate, maxCoordinate, referenceEventString, otherReferenceEventString);
    Flower_SequenceIterator *sequenceIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        if (stSortedSet_search(connectedSequences, (void *) sequence_getHeader(sequence)) == NULL) {
            fprintf(fileHandle, "%s ", sequence_getHeader(sequence));
        }
    }
    flower_destructSequenceIterator(sequenceIt);
    fprintf(fileHandle, "</disconnectedSequences>\n");

    st_logInfo("Finished writing out the disconnected sequences\n");
    fclose(fileHandle);

    flower = NULL;
    cactusDisk_destruct(cactusDisk);

    //while(1);

    return 0;
}
