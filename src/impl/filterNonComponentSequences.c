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

stSortedSet *freeSequences;

static void getDisconnectedSequences(Block *block, FILE *fileHandle) {
    Block_InstanceIterator *instanceIterator1 = block_getInstanceIterator(block);
    Segment *segment;
    bool b = 1;
    while ((segment = block_getNext(instanceIterator1)) != NULL) {
        Event *event = segment_getEvent(segment);
        if (strcmp(event_getHeader(event), referenceEventString) == 0) {
            if (segment_getSequence(segment) != NULL && sequence_getLength(segment_getSequence(segment)) >= 1000000) {
                b = 0;
                break;
            }
        }
    }
    block_destructInstanceIterator(instanceIterator1);
    if(b) {
        instanceIterator1 = block_getInstanceIterator(block);
        while ((segment = block_getNext(instanceIterator1)) != NULL) {
            Event *event = segment_getEvent(segment);
            if (strcmp(event_getHeader(event), referenceEventString) != 0) {
                if (segment_getSequence(segment) != NULL) {
                    stSortedSet_insert(freeSequences, (void *)sequence_getHeader(segment_getSequence(segment)));
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
    if (otherReferenceEventString == NULL) {
        otherReferenceEventString = stString_copy(referenceEventString);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<disconnectedSequences>\n");
    freeSequences = stSortedSet_construct3((int (*)(const void *, const void *))strcmp, NULL);
    getMAFs(flower, fileHandle, getDisconnectedSequences);
    stSortedSetIterator *it = stSortedSet_getIterator(freeSequences);
    char *cA;
    while((cA = stSortedSet_getNext(it)) != NULL) {
        fprintf(fileHandle, "%s ", cA);
    }
    stSortedSet_destructIterator(it);
    fprintf(fileHandle, "</disconnectedSequences>\n");

    st_logInfo("Finished writing out the disconnected sequences\n");
    fclose(fileHandle);

    flower = NULL;
    cactusDisk_destruct(cactusDisk);

    //while(1);

    return 0;
}