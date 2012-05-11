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

stHash *sequencesToAdjacencyLengths;

static void removeSegmentLengths(Block *block, FILE *fileHandle) {
    //Now get the column
    Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
    Segment *segment;
    while ((segment = block_getNext(instanceIterator)) != NULL) {
        Sequence *sequence = segment_getSequence(segment);
        if(sequence != NULL) {
            MetaSequence *metaSequence = sequence_getMetaSequence(sequence);
            stList *list = stHash_search(sequencesToAdjacencyLengths, metaSequence);
            int32_t *i = stList_get(list, 0);
            assert(i != NULL);
            i[0] -= segment_getLength(segment);
            assert(i[0] >= 0);
            stList_append(list, segment);
        }
    }
    block_destructInstanceIterator(instanceIterator);
}

Segment *getOtherSegment(Segment *segment) {
    Block_InstanceIterator *it = block_getInstanceIterator(segment_getBlock(segment));
    Segment *otherSegment;
    while((otherSegment = block_getNext(it)) != NULL) {
        if(otherSegment != segment) {
            assert(segment_getReverse(segment) != otherSegment);
            block_destructInstanceIterator(it);
            return otherSegment;
        }
    }
    st_errAbort("Did not find a segment");
    block_destructInstanceIterator(it);
    return NULL;
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "sequenceCoverages");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    Flower_SequenceIterator *sequenceIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    sequencesToAdjacencyLengths = stHash_construct();
    while((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        MetaSequence *metaSequence = sequence_getMetaSequence(sequence);
        int32_t *i = st_malloc(sizeof(int32_t));
        i[0] = metaSequence_getLength(metaSequence);
        stList *list = stList_construct();
        assert(i[0] >= 0);
        stList_append(list, i);
        stHash_insert(sequencesToAdjacencyLengths, metaSequence, list);
    }
    flower_destructSequenceIterator(sequenceIt);
    FILE *fileHandle = fopen(outputFile, "w");
    getMAFs(flower, fileHandle, removeSegmentLengths);
    sequenceIt = flower_getSequenceIterator(flower);
    while((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        MetaSequence *metaSequence = sequence_getMetaSequence(sequence);
        stList *list = stHash_search(sequencesToAdjacencyLengths, metaSequence);
        assert(list != NULL);
        int32_t *i = stList_get(list, 0);
        assert(i != NULL);
        assert(i[0] >= 0);
        if(i[0] > 0 && sequence_getLength(sequence)-i[0] < 400 && sequence_getLength(sequence)-i[0] > 0) {
            fprintf(fileHandle, "%i\t%i\t%s\t%s\t%i", sequence_getLength(sequence)-i[0], sequence_getLength(sequence), event_getHeader(sequence_getEvent(sequence)), sequence_getHeader(sequence), stList_length(list)-1);
            for(int32_t j=1; j<stList_length(list); j++) {
                Segment *segment = stList_get(list, j);
                fprintf(fileHandle, "\t%s", segment_getString(segment));
                Segment *otherSegment = getOtherSegment(segment);
                fprintf(fileHandle, "\t%s", segment_getString(otherSegment));
                char *cA = segment_getString(segment);
                char *cA2 = segment_getString(otherSegment);
                int32_t l = 0;
                for(int32_t k=0; k<segment_getLength(segment); k++) {
                    if(toupper(cA[k]) != toupper(cA2[k]) && toupper(cA[k]) != 'N' && toupper(cA2[k]) != 'N') {
                        l++;
                    }
                }
                fprintf(fileHandle, "\t%i\t%i", block_getInstanceNumber(segment_getBlock(segment)), l);
            }
            fprintf(fileHandle, "\n");
        }
    }
    flower_destructSequenceIterator(sequenceIt);
    fclose(fileHandle);
    flower = NULL;
    cactusDisk_destruct(cactusDisk);

    return 0;
}
