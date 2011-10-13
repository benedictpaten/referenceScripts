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
            int32_t *i = stHash_search(sequencesToAdjacencyLengths, metaSequence);
            assert(i != NULL);
            i[0] -= segment_getLength(segment);
            assert(i[0] >= 0);
        }
    }
    block_destructInstanceIterator(instanceIterator);
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
        assert(i[0] >= 0);
        stHash_insert(sequencesToAdjacencyLengths, metaSequence, i);
    }
    flower_destructSequenceIterator(sequenceIt);
    FILE *fileHandle = fopen(outputFile, "w");
    getMAFs(flower, fileHandle, removeSegmentLengths);
    sequenceIt = flower_getSequenceIterator(flower);
    while((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        MetaSequence *metaSequence = sequence_getMetaSequence(sequence);
        int32_t *i = stHash_search(sequencesToAdjacencyLengths, metaSequence);
        assert(i != NULL);
        assert(i[0] >= 0);
        fprintf(fileHandle, "%i\t%i\t%s\t%s\n", i[0], sequence_getLength(sequence), event_getHeader(sequence_getEvent(sequence)), sequence_getHeader(sequence));
    }
    flower_destructSequenceIterator(sequenceIt);
    fclose(fileHandle);
    flower = NULL;
    cactusDisk_destruct(cactusDisk);

    return 0;
}
