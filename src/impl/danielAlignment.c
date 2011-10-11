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

const char *sampleEventString;

static void printDanielBlock(Block *block, FILE *fileHandle) {
    //Now get the column
    Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
    Segment *segment;
    Segment *referenceSegment = NULL;
    Segment *sampleSegment = NULL;
    while ((segment = block_getNext(instanceIterator)) != NULL) {
        if (strcmp(event_getHeader(segment_getEvent(segment)), referenceEventString) == 0 && segment_getSequence(
                segment) != NULL) {
            if (referenceSegment != NULL) {
                goto end;
            }
            referenceSegment = segment;
        }
        if (strcmp(event_getHeader(segment_getEvent(segment)), sampleEventString) == 0 && segment_getSequence(segment)
                != NULL) {
            if (sampleSegment != NULL) {
                goto end;
            }
            sampleSegment = segment;
        }
    }

    if (referenceSegment != NULL && sampleSegment != NULL) {
        bool b = segment_getStrand(referenceSegment);
        referenceSegment = b ? referenceSegment : segment_getReverse(referenceSegment);
        sampleSegment = b ? sampleSegment : segment_getReverse(sampleSegment);
        char *referenceSeq = segment_getString(referenceSegment);
        char *sampleSeq = segment_getString(sampleSegment);
        fprintf(fileHandle, "%s %i 1 %s %i %i %s %s\n", sequence_getHeader(segment_getSequence(referenceSegment)),
                segment_getStart(referenceSegment), sequence_getHeader(segment_getSequence(sampleSegment)),
                segment_getStart(sampleSegment), segment_getStrand(sampleSegment), referenceSeq, sampleSeq);
        free(referenceSeq);
        free(sampleSeq);
    }
    end: block_destructInstanceIterator(instanceIterator);
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "danielAlignment");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    sampleEventString = otherReferenceEventString;
    getMAFs(flower, fileHandle, printDanielBlock);

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    flower = NULL;
    cactusDisk_destruct(cactusDisk);

    return 0;
}
