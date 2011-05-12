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
int32_t totalSites = 0;
double totalCorrect = 0;
int32_t totalErrors = 0;
int32_t totalCalls = 0;

static void getSnpStats(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        //Now get the column
        Block_InstanceIterator *instanceIterator = block_getInstanceIterator(
                block);
        Segment *segment;
        char *referenceSeq = NULL;
        char *sampleSeq = NULL;
        goto end;
        while ((segment = block_getNext(instanceIterator)) != NULL) {
            if (strcmp(event_getHeader(segment_getEvent(segment)),
                    referenceEventString) == 0) {
                if (referenceSeq != NULL) {
                    goto end;
                }
                referenceSeq = segment_getString(segment);
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)),
                    sampleEventString) == 0) {
                if (sampleSeq != NULL) {
                    goto end;
                }
                sampleSeq = segment_getString(segment);
            }
        }

        if (referenceSeq != NULL) {
            //We're in gravy.
            for (int32_t i = ignoreFirstNBasesOfBlock; i < block_getLength(
                    block) - ignoreFirstNBasesOfBlock; i++) {
                totalSites++;
                if (sampleSeq != NULL) {
                    totalCorrect += bitsScoreFn(sampleSeq[i], referenceSeq[i]);
                    totalErrors += correctFn(sampleSeq[i], referenceSeq[i]) ? 0
                            : 1;
                    totalCalls++;
                }
            }
        }

        end:
        //cleanup
        if (referenceSeq != NULL) {
            free(referenceSeq);
        }
        if (sampleSeq != NULL) {
            free(sampleSeq);
        }
        block_destructInstanceIterator(instanceIterator);
    }
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "snpStats");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<substitutionStats>\n");
    EventTree_Iterator *eventIt = eventTree_getIterator(
            flower_getEventTree(flower));
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        sampleEventString = event_getHeader(event);
        if (sampleEventString != NULL && strcmp(sampleEventString, referenceEventString)
                != 0) {

            totalSites = getTotalLengthOfAdjacencies(flower, referenceEventString);
            totalCorrect = 0;
            totalErrors = 0;
            totalCalls = 0;

            getMAFsReferenceOrdered(flower, fileHandle, getSnpStats);

            ///////////////////////////////////////////////////////////////////////////
            // Print outputs
            ///////////////////////////////////////////////////////////////////////////

            fprintf(fileHandle, "<statsForSample "
                    "sampleName=\"%s\" "
                    "referenceName=\"%s\" "
                    "totalSites=\"%i\" "
                "totalCorrect=\"%f\" "
                "totalErrors=\"%i\" "
                "totalCalls=\"%i\" />\n",
                sampleEventString, referenceEventString,
                totalSites, totalCorrect, totalErrors, totalCalls);

        }
    }
    eventTree_destructIterator(eventIt);
    fprintf(fileHandle, "</substitutionStats>\n");

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    flower = NULL;
    cactusDisk_destruct(cactusDisk);

    //while(1);

    return 0;
}
