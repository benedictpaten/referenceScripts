/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus.h"
#include "cactusMafs.h"
#include "contigPaths.h"
#include "adjacencyTraversal.h"
#include "adjacencyClassification.h"
#include "scaffoldPaths.h"
#include "referenceCommon.h"

static stHash *setOfPairs;
const char *sampleEventString;

void addToPairs(int64_t referenceNumber, int64_t sampleNumber, int64_t totalLength) {
    stIntTuple *key = stIntTuple_construct2( referenceNumber, sampleNumber);
    int64_t *value;
    if ((value = stHash_search(setOfPairs, key)) == NULL) {
        value = st_malloc(sizeof(int64_t));
        value[0] = 0;
        stHash_insert(setOfPairs, key, value);
    } else {
        stIntTuple_destruct(key);
    }
    value[0] += totalLength;
}

static void getMAFBlock2(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        Segment *segment;
        Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
        int64_t referenceNumber = 0, sampleNumber = 0;
        while ((segment = block_getNext(instanceIt)) != NULL) {
            const char *segmentEvent = event_getHeader(segment_getEvent(segment));
            if (strcmp(segmentEvent, referenceEventString) == 0) { //Establish if we need a line..
                referenceNumber++;
            } else if (strcmp(segmentEvent, sampleEventString) == 0) {
                sampleNumber++;
            }
        }
        block_destructInstanceIterator(instanceIt);
        if (referenceNumber > 0 || sampleNumber > 0) {
            addToPairs(referenceNumber, sampleNumber, block_getLength(block));
        }
    }
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "copyNumberStats");

    ///////////////////////////////////////////////////////////////////////////
    // Now use the MAF printing code to generate the results..
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<copyNumberStats>");
    EventTree_Iterator *eventIt = eventTree_getIterator(flower_getEventTree(flower));
    Event *event;

    while ((event = eventTree_getNext(eventIt)) != NULL) {
        sampleEventString = event_getHeader(event);
        if (sampleEventString != NULL && strcmp(sampleEventString, referenceEventString) != 0) {
            //The pairs to represent the mafs.
            setOfPairs = stHash_construct3((uint64_t(*)(const void *)) stIntTuple_hashKey,
                    (int(*)(const void *, const void *)) stIntTuple_equalsFn, (void(*)(void *)) stIntTuple_destruct,
                    free);
            //Add in the score from the adjacencies
            addToPairs(0, 1, getTotalLengthOfAdjacencies(flower, sampleEventString));

            //Pass over the blocks.
            getMAFs(flower, fileHandle, getMAFBlock2);
            //Now calculate the linkage stats
            stList *copyNumbers = stHash_getKeys(setOfPairs);
            stList_sort(copyNumbers, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
            int64_t totalCopyNumberDeficientColumns = 0;
            int64_t totalCopyNumberDeficientBases = 0;
            int64_t totalCopyNumberDeficientColumnsGreaterThanZero = 0;
            int64_t totalCopyNumberDeficientBasesGreaterThanZero = 0;
            int64_t totalCopyNumberExcessColumns = 0;
            int64_t totalCopyNumberExcessBases = 0;
            int64_t totalColumnCount = 0;
            int64_t totalBaseCount = 0;
            for (int64_t i = 0; i < stList_length(copyNumbers); i++) {
                stIntTuple *copyNumber = stList_get(copyNumbers, i);
                int64_t *columnCount = stHash_search(setOfPairs, copyNumber);
                totalColumnCount += columnCount[0];
                totalBaseCount += stIntTuple_getPosition(copyNumber, 1) * columnCount[0];
            }
            fprintf(
                    fileHandle,
                    "<statsForSample sampleName=\"%s\" referenceName=\"%s\" minimumBlockLength=\"%" PRIi64 "\" totalColumnCount=\"%" PRIi64 "\" totalBaseCount=\"%" PRIi64 "\">\n",
                    sampleEventString, referenceEventString, minimumBlockLength, totalColumnCount, totalBaseCount);
            for (int64_t i = 0; i < stList_length(copyNumbers); i++) {
                stIntTuple *copyNumber = stList_get(copyNumbers, i);
                int64_t *columnCount = stHash_search(setOfPairs, copyNumber);
                int64_t referenceNumber = stIntTuple_getPosition(copyNumber, 0);
                int64_t assemblyNumber = stIntTuple_getPosition(copyNumber, 1);
                assert(assemblyNumber >= 0);
                assert(columnCount != NULL);
                assert(columnCount[0] >= 0);
                fprintf(
                        fileHandle,
                        "<copyNumberCategory referenceCopyNumber=\"%" PRIi64 "\" assemblyCopyNumber=\"%" PRIi64 "\" columnCount=\"%" PRIi64 "\"/>\n",
                        referenceNumber, assemblyNumber, columnCount[0]);
                if (assemblyNumber < referenceNumber) {
                    totalCopyNumberDeficientColumns += columnCount[0];
                    totalCopyNumberDeficientBases += columnCount[0] * (referenceNumber - assemblyNumber);
                    if (assemblyNumber > 0) {
                        totalCopyNumberDeficientColumnsGreaterThanZero += columnCount[0];
                        totalCopyNumberDeficientBasesGreaterThanZero += columnCount[0] * (referenceNumber
                                - assemblyNumber);
                    }
                } else if (assemblyNumber > referenceNumber) {
                    totalCopyNumberExcessColumns += columnCount[0];
                    totalCopyNumberExcessBases += columnCount[0] * (assemblyNumber - referenceNumber);
                }
            }
            fprintf(
                    fileHandle,
                    "<deficientCopyNumberCounts totalColumns=\"%" PRIi64 "\" totalBases=\"%" PRIi64 "\" totalProportionOfColumns=\"%f\" totalProportionOfBases=\"%f\"/>",
                    totalCopyNumberDeficientColumns, totalCopyNumberDeficientBases,
                    ((float) totalCopyNumberDeficientColumns) / totalColumnCount,
                    ((float) totalCopyNumberDeficientBases) / totalBaseCount);
            fprintf(
                    fileHandle,
                    "<deficientCopyNumberCountsGreaterThanZero totalColumns=\"%" PRIi64 "\" totalBases=\"%" PRIi64 "\" totalProportionOfColumns=\"%f\" totalProportionOfBases=\"%f\"/>",
                    totalCopyNumberDeficientColumnsGreaterThanZero, totalCopyNumberDeficientBasesGreaterThanZero,
                    ((float) totalCopyNumberDeficientColumnsGreaterThanZero) / totalColumnCount,
                    ((float) totalCopyNumberDeficientBasesGreaterThanZero) / totalBaseCount);
            fprintf(
                    fileHandle,
                    "<excessCopyNumberCounts totalColumns=\"%" PRIi64 "\" totalBases=\"%" PRIi64 "\" totalProportionOfColumns=\"%f\" totalProportionOfBases=\"%f\"/>",
                    totalCopyNumberExcessColumns, totalCopyNumberExcessBases,
                    ((float) totalCopyNumberExcessColumns) / totalColumnCount,
                    ((float) totalCopyNumberExcessBases) / totalBaseCount);
            fprintf(fileHandle, "</statsForSample>\n");
        }
    }
    fprintf(fileHandle, "</copyNumberStats>\n");
    fclose(fileHandle);

    st_logInfo("Got the copy number counts in %" PRIi64 " seconds/\n", time(NULL) - startTime);

    return 0;
}
