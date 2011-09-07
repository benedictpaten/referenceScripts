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
int32_t *baseCoverages;
int32_t *totalBaseCoverages;
int32_t eventNumber;
int32_t totalReferenceBases;
int32_t totalOtherReferenceBases;

static void getMAFBlock2(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        stSortedSet *otherSampleEvents = stSortedSet_construct3(
                (int(*)(const void *, const void *)) strcmp, NULL);
        Segment *segment;
        Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
        int32_t sampleNumber = 0;
        while ((segment = block_getNext(instanceIt)) != NULL) {
            const char *segmentEvent = event_getHeader(
                    segment_getEvent(segment));
            if (strcmp(segmentEvent, sampleEventString) == 0) {
                sampleNumber++;
            } else if (strcmp(segmentEvent, referenceEventString) != 0) {
                stSortedSet_insert(otherSampleEvents, (void *) segmentEvent);
            }
        }
        block_destructInstanceIterator(instanceIt);
        baseCoverages[stSortedSet_size(otherSampleEvents)] += block_getLength(
                block) * sampleNumber;
        stSortedSet_destruct(otherSampleEvents);
        //Calculate bases in the reference and other reference sequence
        instanceIt = block_getInstanceIterator(block);
        bool includesReference = 0;
        bool includesOtherReference = 0;
        while ((segment = block_getNext(instanceIt)) != NULL) {
            const char *segmentEvent = event_getHeader(
                    segment_getEvent(segment));
            if (strcmp(segmentEvent, referenceEventString) == 0) {
                includesReference = 1;
            } else if (strcmp(segmentEvent, otherReferenceEventString) == 0) {
                includesOtherReference = 1;
            }
        }
        block_destructInstanceIterator(instanceIt);
        totalReferenceBases += includesReference ? block_getLength(
                block) * sampleNumber : 0;
        totalOtherReferenceBases += includesOtherReference ? block_getLength(
                        block) * sampleNumber : 0;
    }
}

void printStatsForSample(bool addToTotalBaseCoverage, FILE *fileHandle) {
    fprintf(fileHandle, "<statsForSample "
        "sampleName=\"%s\" "
        "referenceName=\"%s\" "
        "otherReferenceName=\"%s\" "
        "referenceBasesMapped=\"%i\" "
        "otherReferenceBasesMapped=\"%i\" "
        "baseCoverages=\"", sampleEventString, referenceEventString, otherReferenceEventString, totalReferenceBases, totalOtherReferenceBases);
    for (int32_t i = 0; i < eventNumber; i++) {
        fprintf(fileHandle, "%i ", baseCoverages[i]);
        if(addToTotalBaseCoverage) {
            totalBaseCoverages[i] += baseCoverages[i];
        }
    }
    fprintf(fileHandle, "\"/>\n");
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "coverageStats");
    assert(referenceEventString != NULL);
    assert(otherReferenceEventString != NULL);

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<coverageStats>\n");
    EventTree_Iterator *eventIt = eventTree_getIterator(
            flower_getEventTree(flower));
    eventNumber = eventTree_getEventNumber(
                        flower_getEventTree(flower));
    Event *event;
    totalBaseCoverages = st_calloc(sizeof(int32_t), eventNumber + 1);
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        sampleEventString = event_getHeader(event);
        if (sampleEventString != NULL && strcmp(sampleEventString,
                referenceEventString) != 0) {

            baseCoverages = st_calloc(sizeof(int32_t), eventNumber + 1);

            baseCoverages[0] = getTotalLengthOfAdjacencies(flower,
                    sampleEventString);

            getMAFs(flower, fileHandle, getMAFBlock2);

            printStatsForSample(1, fileHandle);

            free(baseCoverages);
        }
    }

    //Do aggregate base coverages..
    sampleEventString = "aggregate";
    baseCoverages = totalBaseCoverages;
    printStatsForSample(0, fileHandle);

    eventTree_destructIterator(eventIt);
    fprintf(fileHandle, "</coverageStats>\n");

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}

