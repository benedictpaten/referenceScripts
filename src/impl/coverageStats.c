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
int32_t referenceBases;
int32_t otherReferenceBases;
int32_t totalReferenceBases;
int32_t totalOtherReferenceBases;
bool ignoreOtherReferenceBlocks;

static void getMAFBlock2(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {

        //Calculate bases in the reference and other reference sequence
        Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
        bool includesReference = 0;
        bool includesOtherReference = 0;
        Segment *segment;
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
        if (ignoreOtherReferenceBlocks && includesOtherReference) {
            return;
        }

        stSortedSet *otherSampleEvents = stSortedSet_construct3(
                (int(*)(const void *, const void *)) strcmp, NULL);
        instanceIt = block_getInstanceIterator(block);
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

        referenceBases += includesReference ? block_getLength(block)
                * sampleNumber : 0;
        otherReferenceBases += includesOtherReference ? block_getLength(block)
                * sampleNumber : 0;
    }
}

void printStatsForSample(bool addToTotalBaseCoverage, FILE *fileHandle,
        int32_t denominator) {
    fprintf(fileHandle, "<statsForSample "
        "sampleName=\"%s\" "
        "referenceName=\"%s\" "
        "otherReferenceName=\"%s\" "
        "referenceBasesMapped=\"%i\" "
        "otherReferenceBasesMapped=\"%i\" "
        "baseCoverages=\"", sampleEventString, referenceEventString,
            otherReferenceEventString, (referenceBases / denominator),
            (otherReferenceBases / denominator));
    for (int32_t i = 0; i < eventNumber; i++) {
        fprintf(fileHandle, "%i ", (baseCoverages[i] / denominator));
        if (addToTotalBaseCoverage) {
            totalBaseCoverages[i] += baseCoverages[i];
        }
    }
    if (addToTotalBaseCoverage) {
        totalReferenceBases += referenceBases;
        totalOtherReferenceBases += otherReferenceBases;
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
    assert(outgroupEventString != NULL);

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    Sequence *referenceSequence = NULL;
    Sequence *otherReferenceSequence = NULL;
    Flower_SequenceIterator *sequenceIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        const char *eventHeader = event_getHeader(sequence_getEvent(sequence));
        if (eventHeader != NULL && strcmp(eventHeader, referenceEventString)
                == 0) {
            if (referenceSequence == NULL || sequence_getLength(sequence)
                    >= sequence_getLength(referenceSequence)) {
                referenceSequence = sequence;
            }
        }
        if (eventHeader != NULL && strcmp(eventHeader,
                otherReferenceEventString) == 0) {
            if (otherReferenceSequence == NULL || sequence_getLength(sequence)
                    >= sequence_getLength(otherReferenceSequence)) {
                otherReferenceSequence = sequence;
            }
        }
    }
    flower_destructSequenceIterator(sequenceIt);
    assert(referenceSequence != NULL);
    assert(otherReferenceSequence != NULL);

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(
            fileHandle,
            "<coverageStats referenceSequenceLength=\"%i\" otherReferenceSequenceLength=\"%i\">\n",
            sequence_getLength(referenceSequence),
            sequence_getLength(otherReferenceSequence));
    EventTree_Iterator *eventIt = eventTree_getIterator(
            flower_getEventTree(flower));
    eventNumber = eventTree_getEventNumber(flower_getEventTree(flower));
    Event * event;
    totalBaseCoverages = st_calloc(sizeof(int32_t), eventNumber + 1);
    totalReferenceBases = 0;
    totalOtherReferenceBases = 0;
    int32_t totalSamples = 0;
    ignoreOtherReferenceBlocks = 0;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        sampleEventString = event_getHeader(event);
        if (sampleEventString != NULL && strcmp(sampleEventString, "ROOT")
                != 0 && strcmp(sampleEventString, "") != 0) {

            baseCoverages = st_calloc(sizeof(int32_t), eventNumber + 1);

            baseCoverages[0] = getTotalLengthOfAdjacencies(flower,
                    sampleEventString);

            referenceBases = 0;
            otherReferenceBases = 0;

            getMAFs(flower, fileHandle, getMAFBlock2);

            printStatsForSample(
                    strcmp(sampleEventString, referenceEventString) != 0 && strcmp(sampleEventString, outgroupEventString) != 0,
                    fileHandle, 1);

            free(baseCoverages);

            totalSamples += (strcmp(sampleEventString, referenceEventString)
                    != 0 && strcmp(sampleEventString, outgroupEventString) != 0) ? 1 : 0;
        }
    }
    eventTree_destructIterator(eventIt);

    //Do average base coverages..
    sampleEventString = "average";
    baseCoverages = totalBaseCoverages;
    referenceBases = totalReferenceBases;
    otherReferenceBases = totalOtherReferenceBases;
    printStatsForSample(0, fileHandle, totalSamples);

    //Do blocks without other reference
    sampleEventString = referenceEventString;
    baseCoverages = st_calloc(sizeof(int32_t), eventNumber + 1);
    baseCoverages[0] = totalBaseCoverages[0];
    referenceBases = 0;
    getMAFs(flower, fileHandle, getMAFBlock2);
    otherReferenceBases = sequence_getLength(otherReferenceSequence);
    sampleEventString = "all";
    printStatsForSample(0, fileHandle, 1);
    free(baseCoverages);

    //Do blocks without other reference
    ignoreOtherReferenceBlocks = 1;
    sampleEventString = referenceEventString;
    baseCoverages = st_calloc(sizeof(int32_t), eventNumber + 1);
    baseCoverages[0] = totalBaseCoverages[0];
    referenceBases = 0;
    otherReferenceBases = 0;
    getMAFs(flower, fileHandle, getMAFBlock2);
    sampleEventString = "minusOtherReference";
    printStatsForSample(0, fileHandle, 1);
    free(baseCoverages);

    fprintf(fileHandle, "</coverageStats>\n");

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}

