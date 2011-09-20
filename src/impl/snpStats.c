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
stList *events = NULL;

stHash *eventPairsToDistanceMatrix;

stInt64Tuple *getKey(Event *event1, Event *event2) {
    assert(event1 != NULL);
    assert(event2 != NULL);
    if (cactusMisc_nameCompare(event_getName(event1), event_getName(event2)) <= 0) {
        return stInt64Tuple_construct(2, event_getName(event1), event_getName(event2));
    }
    return stInt64Tuple_construct(2, event_getName(event2), event_getName(event1));
}

static void getSnpStatsDistanceMatrix(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        Block_InstanceIterator *instanceIterator1 = block_getInstanceIterator(
                block);
        Segment *segment1;
        while ((segment1 = block_getNext(instanceIterator1)) != NULL) {
            if (segment_getSequence(segment1) != NULL) {
                char *sequence1 = segment_getString(segment1);
                assert(sequence1 != NULL);
                Block_InstanceIterator *instanceIterator2 =
                        block_getInstanceIterator(block);
                Segment *segment2;
                while ((segment2 = block_getNext(instanceIterator2)) != NULL) {
                    if (segment_getSequence(segment2) != NULL && segment1 != segment2 && event_getName(segment_getEvent(segment1)) != event_getName(segment_getEvent(segment2))) {
                        char *sequence2 = segment_getString(segment2);
                        assert(sequence2 != NULL);
                        assert(strlen(sequence1) == strlen(sequence2));
                        assert(segment_getLength(segment1) == strlen(sequence1));
                        assert(segment_getLength(segment2) == strlen(sequence2));
                        for (int32_t i = ignoreFirstNBasesOfBlock; i < block_getLength(
                                            block) - ignoreFirstNBasesOfBlock; i++) {
                            stInt64Tuple *key = getKey(segment_getEvent(segment1), segment_getEvent(segment2));
                            int32_t *iA = stHash_search(eventPairsToDistanceMatrix, key);
                            if(iA == NULL) {
                                iA = st_calloc(2, sizeof(int32_t));
                                stHash_insert(eventPairsToDistanceMatrix, getKey(segment_getEvent(segment1), segment_getEvent(segment2)), iA);
                            }
                            if(!correctFn(sequence1[i], sequence2[i])) {
                                iA[0]++;
                            }
                            iA[1]++;
                            stInt64Tuple_destruct(key);
                        }
                        free(sequence2);
                    }
                }
                block_destructInstanceIterator(instanceIterator2);
                free(sequence1);
            }
        }
        block_destructInstanceIterator(instanceIterator1);
    }
}

static void getSnpStats(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        //Now get the column
        Block_InstanceIterator *instanceIterator = block_getInstanceIterator(
                block);
        Segment *segment;
        char *referenceSeq = NULL;
        Segment *referenceSegment = NULL;
        char *otherReferenceSeq = NULL;
        char *sampleSeq = NULL;
        Segment *sampleSegment = NULL;
        while ((segment = block_getNext(instanceIterator)) != NULL) {
            if (strcmp(event_getHeader(segment_getEvent(segment)),
                    referenceEventString) == 0) {
                if (referenceSeq != NULL) {
                    goto end;
                }
                referenceSeq = segment_getString(segment);
                referenceSegment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)),
                    sampleEventString) == 0) {
                if (sampleSeq != NULL) {
                    goto end;
                }
                sampleSeq = segment_getString(segment);
                sampleSegment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)),
                    otherReferenceEventString) == 0) {
                if (otherReferenceSeq != NULL) {
                    goto end;
                }
                otherReferenceSeq = segment_getString(segment);
            }
        }

        if (referenceSeq != NULL && otherReferenceSeq != NULL) {
            //We're in gravy.
            for (int32_t i = ignoreFirstNBasesOfBlock; i < block_getLength(
                    block) - ignoreFirstNBasesOfBlock; i++) {
                totalSites++;
                if (sampleSeq != NULL) {
                    totalCorrect += bitsScoreFn(sampleSeq[i], referenceSeq[i]);
                    bool b = correctFn(sampleSeq[i], referenceSeq[i]);
                    if (!b) {
                        assert(sampleSegment != NULL);
                        stList_append(events, sampleSegment);
                        assert(referenceSegment != NULL);
                        stList_append(events, referenceSegment);
                        stList_append(
                                events,
                                stIntTuple_construct(3, i,
                                        (int32_t) sampleSeq[i],
                                        (int32_t) referenceSeq[i]));
                    }
                    totalErrors += b ? 0 : 1;
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
        if (otherReferenceSeq != NULL) {
            free(otherReferenceSeq);
        }
        block_destructInstanceIterator(instanceIterator);
    }
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "snpStats");
    if (otherReferenceEventString == NULL) {
        otherReferenceEventString = stString_copy(referenceEventString);
    }

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
        if (sampleEventString != NULL && strcmp(sampleEventString,
                referenceEventString) != 0) {

            totalSites = getTotalLengthOfAdjacencies(flower,
                    referenceEventString);
            totalCorrect = 0;
            totalErrors = 0;
            totalCalls = 0;
            events = stList_construct();

            getMAFs(flower, fileHandle, getSnpStats);

            ///////////////////////////////////////////////////////////////////////////
            // Print outputs
            ///////////////////////////////////////////////////////////////////////////

            fprintf(fileHandle, "<statsForSample "
                "sampleName=\"%s\" "
                "referenceName=\"%s\" "
                "otherReferenceName=\"%s\" "
                "totalSites=\"%i\" "
                "totalCorrect=\"%f\" "
                "totalErrors=\"%i\" "
                "totalCalls=\"%i\">\n", sampleEventString, referenceEventString,
                    otherReferenceEventString, totalSites, totalCorrect,
                    totalErrors, totalCalls);

            for (int32_t i = 0; i < stList_length(events); i += 3) {
                Segment *sampleSegment = stList_get(events, i);
                assert(segment_getSequence(sampleSegment) != NULL);
                Segment *referenceSegment = stList_get(events, i + 1);
                assert(segment_getSequence(referenceSegment) != NULL);
                int32_t coordinate = stIntTuple_getPosition(
                        stList_get(events, i + 2), 0);
                char base1 = stIntTuple_getPosition(stList_get(events, i + 2),
                        1);
                char base2 = stIntTuple_getPosition(stList_get(events, i + 2),
                        2);
                fprintf(
                        fileHandle,
                        "%s %i %c %s %i %c\n",
                        sequence_getHeader(segment_getSequence(sampleSegment)),
                        segment_getStart(sampleSegment) + (segment_getStrand(
                                sampleSegment) ? coordinate : -coordinate)
                                - sequence_getStart(
                                        segment_getSequence(sampleSegment)),
                        base1,
                        sequence_getHeader(
                                segment_getSequence(referenceSegment)),
                        segment_getStart(referenceSegment)
                                + (segment_getStrand(referenceSegment) ? coordinate
                                        : -coordinate) - sequence_getStart(
                                segment_getSequence(referenceSegment)), base2);
                stIntTuple_destruct(stList_get(events, i + 2));
            }
            stList_destruct(events);
            fprintf(fileHandle, "</statsForSample>\n");
        }
    }
    eventTree_destructIterator(eventIt);

    //Now get snp distance matrix
    eventPairsToDistanceMatrix = stHash_construct3((uint32_t (*)(const void *))stInt64Tuple_hashKey, (int (*)(const void *, const void *))stInt64Tuple_equalsFn, (void (*)(void *))stInt64Tuple_destruct, free);
    getMAFs(flower, fileHandle, getSnpStatsDistanceMatrix);
    stHashIterator *hashIt = stHash_getIterator(eventPairsToDistanceMatrix);
    stInt64Tuple *eventPair;
    while((eventPair = stHash_getNext(hashIt)) != NULL) {
        Event *event1 = eventTree_getEvent(flower_getEventTree(flower), stInt64Tuple_getPosition(eventPair, 0));
        Event *event2 = eventTree_getEvent(flower_getEventTree(flower), stInt64Tuple_getPosition(eventPair, 1));
        assert(event1 != NULL);
        assert(event2 != NULL);
        assert(event1 != event2);
        assert(event_getHeader(event1) != NULL);
        assert(event_getHeader(event2) != NULL);
        int32_t *iA = stHash_search(eventPairsToDistanceMatrix, eventPair);
        fprintf(fileHandle, "<distancesForSamples eventName1=\"%s\" eventName2=\"%s\" substitutionNumber=\"%i\" sampleNumber=\"%i\" substitutionRate=\"%f\"/>\n", event_getHeader(event1), event_getHeader(event2), iA[0], iA[1], ((double)iA[0])/iA[1]);
    }
    stHash_destructIterator(hashIt);
    stHash_destruct(eventPairsToDistanceMatrix);

    fprintf(fileHandle, "</substitutionStats>\n");

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    flower = NULL;
    cactusDisk_destruct(cactusDisk);

    //while(1);

    return 0;
}
