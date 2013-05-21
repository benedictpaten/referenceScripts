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
int64_t totalSites = 0;
double totalCorrect = 0;
int64_t totalErrors = 0;
int64_t totalCalls = 0;
stList *events = NULL;

stHash *eventPairsToDistanceMatrix;

stIntTuple *getKey(Event *event1, Event *event2) {
    assert(event1 != NULL);
    assert(event2 != NULL);
    if (cactusMisc_nameCompare(event_getName(event1), event_getName(event2)) <= 0) {
        return stIntTuple_construct2( event_getName(event1), event_getName(event2));
    }
    return stIntTuple_construct2( event_getName(event2), event_getName(event1));
}

static void getSnpStatsDistanceMatrix(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        Block_InstanceIterator *instanceIterator1 = block_getInstanceIterator(block);
        Segment *segment1;
        while ((segment1 = block_getNext(instanceIterator1)) != NULL) {
            if (segment_getSequence(segment1) != NULL) {
                char *sequence1 = segment_getString(segment1);
                assert(sequence1 != NULL);
                Block_InstanceIterator *instanceIterator2 = block_getInstanceIterator(block);
                Segment *segment2;
                while ((segment2 = block_getNext(instanceIterator2)) != NULL) {
                    if (segment_getSequence(segment2) != NULL && segment1 != segment2 && event_getName(
                            segment_getEvent(segment1)) != event_getName(segment_getEvent(segment2))) {
                        char *sequence2 = segment_getString(segment2);
                        assert(sequence2 != NULL);
                        assert(strlen(sequence1) == strlen(sequence2));
                        assert(segment_getLength(segment1) == strlen(sequence1));
                        assert(segment_getLength(segment2) == strlen(sequence2));
                        for (int64_t i = ignoreFirstNBasesOfBlock; i < block_getLength(block)
                                - ignoreFirstNBasesOfBlock; i++) {
                            stIntTuple *key = getKey(segment_getEvent(segment1), segment_getEvent(segment2));
                            int64_t *iA = stHash_search(eventPairsToDistanceMatrix, key);
                            if (iA == NULL) {
                                iA = st_calloc(2, sizeof(int64_t));
                                stHash_insert(eventPairsToDistanceMatrix,
                                        getKey(segment_getEvent(segment1), segment_getEvent(segment2)), iA);
                            }
                            if (toupper(sequence1[i]) != 'N' && toupper(sequence2[i]) != 'N' && toupper(sequence1[i])
                                    != toupper(sequence2[i])) {
                                iA[0]++;
                            }
                            iA[1]++;
                            stIntTuple_destruct(key);
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
        Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
        Segment *segment;
        Segment *referenceSegment = NULL;
        Segment *sampleSegment = NULL;
        Segment *otherReferenceSegment = NULL;
        stList *otherSegments = stList_construct();
        while ((segment = block_getNext(instanceIterator)) != NULL) {
            if (segment_getSequence(segment) == NULL) {
                continue;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)), referenceEventString) == 0) {
                if (referenceSegment != NULL) {
                    goto end;
                }
                referenceSegment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)), sampleEventString) == 0) {
                if (sampleSegment != NULL) {
                    goto end;
                }
                sampleSegment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)), otherReferenceEventString) == 0) {
                if (otherReferenceSegment != NULL || ignoreSitesWithOtherReferencePresent) {
                    goto end;
                }
                otherReferenceSegment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)), referenceEventString) != 0 && strcmp(
                    event_getHeader(segment_getEvent(segment)), sampleEventString) != 0 && strcmp(
                    event_getHeader(segment_getEvent(segment)), otherReferenceEventString) != 0) {
                stList_append(otherSegments, segment);
            }
        }

        if (referenceSegment != NULL && otherReferenceSegment != NULL && sampleSegment != NULL) {
            bool b = segment_getStrand(referenceSegment);
            referenceSegment = b ? referenceSegment : segment_getReverse(referenceSegment);
            otherReferenceSegment = b ? otherReferenceSegment : segment_getReverse(otherReferenceSegment);
            sampleSegment = b ? sampleSegment : segment_getReverse(sampleSegment);
            stList_setDestructor(otherSegments, free);
            for (int64_t i = 0; i < stList_length(otherSegments); i++) {
                stList_set(
                        otherSegments,
                        i,
                        segment_getString(
                                b ? stList_get(otherSegments, i) : segment_getReverse(stList_get(otherSegments, i))));
            }
            char *referenceSeq = segment_getString(referenceSegment);
            char *otherReferenceSeq = segment_getString(otherReferenceSegment);
            char *sampleSeq = segment_getString(sampleSegment);
            //We're in gravy.
            for (int64_t i = ignoreFirstNBasesOfBlock; i < block_getLength(block) - ignoreFirstNBasesOfBlock; i++) {
                totalSites++;
                if (sampleSeq != NULL) {
                    totalCorrect += bitsScoreFn(sampleSeq[i], referenceSeq[i]);
                    bool b = correctFn(sampleSeq[i], referenceSeq[i]);
                    if (!b) {
                        assert(sampleSegment != NULL);
                        stList_append(events, sampleSegment);
                        assert(referenceSegment != NULL);
                        stList_append(events, referenceSegment);
                        int64_t k = 0;
                        for (int64_t j = 0; j < stList_length(otherSegments); j++) {
                            char *seq = stList_get(otherSegments, j);
                            if (bitsScoreFn(sampleSeq[i], seq[i]) != 0) {
                                k++;
                            }
                        }
                        stList_append(events,
                                stIntTuple_construct4( i, (int64_t) sampleSeq[i], (int64_t) referenceSeq[i], k));
                    }
                    totalErrors += b ? 0 : 1;
                    totalCalls++;
                }
            }
            free(referenceSeq);
            free(otherReferenceSeq);
            free(sampleSeq);
        }
        end: stList_destruct(otherSegments);
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
    EventTree_Iterator *eventIt = eventTree_getIterator(flower_getEventTree(flower));
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        sampleEventString = event_getHeader(event);
        if (sampleEventString != NULL && strcmp(sampleEventString, referenceEventString) != 0) {

            totalSites = getTotalLengthOfAdjacencies(flower, referenceEventString);
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
                "totalSites=\"%" PRIi64 "\" "
                "totalCorrect=\"%f\" "
                "totalErrors=\"%" PRIi64 "\" "
                "totalCalls=\"%" PRIi64 "\">\n", sampleEventString, referenceEventString, otherReferenceEventString, totalSites,
                    totalCorrect, totalErrors, totalCalls);

            for (int64_t i = 0; i < stList_length(events); i += 3) {
                Segment *sampleSegment = stList_get(events, i);
                assert(segment_getSequence(sampleSegment) != NULL);
                Segment *referenceSegment = stList_get(events, i + 1);
                assert(segment_getSequence(referenceSegment) != NULL);
                int64_t coordinate = stIntTuple_get(stList_get(events, i + 2), 0);
                char base1 = stIntTuple_get(stList_get(events, i + 2), 1);
                char base2 = stIntTuple_get(stList_get(events, i + 2), 2);
                int64_t recurrence = stIntTuple_get(stList_get(events, i + 2), 3);
                if (recurrence + 1 >= minimumRecurrence) {
                    fprintf(
                            fileHandle,
                            "%s %" PRIi64 " %c %s %" PRIi64 " %c\n",
                            sequence_getHeader(segment_getSequence(sampleSegment)),
                            segment_getStart(sampleSegment) + (segment_getStrand(sampleSegment) ? coordinate
                                    : -coordinate) - sequence_getStart(segment_getSequence(sampleSegment)),
                            base1,
                            sequence_getHeader(segment_getSequence(referenceSegment)),
                            segment_getStart(referenceSegment) + (segment_getStrand(referenceSegment) ? coordinate
                                    : -coordinate) - sequence_getStart(segment_getSequence(referenceSegment)), base2);
                }
                stIntTuple_destruct(stList_get(events, i + 2));
            }
            stList_destruct(events);
            fprintf(fileHandle, "</statsForSample>\n");
        }
    }
    eventTree_destructIterator(eventIt);

    //Now get snp distance matrix
    eventPairsToDistanceMatrix = stHash_construct3((uint64_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, (void(*)(void *)) stIntTuple_destruct, free);
    getMAFs(flower, fileHandle, getSnpStatsDistanceMatrix);
    stHashIterator *hashIt = stHash_getIterator(eventPairsToDistanceMatrix);
    stIntTuple *eventPair;
    while ((eventPair = stHash_getNext(hashIt)) != NULL) {
        Event *event1 = eventTree_getEvent(flower_getEventTree(flower), stIntTuple_get(eventPair, 0));
        Event *event2 = eventTree_getEvent(flower_getEventTree(flower), stIntTuple_get(eventPair, 1));
        assert(event1 != NULL);
        assert(event2 != NULL);
        assert(event1 != event2);
        assert(event_getHeader(event1) != NULL);
        assert(event_getHeader(event2) != NULL);
        int64_t *iA = stHash_search(eventPairsToDistanceMatrix, eventPair);
        fprintf(
                fileHandle,
                "<distancesForSamples eventName1=\"%s\" eventName2=\"%s\" substitutionNumber=\"%" PRIi64 "\" sampleNumber=\"%" PRIi64 "\" substitutionRate=\"%f\"/>\n",
                event_getHeader(event1), event_getHeader(event2), iA[0], iA[1], ((double) iA[0]) / iA[1]);
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
