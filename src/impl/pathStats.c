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

typedef struct _indelEvent {
    Cap *cap, *otherCap;
    int32_t insertionLength;
    int32_t deletionLength;
} IndelEvent;

IndelEvent *indelEvent_construct(Cap *cap, Cap *otherCap,
        int32_t insertionLength, int32_t deletionLength) {
    IndelEvent *indelEvent = st_malloc(sizeof(IndelEvent));
    cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
    otherCap = cap_getStrand(otherCap) ? otherCap : cap_getReverse(otherCap);
    assert(cap_getSide(cap) != cap_getSide(otherCap));
    if (cap_getSide(cap)) {
        Cap *cap2 = otherCap;
        otherCap = cap;
        cap = cap2;
    }
    assert(cap_getCoordinate(cap) < cap_getCoordinate(otherCap));
    indelEvent->cap = cap;
    indelEvent->otherCap = otherCap;
    indelEvent->insertionLength = insertionLength;
    indelEvent->deletionLength = deletionLength;
    return indelEvent;
}

int indelEvent_cmpFn(IndelEvent *indelEvent1, IndelEvent *indelEvent2) {
    int i = cactusMisc_nameCompare(cap_getName(indelEvent1->cap),
            cap_getName(indelEvent2->cap));
    if (i != 0) {
        return i;
    }
    assert(indelEvent1->cap == indelEvent2->cap);
    assert(indelEvent1->otherCap == indelEvent2->otherCap);
    assert(indelEvent1->insertionLength == indelEvent2->insertionLength);
    assert(indelEvent1->deletionLength == indelEvent2->deletionLength);
    return 0;
}

typedef struct _sampleStats {
    int32_t totalCleanEnds;
    int32_t totalHangingEndWithsNs;
    int32_t totalScaffoldGaps;
    int32_t totalAmbiguityGaps;

    int32_t totalIntraJoins;
    int32_t totalInterJoins;
    int32_t totalInsertions;
    int32_t totalDeletions;
    int32_t totalInsertionAndDeletions;
    int32_t totalHangingInsertions;

    stList *indelEvents;

    stList *contigPathLengthDistribution;
    stList *scaffoldPathLengthDistribution;
    stList *blockLengthDistribution;
    stList *sampleSequenceLengthDistribution;
} SampleStats;

SampleStats *sampleStats_construct() {
    SampleStats *sampleStats = st_calloc(1, sizeof(SampleStats));
    sampleStats->indelEvents = stList_construct3(0, free);
    sampleStats->contigPathLengthDistribution = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    sampleStats->scaffoldPathLengthDistribution = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    sampleStats->blockLengthDistribution = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    sampleStats->sampleSequenceLengthDistribution = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    return sampleStats;
}

void sampleStats_destruct(SampleStats *sampleStats) {
    stList_destruct(sampleStats->indelEvents);
    stList_destruct(sampleStats->contigPathLengthDistribution);
    stList_destruct(sampleStats->scaffoldPathLengthDistribution);
    stList_destruct(sampleStats->blockLengthDistribution);
    stList_destruct(sampleStats->sampleSequenceLengthDistribution);
    free(sampleStats);
}

void getHaplotypePathStatsP(Cap *cap, stList *referenceEventStrings,
        stList *contaminationEventStrings,
        CapCodeParameters *capCodeParameters, SampleStats *sampleStats) {
    int32_t insertLength, deleteLength;
    Cap *otherCap = NULL;
    switch (getCapCode(cap, &otherCap, referenceEventStrings,
            contaminationEventStrings, &insertLength, &deleteLength,
            capCodeParameters)) {
        case HAP_NOTHING:
            return;
        case CONTIG_END:
            sampleStats->totalCleanEnds++;
            return;
        case CONTIG_END_WITH_SCAFFOLD_GAP:
        case CONTIG_END_WITH_AMBIGUITY_GAP:
            sampleStats->totalHangingEndWithsNs++;
            return;
        case SCAFFOLD_GAP:
        case AMBIGUITY_GAP:
            if (insertLength > 0 && deleteLength == 0) {
                stList_append(sampleStats->indelEvents,
                        indelEvent_construct(cap, otherCap, insertLength, 0));
                sampleStats->totalInsertions++;
            }
            if (deleteLength > 0 && insertLength == 0) {
                stList_append(sampleStats->indelEvents,
                        indelEvent_construct(cap, otherCap, 0, deleteLength));
                sampleStats->totalDeletions++;
            }
            if (insertLength > 0 && deleteLength > 0) {
                stList_append(
                        sampleStats->indelEvents,
                        indelEvent_construct(cap, otherCap, insertLength,
                                deleteLength));
                sampleStats->totalInsertions++;
                sampleStats->totalDeletions++;
            }
            sampleStats->totalScaffoldGaps++;
            return;
        case ERROR_HAP_TO_HAP_SAME_CHROMOSOME:
            sampleStats->totalIntraJoins++;
            return;
        case ERROR_HAP_TO_HAP_DIFFERENT_CHROMOSOMES:
            sampleStats->totalInterJoins++;
            return;
        case ERROR_HAP_TO_INSERT:
            assert(insertLength > 0);
            stList_append(sampleStats->indelEvents,
                    indelEvent_construct(cap, otherCap, insertLength, 0));
            sampleStats->totalInsertions++;
            return;
        case ERROR_HAP_TO_DELETION:
            assert(deleteLength > 0);
            stList_append(sampleStats->indelEvents,
                    indelEvent_construct(cap, otherCap, 0, deleteLength));
            sampleStats->totalDeletions++;
            return;
        case ERROR_HAP_TO_INSERT_AND_DELETION:
            assert(insertLength > 0);
            assert(deleteLength > 0);
            stList_append(
                    sampleStats->indelEvents,
                    indelEvent_construct(cap, otherCap, insertLength,
                            deleteLength));
            sampleStats->totalInsertionAndDeletions++;
            return;
        case ERROR_CONTIG_END_WITH_INSERT:
            sampleStats->totalHangingInsertions++;
            return;
            //These should not occur
        case HAP_SWITCH:
            assert(0);
        case ERROR_HAP_TO_CONTAMINATION:
            assert(0);
        case ERROR_HAP_TO_INSERT_TO_CONTAMINATION:
            assert(0);
    }
}

void printIndel(IndelEvent *indelEvent, stList *eventStrings,
        FILE *fileHandle) {
    Cap *cap3 = NULL, *cap4 = NULL;
    int32_t i = 0;
    bool b = endsAreAdjacent2(cap_getEnd(indelEvent->cap),
            cap_getEnd(indelEvent->otherCap), &cap3, &cap4, &i, eventStrings);
    cap3 = cap_getStrand(cap3) ? cap3 : cap_getReverse(cap3);
    cap4 = cap_getStrand(cap4) ? cap4 : cap_getReverse(cap4);
    if (cap_getSide(cap3)) {
        assert(cap_getCoordinate(cap3) - cap_getCoordinate(cap4) > 0);
    } else {
        assert(cap_getCoordinate(cap4) - cap_getCoordinate(cap3) > 0);
    }
    assert(b);
    assert(i == indelEvent->deletionLength);
    assert(
            cap_getCoordinate(indelEvent->otherCap) - cap_getCoordinate(
                    indelEvent->cap) - 1 == indelEvent->insertionLength);
    assert(cap_getSequence(indelEvent->cap) != NULL);
    assert(cap_getSequence(cap3) != NULL);

    fprintf(fileHandle,
            "SEQ1: %s START1: %i LENGTH1: %i STRAND1: 1 SEQ2: %s START2: %i LENGTH2: %i STRAND2: %i\n",
            sequence_getHeader(cap_getSequence(indelEvent->cap)),
            cap_getCoordinate(indelEvent->cap) + 1 - sequence_getStart(
                    cap_getSequence(indelEvent->cap)),
            indelEvent->insertionLength,
            sequence_getHeader(cap_getSequence(cap3)),
            cap_getSide(cap3) ? cap_getCoordinate(cap3) - 1
                    - sequence_getStart(cap_getSequence(cap3))
                    : cap_getCoordinate(cap3) + 1 - sequence_getStart(
                            cap_getSequence(cap3)), indelEvent->deletionLength,
            !cap_getSide(cap3));
}

void getSequenceLengths(Flower *flower, const char *eventString,
        stList *sequenceLengths) {
    Flower_SegmentIterator *sequenceIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        if (strcmp(event_getHeader(sequence_getEvent(sequence)), eventString)
                == 0) {
            stList_append(sequenceLengths,
                    stIntTuple_construct(1, sequence_getLength(sequence)));
        }
    }
    flower_destructSequenceIterator(sequenceIt);
}

void getBlockLengths(Flower *flower, const char *sampleEventString,
        const char *referenceEventString, stList *blockLengths) {
    //Find blocks..
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIt)) != NULL) {
        if (hasCapInEvent(block_get5End(block), sampleEventString)
                && hasCapInEvent(block_get5End(block), referenceEventString)) {
            stList_append(blockLengths,
                    stIntTuple_construct(1, block_getLength(block)));
        }
    }
    flower_destructBlockIterator(blockIt);

    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (!group_isLeaf(group)) {
            getBlockLengths(group_getNestedFlower(group), sampleEventString,
                    referenceEventString, blockLengths);
        }
    }
    flower_destructGroupIterator(groupIt);
}

void getScaffoldPathsLengths(stHash *scaffoldPaths, stList *scaffoldPathLengths) {
    stList *scaffoldPathsList = stHash_getValues(scaffoldPaths);
    for (int32_t i = 0; i < stList_length(scaffoldPathsList); i++) {
        stSortedSet *scaffoldPath = stList_get(scaffoldPathsList, i);
        stSortedSetIterator *contigPathIt = stSortedSet_getIterator(
                scaffoldPath);
        stList *contigPath;
        int32_t j = 0;
        while ((contigPath = stSortedSet_getNext(contigPathIt)) != NULL) {
            j += contigPathLength(contigPath);
        }
        stSortedSet_destructIterator(contigPathIt);
        stList_append(scaffoldPathLengths, stIntTuple_construct(1, j));
    }
    stList_destruct(scaffoldPathsList);
}

SampleStats *getSamplePathStats(Flower *flower, const char *sampleEventString,
        const char *referenceEventString, CapCodeParameters *capCodeParameters) {
    /*
     * Gets stats for a given sample.
     */

    SampleStats *sampleStats = sampleStats_construct();

    stList *referenceEventStrings = stList_construct(); //This is the holder of the event strings
    stList_append(referenceEventStrings, (char *) referenceEventString);
    stList *emptyList = stList_construct();

    //Calculate the haplotype paths.
    stList *contigPaths = getContigPaths(flower, sampleEventString,
            referenceEventStrings);
    stHash *scaffoldPaths = getScaffoldPaths(contigPaths,
            referenceEventStrings, emptyList, capCodeParameters);

    //Now calculate the stats for the ends of the contig paths.
    for (int32_t i = 0; i < stList_length(contigPaths); i++) {
        stList *contigPath = stList_get(contigPaths, i);
        stList_append(sampleStats->contigPathLengthDistribution,
                stIntTuple_construct(1, contigPathLength(contigPath)));
        for (int32_t j = 0; j < stList_length(contigPath); j++) {
            Segment *segment = stList_get(contigPath, j);
            getHaplotypePathStatsP(segment_get5Cap(segment),
                    referenceEventStrings, emptyList, capCodeParameters,
                    sampleStats);
            getHaplotypePathStatsP(segment_get3Cap(segment),
                    referenceEventStrings, emptyList, capCodeParameters,
                    sampleStats);
        }
    }
    //Normalise stats
    assert(sampleStats->totalScaffoldGaps % 2 == 0);
    assert(sampleStats->totalAmbiguityGaps == 0);
    assert(sampleStats->totalInsertions % 2 == 0);
    assert(sampleStats->totalDeletions % 2 == 0);
    assert(sampleStats->totalInsertionAndDeletions % 2 == 0);
    assert(sampleStats->totalIntraJoins % 2 == 0);
    assert(sampleStats->totalInterJoins % 2 == 0);

    sampleStats->totalScaffoldGaps /= 2;
    sampleStats->totalAmbiguityGaps /= 2;
    sampleStats->totalInsertions /= 2;
    sampleStats->totalDeletions /= 2;
    sampleStats->totalInsertionAndDeletions /= 2;
    sampleStats->totalIntraJoins /= 2;
    sampleStats->totalInterJoins /= 2;
    sampleStats->totalInsertions += sampleStats->totalInsertionAndDeletions;
    sampleStats->totalDeletions += sampleStats->totalInsertionAndDeletions;

    //Add the lengths for the scaffold paths.
    getScaffoldPathsLengths(scaffoldPaths,
            sampleStats->scaffoldPathLengthDistribution);
    getBlockLengths(flower, sampleEventString, referenceEventString,
            sampleStats->blockLengthDistribution);
    getSequenceLengths(flower, sampleEventString,
            sampleStats->sampleSequenceLengthDistribution);

    //Cleanup
    stHash_destruct(scaffoldPaths);
    stList_destruct(contigPaths);
    stList_destruct(referenceEventStrings);
    stList_destruct(emptyList);

    return sampleStats;
}

int32_t getN50(stList *lengths, int32_t genomeLength) {
    int32_t totalLength = 0;
    int32_t pJ = INT32_MAX;
    stList_sort(lengths, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    for (int32_t i = stList_length(lengths) - 1; i >= 0; i--) {
        int32_t j = stIntTuple_getPosition(stList_get(lengths, i), 0);
        assert(j <= pJ);
        pJ = j;
        totalLength += j;
        if (totalLength >= genomeLength / 2) {
            return j;
        }
    }
    return 0; //this should not happen!
}

int32_t getSum(stList *lengths) {
    int32_t j = 0;
    for (int32_t i = 0; i < stList_length(lengths); i++) {
        j += stIntTuple_getPosition(stList_get(lengths, i), 0);
    }
    return j;
}

int32_t getGenomeLength(Flower *flower, const char *eventString) {
    stList *sequences = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    getSequenceLengths(flower, eventString, sequences);
    int32_t i = getSum(sequences);
    stList_destruct(sequences);
    return i;
}

void reportPathStatsForReference(Flower *flower, FILE *fileHandle,
        const char *referenceEventString, CapCodeParameters *capCodeParameters) {

    int32_t referenceGenomeLength = getGenomeLength(flower,
            referenceEventString);

    EventTree_Iterator *eventIt = eventTree_getIterator(
            flower_getEventTree(flower));
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        const char *eventString = event_getHeader(event);
        if (eventString != NULL && strcmp(eventString, referenceEventString)
                != 0) {
            stList *eventStrings = stList_construct();
            stList_append(eventStrings, (void *)eventString);
            SampleStats *sampleStats = getSamplePathStats(flower,
                    referenceEventString, eventString, capCodeParameters);
            int32_t minimumNCount = capCodeParameters->minimumNCount;
            capCodeParameters->minimumNCount = 0;
            SampleStats *sampleStatsFreeIndels = getSamplePathStats(flower,
                    referenceEventString, eventString, capCodeParameters);
            capCodeParameters->minimumNCount = minimumNCount;

            int32_t genomeLength = getGenomeLength(flower, eventString);
            int32_t alignedGenomeLength = getSum(
                    sampleStats->blockLengthDistribution);

            stSortedSet *indelEvents = stList_getSortedSet(
                    sampleStats->indelEvents,
                    (int(*)(const void *, const void *)) indelEvent_cmpFn);
            assert(
                    stList_length(sampleStats->indelEvents) == 2
                            * stSortedSet_size(indelEvents));

            int32_t totalEvents = sampleStats->totalIntraJoins
                    + sampleStats->totalInterJoins
                    + sampleStats->totalInsertions
                    + sampleStats->totalDeletions
                    + sampleStats->totalInsertionAndDeletions
                    + sampleStats->totalHangingInsertions;

            fprintf(
                    fileHandle,
                    "<statsForSample "
                        "sampleName=\"%s\" "
                        "referenceName=\"%s\" "
                        "totalScaffoldGaps=\"%i\" "
                        "totalContigEnds=\"%i\" "
                        "totalContigEndsWithNs=\"%i\" "
                        "totalIntraJoin=\"%i\" "
                        "totalInterJoin=\"%i\" "
                        "totalInsertion=\"%i\" "
                        "totalDeletion=\"%i\" "
                        "totalInsertionAndDeletion=\"%i\" "
                        "totalHangingInsertion=\"%i\" "
                        "totalDifferences=\"%i\" "

                        "totalScaffoldGapsPerAlignedBase=\"%f\" "
                        "totalContigEndsPerAlignedBase=\"%f\" "
                        "totalContigEndsWithNsPerAlignedBase=\"%f\" "
                        "totalIntraJoinPerAlignedBase=\"%f\" "
                        "totalInterJoinPerAlignedBase=\"%f\" "
                        "totalInsertionPerAlignedBase=\"%f\" "
                        "totalDeletionPerAlignedBase=\"%f\" "
                        "totalInsertionAndDeletionPerAlignedBase=\"%f\" "
                        "totalHangingInsertionPerAlignedBase=\"%f\" "
                        "totalDifferencesPerAlignedBase=\"%f\" "

                        "totalSampleGenomeLength=\"%i\" "
                        "totalReferenceGenomeLength=\"%i\" "
                        "totalAlignedGenomeLength=\"%i\" "
                        "sequenceN50=\"%i\" "
                        "blockN50=\"%i\" "
                        "contigPathN50=\"%i\" "
                        "scaffoldPathN50=\"%i\" "
                        "contigPathN50Check=\"%i\" "
                        "indelPathN50=\"%i\" "
                        "totalSequenceNumber=\"%i\" "
                        "totalBlockNumber=\"%i\" "
                        "totalContigPaths=\"%i\" "
                        "totalScaffoldPaths=\"%i\" "
                        "totalIndelPaths=\"%i\" ",
                    eventString,
                    referenceEventString,
                    sampleStats->totalScaffoldGaps
                            + sampleStats->totalAmbiguityGaps,
                    sampleStats->totalCleanEnds,
                    sampleStats->totalHangingEndWithsNs,
                    sampleStats->totalIntraJoins,
                    sampleStats->totalInterJoins,
                    sampleStats->totalInsertions,
                    sampleStats->totalDeletions,
                    sampleStats->totalInsertionAndDeletions,
                    sampleStats->totalHangingInsertions,
                    totalEvents,

                    (double) (sampleStats->totalScaffoldGaps
                            + sampleStats->totalAmbiguityGaps)
                            / alignedGenomeLength,
                    (double) sampleStats->totalCleanEnds / alignedGenomeLength,
                    (double) sampleStats->totalHangingEndWithsNs
                            / alignedGenomeLength,
                    (double) sampleStats->totalIntraJoins / alignedGenomeLength,
                    (double) sampleStats->totalInterJoins / alignedGenomeLength,
                    (double) sampleStats->totalInsertions / alignedGenomeLength,
                    (double) sampleStats->totalDeletions / alignedGenomeLength,
                    (double) sampleStats->totalInsertionAndDeletions
                            / alignedGenomeLength,
                    (double) sampleStats->totalHangingInsertions
                            / alignedGenomeLength,
                    (double) totalEvents / alignedGenomeLength,

                    genomeLength,
                    referenceGenomeLength,
                    alignedGenomeLength,
                    getN50(sampleStats->sampleSequenceLengthDistribution,
                            genomeLength),
                    getN50(sampleStats->blockLengthDistribution, genomeLength),
                    getN50(sampleStats->contigPathLengthDistribution,
                            genomeLength),
                    getN50(sampleStats->scaffoldPathLengthDistribution,
                            genomeLength),
                    getN50(sampleStatsFreeIndels->contigPathLengthDistribution,
                            genomeLength),
                    getN50(
                            sampleStatsFreeIndels->scaffoldPathLengthDistribution,
                            genomeLength),
                    stList_length(sampleStats->sampleSequenceLengthDistribution),
                    stList_length(sampleStats->blockLengthDistribution),
                    stList_length(sampleStats->contigPathLengthDistribution),
                    stList_length(sampleStats->scaffoldPathLengthDistribution),
                    stList_length(
                            sampleStatsFreeIndels->scaffoldPathLengthDistribution));

            stSortedSetIterator *it = stSortedSet_getIterator(indelEvents);
            IndelEvent *indelEvent;
            fprintf(fileHandle, "insertionLengthDistribution=\"");
            while ((indelEvent = stSortedSet_getNext(it)) != NULL) {
                if (indelEvent->insertionLength > 0) {
                    fprintf(fileHandle, "%i ", indelEvent->insertionLength);
                }
            }
            fprintf(fileHandle, "\" deletionLengthDistribution=\"");
            while ((indelEvent = stSortedSet_getPrevious(it)) != NULL) {
                if (indelEvent->deletionLength > 0) {
                    fprintf(fileHandle, "%i ", indelEvent->deletionLength);
                }
            }
            fprintf(fileHandle, "\">\n");
            while ((indelEvent = stSortedSet_getNext(it)) != NULL) {
                printIndel(indelEvent, eventStrings, fileHandle);
            }
            stSortedSet_destructIterator(it);
            stSortedSet_destruct(indelEvents);
            fprintf(fileHandle, "</statsForSample>\n");
            sampleStats_destruct(sampleStats);
        }
    }
    eventTree_destructIterator(eventIt);
}

void reportDistanceMatrix(Flower *flower, FILE *fileHandle,
        CapCodeParameters *capCodeParameters) {
    EventTree_Iterator *eventIt = eventTree_getIterator(
            flower_getEventTree(flower));
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        const char *eventString = event_getHeader(event);
        if (eventString != NULL && strcmp(eventString, "ROOT") != 0 && strcmp(
                eventString, "") != 0) {
            EventTree_Iterator *eventIt2 = eventTree_getIterator(
                    flower_getEventTree(flower));
            Event *event2;
            while ((event2 = eventTree_getNext(eventIt2)) != NULL) {
                const char *eventString2 = event_getHeader(event2);
                if (eventString2 != NULL && strcmp(eventString2, "ROOT") != 0
                        && strcmp(eventString2, "") != 0) {
                    SampleStats *sampleStats = getSamplePathStats(flower,
                            eventString, eventString2, capCodeParameters);
                    //Return total length of contig paths, total number of insertion, total number of deletions, total number of indels and total number of indels per base.
                    int32_t alignedGenomeLength = getSum(
                            sampleStats->blockLengthDistribution);

                    fprintf(
                            fileHandle,
                            "<indelDistanceForEvents "
                                "eventName1=\"%s\" "
                                "eventName2=\"%s\" "
                                "alignmentLength=\"%i\" "
                                "totalInsertions=\"%i\" "
                                "totalDeletions=\"%i\" "
                                "totalIndels=\"%i\" "
                                "indelsPerBase=\"%f\" "
                                "/>\n",
                            eventString,
                            eventString2,
                            alignedGenomeLength,
                            sampleStats->totalInsertions,
                            sampleStats->totalDeletions,
                            sampleStats->totalInsertions
                                    + sampleStats->totalDeletions,
                            ((double) sampleStats->totalInsertions
                                    + sampleStats->totalDeletions)
                                    / alignedGenomeLength);
                    sampleStats_destruct(sampleStats);
                }
            }
            eventTree_destructIterator(eventIt2);
        }
    }
    eventTree_destructIterator(eventIt);
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "pathStats");

    ///////////////////////////////////////////////////////////////////////////
    // Now print the haplotype path stats.
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<pathStats>\n");
    reportPathStatsForReference(flower, fileHandle, referenceEventString,
            capCodeParameters);
    reportDistanceMatrix(flower, fileHandle, capCodeParameters);
    fprintf(fileHandle, "</pathStats>\n");
    fclose(fileHandle);
    st_logInfo("Got the stats in %i seconds/\n", time(NULL) - startTime);

    return 0;
}
