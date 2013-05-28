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
    int64_t insertionLength;
    int64_t deletionLength;
} IndelEvent;

IndelEvent *indelEvent_construct(Cap *cap, Cap *otherCap, int64_t insertionLength, int64_t deletionLength) {
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
    int i = cactusMisc_nameCompare(cap_getName(indelEvent1->cap), cap_getName(indelEvent2->cap));
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
    int64_t totalCleanEnds;
    int64_t totalHangingEndWithsNs;
    int64_t totalScaffoldGaps;
    int64_t totalAmbiguityGaps;

    int64_t totalIntraJoins;
    int64_t totalInterJoins;
    int64_t totalInsertions;
    int64_t totalDeletions;
    int64_t totalInsertionAndDeletions;
    int64_t totalHangingInsertions;

    stList *indelEvents;

    stList *contigPathLengthDistribution;
    stList *scaffoldPathLengthDistribution;
    stList *blockLengthDistribution;
    stList *sampleSequenceLengthDistribution;
    stList *nonLinearRearrangements;
} SampleStats;

SampleStats *sampleStats_construct() {
    SampleStats *sampleStats = st_calloc(1, sizeof(SampleStats));
    sampleStats->indelEvents = stList_construct3(0, free);
    sampleStats->contigPathLengthDistribution = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    sampleStats->scaffoldPathLengthDistribution = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    sampleStats->blockLengthDistribution = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    sampleStats->sampleSequenceLengthDistribution = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    sampleStats->nonLinearRearrangements = stList_construct();
    return sampleStats;
}

void sampleStats_destruct(SampleStats *sampleStats) {
    stList_destruct(sampleStats->indelEvents);
    stList_destruct(sampleStats->contigPathLengthDistribution);
    stList_destruct(sampleStats->scaffoldPathLengthDistribution);
    stList_destruct(sampleStats->blockLengthDistribution);
    stList_destruct(sampleStats->sampleSequenceLengthDistribution);
    stList_destruct(sampleStats->nonLinearRearrangements);
    free(sampleStats);
}

void getHaplotypePathStatsP(Cap *cap, stList *referenceEventStrings, stList *contaminationEventStrings,
        CapCodeParameters *capCodeParameters, SampleStats *sampleStats) {
    int64_t insertLength, deleteLength;
    Cap *otherCap = NULL;
    switch (getCapCode(cap, &otherCap, referenceEventStrings, contaminationEventStrings, &insertLength, &deleteLength,
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
                stList_append(sampleStats->indelEvents, indelEvent_construct(cap, otherCap, insertLength, 0));
                sampleStats->totalInsertions++;
            }
            if (deleteLength > 0 && insertLength == 0) {
                stList_append(sampleStats->indelEvents, indelEvent_construct(cap, otherCap, 0, deleteLength));
                sampleStats->totalDeletions++;
            }
            if (insertLength > 0 && deleteLength > 0) {
                stList_append(sampleStats->indelEvents, indelEvent_construct(cap, otherCap, insertLength, deleteLength));
                sampleStats->totalInsertions++;
                sampleStats->totalDeletions++;
            }
            sampleStats->totalScaffoldGaps++;
            return;
        case ERROR_HAP_TO_HAP_SAME_CHROMOSOME:
            stList_append(sampleStats->nonLinearRearrangements, cap);
            stList_append(sampleStats->nonLinearRearrangements, otherCap);
            sampleStats->totalIntraJoins++;
            return;
        case ERROR_HAP_TO_HAP_DIFFERENT_CHROMOSOMES:
            sampleStats->totalInterJoins++;
            return;
        case ERROR_HAP_TO_INSERT:
            assert(insertLength > 0);
            stList_append(sampleStats->indelEvents, indelEvent_construct(cap, otherCap, insertLength, 0));
            sampleStats->totalInsertions++;
            return;
        case ERROR_HAP_TO_DELETION:
            assert(deleteLength > 0);
            stList_append(sampleStats->indelEvents, indelEvent_construct(cap, otherCap, 0, deleteLength));
            sampleStats->totalDeletions++;
            return;
        case ERROR_HAP_TO_INSERT_AND_DELETION:
            assert(insertLength > 0);
            assert(deleteLength > 0);
            stList_append(sampleStats->indelEvents, indelEvent_construct(cap, otherCap, insertLength, deleteLength));
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

static Cap *getCapForReferenceEventString(End *end, const char *referenceEventString) {
    /*
     * Get the cap for a given event.
     */
    End_InstanceIterator *it = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        if (strcmp(event_getHeader(cap_getEvent(cap)), referenceEventString) == 0) {
            end_destructInstanceIterator(it);
            return cap;
        }
    }
    end_destructInstanceIterator(it);
    assert(0);
    return NULL;
}

void printNonLinearRearrangement(Cap *cap, Cap *otherCap, FILE *fileHandle, const char *referenceEventString) {
    Cap *cap2 = getCapForReferenceEventString(cap_getEnd(cap), referenceEventString);
    Cap *otherCap2 = getCapForReferenceEventString(cap_getEnd(otherCap), referenceEventString);
    fprintf(fileHandle, "SEQ1: %s START1: %" PRIi64 " SEQ2: %s START2: %" PRIi64 " ", sequence_getHeader(cap_getSequence(cap)),
            cap_getCoordinate(cap) + 1 - sequence_getStart(cap_getSequence(cap)),
            sequence_getHeader(cap_getSequence(otherCap)),
            cap_getCoordinate(otherCap) + 1 - sequence_getStart(cap_getSequence(cap)));
    if (cap2 != NULL) {
        fprintf(fileHandle, "SEQ3: %s START3: %" PRIi64 " ", sequence_getHeader(cap_getSequence(cap2)),
                cap_getCoordinate(cap2) + 1 - sequence_getStart(cap_getSequence(cap2)));
    }
    if (otherCap2 != NULL) {
        fprintf(fileHandle, "SEQ4: %s START4: %" PRIi64 " ", sequence_getHeader(cap_getSequence(otherCap2)),
                cap_getCoordinate(otherCap2) + 1 - sequence_getStart(cap_getSequence(cap2)));
    }
    fprintf(fileHandle, "# ");
}

bool indelCanBeDescribedWithRespectToGivenReference(IndelEvent *indelEvent, const char *reference) {
    stList *eventStrings = stList_construct();
    stList_append(eventStrings, (char *)reference);
    Cap *cap3 = NULL, *cap4 = NULL;
    int64_t i = 0;
    bool b = endsAreAdjacent2(cap_getEnd(indelEvent->cap), cap_getEnd(indelEvent->otherCap), &cap3, &cap4, &i,
            eventStrings);
    stList_destruct(eventStrings);
    return b;
}

void printIndel(IndelEvent *indelEvent, const char *referenceEventString,
        const char *otherReferenceEventString, FILE *fileHandle) {
    stList *eventStrings = stList_construct();
    stList_append(eventStrings, (void *) referenceEventString);
    Cap *cap3 = NULL, *cap4 = NULL;
    int64_t i = 0;
    bool b = endsAreAdjacent2(cap_getEnd(indelEvent->cap), cap_getEnd(indelEvent->otherCap), &cap3, &cap4, &i,
            eventStrings);
    stList_destruct(eventStrings);
    (void)b;
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
            ignoreAdjacencies || cap_getCoordinate(indelEvent->otherCap) - cap_getCoordinate(indelEvent->cap) - 1
                    == indelEvent->insertionLength);
    assert(cap_getSequence(indelEvent->cap) != NULL);
    assert(cap_getSequence(cap3) != NULL);

    fprintf(
            fileHandle,
            "SEQ1: %s START1: %" PRIi64 " LENGTH1: %" PRIi64 " STRAND1: 1 SEQ2: %s START2: %" PRIi64 " LENGTH2: %" PRIi64 " STRAND2: %" PRIi64 " PRESENT_IN_OTHER_REFERENCE: %" PRIi64 "\n",
            sequence_getHeader(cap_getSequence(indelEvent->cap)),
            cap_getCoordinate(indelEvent->cap) + 1 - sequence_getStart(cap_getSequence(indelEvent->cap)),
            indelEvent->insertionLength,
            sequence_getHeader(cap_getSequence(cap3)),
            cap_getSide(cap3) ? cap_getCoordinate(cap3) - 1 - sequence_getStart(cap_getSequence(cap3))
                    : cap_getCoordinate(cap3) + 1 - sequence_getStart(cap_getSequence(cap3)),
            indelEvent->deletionLength, (int64_t)!cap_getSide(cap3), (int64_t)indelCanBeDescribedWithRespectToGivenReference(indelEvent, otherReferenceEventString));
}

void getSequenceLengths(Flower *flower, const char *eventString, stList *sequenceLengths) {
    Flower_SegmentIterator *sequenceIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        if (strcmp(event_getHeader(sequence_getEvent(sequence)), eventString) == 0) {
            stList_append(sequenceLengths, stIntTuple_construct1( sequence_getLength(sequence)));
        }
    }
    flower_destructSequenceIterator(sequenceIt);
}

void getBlockLengths(Flower *flower, const char *sampleEventString, const char *referenceEventString,
        stList *blockLengths) {
    //Find blocks..
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIt)) != NULL) {
        if (hasCapInEvent(block_get5End(block), sampleEventString) && hasCapInEvent(block_get5End(block),
                referenceEventString)) {
            stList_append(blockLengths, stIntTuple_construct1( block_getLength(block)));
        }
    }
    flower_destructBlockIterator(blockIt);

    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (!group_isLeaf(group)) {
            getBlockLengths(group_getNestedFlower(group), sampleEventString, referenceEventString, blockLengths);
        }
    }
    flower_destructGroupIterator(groupIt);
}

void getScaffoldPathsLengths(stHash *scaffoldPaths, stList *scaffoldPathLengths) {
    stList *scaffoldPathsList = stHash_getValues(scaffoldPaths);
    for (int64_t i = 0; i < stList_length(scaffoldPathsList); i++) {
        stSortedSet *scaffoldPath = stList_get(scaffoldPathsList, i);
        stSortedSetIterator *contigPathIt = stSortedSet_getIterator(scaffoldPath);
        stList *contigPath;
        int64_t j = 0;
        while ((contigPath = stSortedSet_getNext(contigPathIt)) != NULL) {
            j += contigPathLength(contigPath);
        }
        stSortedSet_destructIterator(contigPathIt);
        stList_append(scaffoldPathLengths, stIntTuple_construct1( j));
    }
    stList_destruct(scaffoldPathsList);
}

SampleStats *getSamplePathStats(Flower *flower, const char *sampleEventString, const char *referenceEventString,
        CapCodeParameters *capCodeParameters) {
    /*
     * Gets stats for a given sample.
     */

    SampleStats *sampleStats = sampleStats_construct();

    stList *referenceEventStrings = stList_construct(); //This is the holder of the event strings
    stList_append(referenceEventStrings, (char *) referenceEventString);
    stList *emptyList = stList_construct();

    //Calculate the haplotype paths.
    stList *contigPaths = getContigPaths(flower, sampleEventString, referenceEventStrings);
    stHash *scaffoldPaths = getScaffoldPaths(contigPaths, referenceEventStrings, emptyList, capCodeParameters);

    //Now calculate the stats for the ends of the contig paths.
    for (int64_t i = 0; i < stList_length(contigPaths); i++) {
        stList *contigPath = stList_get(contigPaths, i);
        stList_append(sampleStats->contigPathLengthDistribution, stIntTuple_construct1( contigPathLength(contigPath)));
        for (int64_t j = 0; j < stList_length(contigPath); j++) {
            Segment *segment = stList_get(contigPath, j);
            getHaplotypePathStatsP(segment_get5Cap(segment), referenceEventStrings, emptyList, capCodeParameters,
                    sampleStats);
            getHaplotypePathStatsP(segment_get3Cap(segment), referenceEventStrings, emptyList, capCodeParameters,
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
    getScaffoldPathsLengths(scaffoldPaths, sampleStats->scaffoldPathLengthDistribution);
    getBlockLengths(flower, sampleEventString, referenceEventString, sampleStats->blockLengthDistribution);
    getSequenceLengths(flower, sampleEventString, sampleStats->sampleSequenceLengthDistribution);

    //Cleanup
    stHash_destruct(scaffoldPaths);
    stList_destruct(contigPaths);
    stList_destruct(referenceEventStrings);
    stList_destruct(emptyList);

    return sampleStats;
}

int64_t getN50(stList *lengths, int64_t genomeLength) {
    int64_t totalLength = 0;
    int64_t pJ = INT64_MAX;
    stList_sort(lengths, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    for (int64_t i = stList_length(lengths) - 1; i >= 0; i--) {
        int64_t j = stIntTuple_get(stList_get(lengths, i), 0);
        assert(j <= pJ);
        pJ = j;
        totalLength += j;
        if (totalLength >= genomeLength / 2) {
            return j;
        }
    }
    return 0; //this should not happen!
}

int64_t getSum(stList *lengths) {
    int64_t j = 0;
    for (int64_t i = 0; i < stList_length(lengths); i++) {
        j += stIntTuple_get(stList_get(lengths, i), 0);
    }
    return j;
}

int64_t getGenomeLength(Flower *flower, const char *eventString) {
    stList *sequences = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    getSequenceLengths(flower, eventString, sequences);
    int64_t i = getSum(sequences);
    stList_destruct(sequences);
    return i;
}

void reportPathStatsForReference(Flower *flower, FILE *fileHandle, const char *referenceEventString, const char *otherReferenceEventString,
        CapCodeParameters *capCodeParameters) {

    int64_t referenceGenomeLength = getGenomeLength(flower, referenceEventString);

    EventTree_Iterator *eventIt = eventTree_getIterator(flower_getEventTree(flower));
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        const char *eventString = event_getHeader(event);
        if (eventString != NULL && strcmp(eventString, referenceEventString) != 0) {
            SampleStats *sampleStats = getSamplePathStats(flower, eventString, referenceEventString, capCodeParameters);
            int64_t minimumNCount = capCodeParameters->minimumNCount;
            capCodeParameters->minimumNCount = 0;
            SampleStats *sampleStatsFreeIndels = getSamplePathStats(flower, eventString, referenceEventString,
                    capCodeParameters);
            capCodeParameters->minimumNCount = minimumNCount;

            int64_t genomeLength = getGenomeLength(flower, eventString);
            int64_t alignedGenomeLength = getSum(sampleStats->blockLengthDistribution);

            stSortedSet *indelEvents = stList_getSortedSet(sampleStats->indelEvents,
                    (int(*)(const void *, const void *)) indelEvent_cmpFn);
            assert(stList_length(sampleStats->indelEvents) == 2 * stSortedSet_size(indelEvents));

            int64_t totalEvents = sampleStats->totalIntraJoins + sampleStats->totalInterJoins
                    + sampleStats->totalInsertions + sampleStats->totalDeletions
                    + sampleStats->totalInsertionAndDeletions + sampleStats->totalHangingInsertions;

            fprintf(fileHandle, "<statsForSample "
                "sampleName=\"%s\" "
                "referenceName=\"%s\" "
                "totalScaffoldGaps=\"%" PRIi64 "\" "
                "totalContigEnds=\"%" PRIi64 "\" "
                "totalContigEndsWithNs=\"%" PRIi64 "\" "
                "totalIntraJoin=\"%" PRIi64 "\" "
                "totalInterJoin=\"%" PRIi64 "\" "
                "totalInsertion=\"%" PRIi64 "\" "
                "totalDeletion=\"%" PRIi64 "\" "
                "totalInsertionAndDeletion=\"%" PRIi64 "\" "
                "totalHangingInsertion=\"%" PRIi64 "\" "
                "totalDifferences=\"%" PRIi64 "\" "

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

                "totalSampleGenomeLength=\"%" PRIi64 "\" "
                "totalReferenceGenomeLength=\"%" PRIi64 "\" "
                "totalAlignedGenomeLength=\"%" PRIi64 "\" "
                "sequenceN50=\"%" PRIi64 "\" "
                "blockN50=\"%" PRIi64 "\" "
                "contigPathN50=\"%" PRIi64 "\" "
                "scaffoldPathN50=\"%" PRIi64 "\" "
                "contigPathN50Check=\"%" PRIi64 "\" "
                "indelPathN50=\"%" PRIi64 "\" "
                "totalSequenceNumber=\"%" PRIi64 "\" "
                "totalBlockNumber=\"%" PRIi64 "\" "
                "totalContigPaths=\"%" PRIi64 "\" "
                "totalScaffoldPaths=\"%" PRIi64 "\" "
                "totalIndelPaths=\"%" PRIi64 "\" ", eventString, referenceEventString,
                    sampleStats->totalScaffoldGaps + sampleStats->totalAmbiguityGaps, sampleStats->totalCleanEnds,
                    sampleStats->totalHangingEndWithsNs, sampleStats->totalIntraJoins, sampleStats->totalInterJoins,
                    sampleStats->totalInsertions, sampleStats->totalDeletions, sampleStats->totalInsertionAndDeletions,
                    sampleStats->totalHangingInsertions, totalEvents,

                    (double) (sampleStats->totalScaffoldGaps + sampleStats->totalAmbiguityGaps) / alignedGenomeLength,
                    (double) sampleStats->totalCleanEnds / alignedGenomeLength,
                    (double) sampleStats->totalHangingEndWithsNs / alignedGenomeLength,
                    (double) sampleStats->totalIntraJoins / alignedGenomeLength,
                    (double) sampleStats->totalInterJoins / alignedGenomeLength,
                    (double) sampleStats->totalInsertions / alignedGenomeLength,
                    (double) sampleStats->totalDeletions / alignedGenomeLength,
                    (double) sampleStats->totalInsertionAndDeletions / alignedGenomeLength,
                    (double) sampleStats->totalHangingInsertions / alignedGenomeLength,
                    (double) totalEvents / alignedGenomeLength,

                    genomeLength, referenceGenomeLength, alignedGenomeLength,
                    getN50(sampleStats->sampleSequenceLengthDistribution, genomeLength),
                    getN50(sampleStats->blockLengthDistribution, genomeLength),
                    getN50(sampleStats->contigPathLengthDistribution, genomeLength),
                    getN50(sampleStats->scaffoldPathLengthDistribution, genomeLength),
                    getN50(sampleStatsFreeIndels->contigPathLengthDistribution, genomeLength),
                    getN50(sampleStatsFreeIndels->scaffoldPathLengthDistribution, genomeLength),
                    stList_length(sampleStats->sampleSequenceLengthDistribution),
                    stList_length(sampleStats->blockLengthDistribution),
                    stList_length(sampleStats->contigPathLengthDistribution),
                    stList_length(sampleStats->scaffoldPathLengthDistribution),
                    stList_length(sampleStatsFreeIndels->scaffoldPathLengthDistribution));

            stSortedSetIterator *it = stSortedSet_getIterator(indelEvents);
            IndelEvent *indelEvent;
            fprintf(fileHandle, "insertionSizeDistribution=\"");
            while ((indelEvent = stSortedSet_getNext(it)) != NULL) {
                if (indelEvent->insertionLength > 0) {
                    fprintf(fileHandle, "%" PRIi64 " ", indelEvent->insertionLength);
                }
            }
            fprintf(fileHandle, "\" deletionSizeDistribution=\"");
            while ((indelEvent = stSortedSet_getPrevious(it)) != NULL) {
                if (indelEvent->deletionLength > 0) {
                    fprintf(fileHandle, "%" PRIi64 " ", indelEvent->deletionLength);
                }
            }
            fprintf(fileHandle, "\" locationOfNonLinearBreakpoints=\"");
            for (int64_t i = 0; i < stList_length(sampleStats->nonLinearRearrangements); i += 2) {
                printNonLinearRearrangement(stList_get(sampleStats->nonLinearRearrangements, i),
                        stList_get(sampleStats->nonLinearRearrangements, i + 1), fileHandle, referenceEventString);
            }
            fprintf(fileHandle, "\">\n");
            while ((indelEvent = stSortedSet_getNext(it)) != NULL) {
                printIndel(indelEvent, referenceEventString, otherReferenceEventString, fileHandle);
            }
            stSortedSet_destructIterator(it);
            stSortedSet_destruct(indelEvents);
            fprintf(fileHandle, "</statsForSample>\n");
            sampleStats_destruct(sampleStats);
        }
    }
    eventTree_destructIterator(eventIt);
}

void reportDistanceMatrix(Flower *flower, FILE *fileHandle, CapCodeParameters *capCodeParameters) {
    EventTree_Iterator *eventIt = eventTree_getIterator(flower_getEventTree(flower));
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        const char *eventString = event_getHeader(event);
        if (eventString != NULL && strcmp(eventString, "ROOT") != 0 && strcmp(eventString, "") != 0) {
            EventTree_Iterator *eventIt2 = eventTree_getIterator(flower_getEventTree(flower));
            Event *event2;
            while ((event2 = eventTree_getNext(eventIt2)) != NULL) {
                const char *eventString2 = event_getHeader(event2);
                if (eventString2 != NULL && strcmp(eventString2, "ROOT") != 0 && strcmp(eventString2, "") != 0) {
                    SampleStats *sampleStats = getSamplePathStats(flower, eventString, eventString2, capCodeParameters);
                    //Return total length of contig paths, total number of insertion, total number of deletions, total number of indels and total number of indels per base.
                    int64_t alignedGenomeLength = getSum(sampleStats->blockLengthDistribution);

                    fprintf(fileHandle, "<indelDistanceForEvents "
                        "eventName1=\"%s\" "
                        "eventName2=\"%s\" "
                        "alignmentLength=\"%" PRIi64 "\" "
                        "totalInsertions=\"%" PRIi64 "\" "
                        "totalDeletions=\"%" PRIi64 "\" "
                        "totalIndels=\"%" PRIi64 "\" "
                        "indelsPerBase=\"%f\" "
                        "/>\n", eventString, eventString2, alignedGenomeLength, sampleStats->totalInsertions,
                            sampleStats->totalDeletions, sampleStats->totalInsertions + sampleStats->totalDeletions,
                            ((double) sampleStats->totalInsertions + sampleStats->totalDeletions) / alignedGenomeLength);
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

    getTerminalAdjacencyLength_ignoreAdjacencies = ignoreAdjacencies;
    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<pathStats>\n");
    if(otherReferenceEventString == NULL) {
        otherReferenceEventString = referenceEventString;
    }
    reportPathStatsForReference(flower, fileHandle, referenceEventString, otherReferenceEventString, capCodeParameters);
    if(makeDistanceMatrix) {
        reportDistanceMatrix(flower, fileHandle, capCodeParameters);
    }
    fprintf(fileHandle, "</pathStats>\n");
    fclose(fileHandle);
    st_logInfo("Got the stats in %" PRIi64 " seconds/\n", time(NULL) - startTime);

    return 0;
}
