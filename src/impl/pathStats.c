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

    stList *insertionDistribution;
    stList *deletionDistribution;

    stList *contigPathLengthDistribution;
    stList *scaffoldPathLengthDistribution;
    stList *blockLengthDistribution;
    stList *sampleSequenceLengthDistribution;
} SampleStats;

SampleStats *sampleStats_construct() {
    SampleStats *sampleStats = st_calloc(1, sizeof(SampleStats));
    sampleStats->insertionDistribution = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    sampleStats->deletionDistribution = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    sampleStats->contigPathLengthDistribution = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    sampleStats->scaffoldPathLengthDistribution = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    sampleStats->blockLengthDistribution = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    sampleStats->sampleSequenceLengthDistribution = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    return sampleStats;
}

void sampleStats_destruct(SampleStats *sampleStats) {
    stList_destruct(sampleStats->insertionDistribution);
    stList_destruct(sampleStats->deletionDistribution);
    stList_destruct(sampleStats->contigPathLengthDistribution);
    stList_destruct(sampleStats->scaffoldPathLengthDistribution);
    stList_destruct(sampleStats->blockLengthDistribution);
    stList_destruct(sampleStats->sampleSequenceLengthDistribution);
    free(sampleStats);
}

void getHaplotypePathStatsP(Cap *cap, stList *referenceEventStrings, stList *contaminationEventStrings, CapCodeParameters *capCodeParameters, SampleStats *sampleStats) {
    int32_t insertLength, deleteLength;
    switch (getCapCode(cap, referenceEventStrings, contaminationEventStrings, &insertLength, &deleteLength, capCodeParameters)) {
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
            sampleStats->totalScaffoldGaps++;
            return;
        case AMBIGUITY_GAP:
            sampleStats->totalAmbiguityGaps++;
            return;
        case ERROR_HAP_TO_HAP_SAME_CHROMOSOME:
            sampleStats->totalIntraJoins++;
            return;
        case ERROR_HAP_TO_HAP_DIFFERENT_CHROMOSOMES:
            sampleStats->totalInterJoins++;
            return;
        case ERROR_HAP_TO_INSERT:
            assert(insertLength > 0);
            stList_append(sampleStats->insertionDistribution, stIntTuple_construct(1, insertLength));
            sampleStats->totalInsertions++;
            return;
        case ERROR_HAP_TO_DELETION:
            assert(deleteLength > 0);
            stList_append(sampleStats->deletionDistribution, stIntTuple_construct(1, deleteLength));
            sampleStats->totalDeletions++;
            return;
        case ERROR_HAP_TO_INSERT_AND_DELETION:
            assert(insertLength > 0);
            assert(deleteLength > 0);
            stList_append(sampleStats->insertionDistribution, stIntTuple_construct(1, insertLength));
            stList_append(sampleStats->deletionDistribution, stIntTuple_construct(1, deleteLength));
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

void getSequenceLengths(Flower *flower, const char *eventString, stList *sequenceLengths) {
    Flower_SegmentIterator *sequenceIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while((sequence = flower_getNextSequence(sequenceIt)) != NULL) {
        if(strcmp(event_getHeader(sequence_getEvent(sequence)), eventString) == 0) {
            stList_append(sequenceLengths, stIntTuple_construct(1, sequence_getLength(sequence)));
        }
    }
    flower_destructSequenceIterator(sequenceIt);
}

void getBlockLengths(Flower *flower,
        const char *sampleEventString,
        const char *referenceEventString,
        stList *blockLengths) {
    //Find blocks..
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIt)) != NULL) {
        if(hasCapInEvent(block_get5End(block), sampleEventString) && hasCapInEvent(block_get5End(block), referenceEventString)) {
            stList_append(blockLengths, stIntTuple_construct(1, block_getLength(block)));
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
    for(int32_t i=0; i<stList_length(scaffoldPathsList);i++) {
        stSortedSet *scaffoldPath = stList_get(scaffoldPathsList, i);
        stSortedSetIterator *contigPathIt = stSortedSet_getIterator(scaffoldPath);
        stList *contigPath;
        int32_t j=0;
        while((contigPath = stSortedSet_getNext(contigPathIt)) != NULL) {
            j += contigPathLength(contigPath);
        }
        stSortedSet_destructIterator(contigPathIt);
        stList_append(scaffoldPathLengths, stIntTuple_construct(1, j));
    }
    stList_destruct(scaffoldPathsList);
}

SampleStats *getSamplePathStats(Flower *flower,
        const char *sampleEventString,
        const char *referenceEventString,
        CapCodeParameters *capCodeParameters) {
    /*
     * Gets stats for a given sample.
     */

    SampleStats *sampleStats = sampleStats_construct();

    stList *referenceEventStrings = stList_construct(); //This is the holder of the event strings
    stList_append(referenceEventStrings, (char *)referenceEventString);
    stList *emptyList = stList_construct();

    //Calculate the haplotype paths.
    stList *contigPaths = getContigPaths(flower, sampleEventString, referenceEventStrings);
    stHash *scaffoldPaths = getScaffoldPaths(contigPaths, referenceEventStrings, emptyList, capCodeParameters);

    //Now calculate the stats for the ends of the contig paths.
    for (int32_t i = 0; i < stList_length(contigPaths); i++) {
        stList *contigPath = stList_get(contigPaths, i);
        stList_append(sampleStats->contigPathLengthDistribution, stIntTuple_construct(1, contigPathLength(contigPath)));
        for (int32_t j = 0; j < stList_length(contigPath); j++) {
            Segment *segment = stList_get(contigPath, j);
            getHaplotypePathStatsP(segment_get5Cap(segment), referenceEventStrings, emptyList, capCodeParameters, sampleStats);
            getHaplotypePathStatsP(segment_get3Cap(segment), referenceEventStrings, emptyList, capCodeParameters, sampleStats);
        }
    }
    //Normalise stats
    assert(sampleStats->totalScaffoldGaps % 2 == 0);
    assert(sampleStats->totalAmbiguityGaps % 2 == 0);
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

int32_t getN50(stList *lengths, int32_t genomeLength) {
    int32_t totalLength = 0;
    int32_t pJ = INT32_MAX;
    stList_sort(lengths, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    for (int32_t i = stList_length(lengths)-1; i >= 0; i--) {
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
    for(int32_t i=0; i<stList_length(lengths); i++) {
        j += stIntTuple_getPosition(stList_get(lengths, i), 0);
    }
    return j;
}

int32_t getGenomeLength(Flower *flower, const char *eventString) {
    stList *sequences = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    getSequenceLengths(flower, eventString, sequences);
    int32_t i = getSum(sequences);
    stList_destruct(sequences);
    return i;
}

char *concatenateList(stList *list) {
    char **cAA = st_malloc(sizeof(char *) * (stList_length(list)/2));
    stList_sort(list, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    assert(stList_length(list) % 2 == 0);
    for (int32_t i = 0; i < stList_length(list); i+=2) {
        assert(stIntTuple_getPosition(stList_get(list, i), 0) == stIntTuple_getPosition(stList_get(list, i+1), 0));
        cAA[i/2] = stString_print("%i", stIntTuple_getPosition(stList_get(list, i), 0));
    }
    char *cA = stString_join(" ", (const char **)cAA, stList_length(list)/2);
    for (int32_t i = 0; i < stList_length(list)/2; i++) {
        free(cAA[i]);
    }
    free(cAA);
    return cA;
}

void reportPathStatsForReference(Flower *flower, FILE *fileHandle,
        const char *referenceEventString,
        CapCodeParameters *capCodeParameters) {

    int32_t referenceGenomeLength = getGenomeLength(flower, referenceEventString);

    EventTree_Iterator *eventIt = eventTree_getIterator(flower_getEventTree(flower));
    Event *event;
    while((event = eventTree_getNext(eventIt)) != NULL) {
        const char *eventString = event_getHeader(event);
        if(eventString != NULL && strcmp(eventString, referenceEventString) != 0) {
            SampleStats *sampleStats = getSamplePathStats(flower, referenceEventString, eventString, capCodeParameters);

            int32_t genomeLength = getGenomeLength(flower, eventString);

            char *insertionDistributionString = concatenateList(sampleStats->insertionDistribution);
            char *deletionDistributionString = concatenateList(sampleStats->deletionDistribution);

            fprintf(fileHandle, "<statsForSample "
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
                    "totalSampleGenomeLength=\"%i\" "
                    "totalReferenceGenomeLength=\"%i\" "
                    "sequenceN50=\"%i\" "
                    "blockN50=\"%i\" "
                    "contigPathN50=\"%i\" "
                    "scaffoldPathN50=\"%i\" "
                    "totalSequenceNumber=\"%i\" "
                    "totalBlockNumber=\"%i\" "
                    "totalContigPaths=\"%i\" "
                    "totalScaffoldPaths=\"%i\" "
                    "insertionSizeDistribution=\"%s\" "
                    "deletionSizeDistribution=\"%s\"/>\n",
                    eventString,
                    referenceEventString,
                    sampleStats->totalScaffoldGaps + sampleStats->totalAmbiguityGaps,
                    sampleStats->totalCleanEnds,
                    sampleStats->totalHangingEndWithsNs,
                    sampleStats->totalIntraJoins,
                    sampleStats->totalInterJoins,
                    sampleStats->totalInsertions,
                    sampleStats->totalDeletions,
                    sampleStats->totalInsertionAndDeletions,
                    sampleStats->totalHangingInsertions,
                    sampleStats->totalIntraJoins +
                    sampleStats->totalInterJoins +
                    sampleStats->totalInsertions +
                    sampleStats->totalDeletions +
                    sampleStats->totalInsertionAndDeletions +
                    sampleStats->totalHangingInsertions,
                    genomeLength,
                    referenceGenomeLength,
                    getN50(sampleStats->sampleSequenceLengthDistribution, genomeLength),
                    getN50(sampleStats->blockLengthDistribution, genomeLength),
                    getN50(sampleStats->contigPathLengthDistribution, genomeLength),
                    getN50(sampleStats->scaffoldPathLengthDistribution, genomeLength),
                    stList_length(sampleStats->sampleSequenceLengthDistribution),
                    stList_length(sampleStats->blockLengthDistribution),
                    stList_length(sampleStats->contigPathLengthDistribution),
                    stList_length(sampleStats->scaffoldPathLengthDistribution),
                    insertionDistributionString,
                    deletionDistributionString);
            sampleStats_destruct(sampleStats);
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
    reportPathStatsForReference(flower, fileHandle, referenceEventString, capCodeParameters);
    fprintf(fileHandle, "</pathStats>\n");
    fclose(fileHandle);
    st_logInfo("Got the stats in %i seconds/\n", time(NULL) - startTime);

    return 0;
}
