/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "cactus.h"
#include "cactusMafs.h"
#include "adjacencyTraversal.h"
#include "referenceCommon.h"
#include "linkage.h"

int64_t sumMetaSequenceLengths(stSortedSet *sequences) {
    int64_t j = 0;
    stSortedSetIterator *it = stSortedSet_getIterator(sequences);
    MetaSequence *i;
    while ((i = stSortedSet_getNext(it)) != NULL) {
        j += metaSequence_getLength(i);
    }
    stSortedSet_destructIterator(it);
    return j;
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "linkageStats");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    //Setup the buckets for the experiment
    double bucketSize = bucketNumber / log10(upperLinkageBound);
    int64_t *correct = st_malloc(sizeof(int64_t) * bucketNumber);
    int64_t *aligned = st_malloc(sizeof(int64_t) * bucketNumber);
    int64_t *samples = st_malloc(sizeof(int64_t) * bucketNumber);

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<contiguityStats>");

    stSortedSet *sortedSegments = getOrderedSegments(flower);

    EventTree_Iterator *eventIt = eventTree_getIterator(flower_getEventTree(flower));
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        const char *eventString = event_getHeader(event);
        if (eventString != NULL && strcmp(eventString, referenceEventString) != 0) {
            //Set the buckets up.
            for (int64_t i = 0; i < bucketNumber; i++) {
                correct[i] = 0;
                aligned[i] = 0;
                samples[i] = 0;
            }
            stList *eventStrings = stList_construct(); //This is the holder of the event strings
            stList_append(eventStrings, (char *) eventString);
            stSortedSet *sequences = getMetaSequencesForEvents(flower, eventStrings);
            int64_t metaSequencesLength = sumMetaSequenceLengths(sequences);
            stSortedSetIterator *it = stSortedSet_getIterator(sequences);
            MetaSequence *metaSequence;
            while ((metaSequence = stSortedSet_getNext(it)) != NULL) {
                if(metaSequence_getLength(metaSequence) > 20) {
                    int64_t i = sampleNumber * (((double) metaSequence_getLength(metaSequence)) / metaSequencesLength);
                    if(otherReferenceEventString != NULL) {
                        samplePointsWithOtherReference(flower, metaSequence, referenceEventString, otherReferenceEventString, i, correct, aligned, samples, bucketNumber, bucketSize,
                                                        sortedSegments, !doNotSampleDuplicatedPositions, 0.9);
                    }
                    else {
                        samplePoints(flower, metaSequence, referenceEventString, i, correct, aligned, samples, bucketNumber, bucketSize,
                                sortedSegments, !doNotSampleDuplicatedPositions, 0.9);
                    }
                }
            }
            stSortedSet_destructIterator(it);
            stSortedSet_destruct(sequences);
            stList_destruct(eventStrings);

            int64_t cumulativeCorrect = 0;
            int64_t cumulativeAligned = 0;
            int64_t cumulativeSamples = 0;
            for (int64_t i = 0; i < bucketNumber; i++) {
                if (samples[i] > 0) {
                    cumulativeCorrect += correct[i];
                    cumulativeAligned += aligned[i];
                    cumulativeSamples += samples[i];
                }
            }

            //Now report the result
            fprintf(fileHandle, "<statsForSample sampleName=\"%s\" referenceName=\"%s\" totalSamples=\"%" PRIi64 "\" totalCorrect=\"%" PRIi64 "\" totalAligned=\"%" PRIi64 "\" correctPerSample=\"%f\" correctPerAligned=\"%f\">\n", eventString, referenceEventString, cumulativeSamples, cumulativeCorrect, cumulativeAligned, ((double)cumulativeCorrect)/cumulativeSamples, ((double)cumulativeCorrect)/cumulativeAligned);
            int64_t pMaxSize = 1;
            cumulativeCorrect = 0;
            cumulativeAligned = 0;
            cumulativeSamples = 0;
            for (int64_t i = 0; i < bucketNumber; i++) {
                if (samples[i] > 0) {
                    assert(correct[i] <= samples[i]);
                    assert(correct[i] <= aligned[i]);
                    assert(aligned[i] <= samples[i]);
                    assert(correct[i] >= 0);
                    int64_t maxSize = pow(10, i / bucketSize);
                    cumulativeCorrect += correct[i];
                    cumulativeAligned += aligned[i];
                    cumulativeSamples += samples[i];
                    fprintf(
                            fileHandle,
                            "\t<bucket from=\"%" PRIi64 "\" to=\"%" PRIi64 "\" correct=\"%" PRIi64 "\" aligned=\"%" PRIi64 "\" samples=\"%" PRIi64 "\" correctPerSample=\"%f\" correctPerAligned=\"%f\" cumulativeCorrect=\"%" PRIi64 "\" cumulativeAligned=\"%" PRIi64 "\" cumulativeSamples=\"%" PRIi64 "\" cumulativeCorrectPerSample=\"%f\" cumulativeCorrectPerAligned=\"%f\"/>\n",
                            pMaxSize, maxSize, correct[i], aligned[i], samples[i], ((float) correct[i]) / samples[i], ((float) correct[i]) / aligned[i],
                            cumulativeCorrect, cumulativeAligned, cumulativeSamples, ((float) cumulativeCorrect) / cumulativeSamples, ((float) cumulativeCorrect) / cumulativeAligned);
                    pMaxSize = maxSize + 1;
                }
            }
            fprintf(fileHandle, "</statsForSample>");
        }
    }
    eventTree_destructIterator(eventIt);

    fprintf(fileHandle, "\n</contiguityStats>\n");
    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}

