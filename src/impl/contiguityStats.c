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
    int32_t *correct = st_malloc(sizeof(int32_t) * bucketNumber);
    int32_t *samples = st_malloc(sizeof(int32_t) * bucketNumber);

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<contiguityStats>");

    EventTree_Iterator *eventIt = eventTree_getIterator(flower_getEventTree(flower));
    Event *event;
    while((event = eventTree_getNext(eventIt)) != NULL) {
        const char *eventString = event_getHeader(event);
        if(eventString != NULL && strcmp(eventString, referenceEventString) != 0) {
            //Set the buckets up.
            for (int32_t i = 0; i < bucketNumber; i++) {
                correct[i] = 0.0;
                samples[i] = 0.0;
            }

            stList *eventStrings = stList_construct(); //This is the holder of the event strings
            stList_append(eventStrings, (char *)eventString);
            stSortedSet *sequences = getMetaSequencesForEvents(flower, eventStrings);
            stSortedSetIterator *it = stSortedSet_getIterator(sequences);
            MetaSequence *metaSequence;
            while ((metaSequence = stSortedSet_getNext(it)) != NULL) {
                samplePoints(flower, metaSequence, referenceEventString,
                        sampleNumber, correct, samples,
                        bucketNumber, bucketSize);
            }
            stSortedSet_destructIterator(it);
            stSortedSet_destruct(sequences);
            stList_destruct(eventStrings);

            //Now report the result
            fprintf(fileHandle, "<statsForSample sampleName=\"%s\" referenceName=\"%s\">\n", eventString, referenceEventString);
            int32_t pMaxSize = 1;
            int32_t cumulativeCorrect = 0;
            int32_t cumulativeSamples = 0;
            for (int32_t i = 0; i < bucketNumber; i++) {
                if(samples[i] > 0) {
                    assert(correct[i] <= samples[i]);
                    assert(correct[i] >= 0);
                    int32_t maxSize = pow(10, i
                            /bucketSize);
                    cumulativeCorrect += correct[i];
                    cumulativeSamples += samples[i];
                    fprintf(fileHandle, "\t<bucket from=\"%i\" to=\"%i\" correct=\"%i\" samples=\"%i\" correctPerSample=\"%f\" cumulativeCorrect=\"%i\" cumulativeSamples=\"%i\" cumulativeCorrectPerSample=\"%f\"/>\n",
                            pMaxSize, maxSize, correct[i], samples[i], ((float) correct[i])/samples[i], cumulativeCorrect, cumulativeSamples, ((float) cumulativeCorrect)/cumulativeSamples);
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

