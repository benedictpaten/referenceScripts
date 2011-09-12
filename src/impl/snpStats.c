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

static void getSnpStats(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        //Now get the column
        Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
        Segment *segment;
        char *referenceSeq = NULL;
        Segment *referenceSegment = NULL;
        char *otherReferenceSeq = NULL;
        char *sampleSeq = NULL;
        Segment *sampleSegment = NULL;
        while ((segment = block_getNext(instanceIterator)) != NULL) {
            if (strcmp(event_getHeader(segment_getEvent(segment)), referenceEventString) == 0) {
                if (referenceSeq != NULL) {
                    goto end;
                }
                referenceSeq = segment_getString(segment);
                referenceSegment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)), sampleEventString) == 0) {
                if (sampleSeq != NULL) {
                    goto end;
                }
                sampleSeq = segment_getString(segment);
                sampleSegment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)), otherReferenceEventString) == 0) {
                if (otherReferenceSeq != NULL) {
                    goto end;
                }
                otherReferenceSeq = segment_getString(segment);
            }
        }

        if (referenceSeq != NULL && otherReferenceSeq != NULL) {
            //We're in gravy.
            for (int32_t i = ignoreFirstNBasesOfBlock; i < block_getLength(block) - ignoreFirstNBasesOfBlock; i++) {
                totalSites++;
                if (sampleSeq != NULL) {
                    totalCorrect += bitsScoreFn(sampleSeq[i], referenceSeq[i]);
                    bool b = correctFn(sampleSeq[i], referenceSeq[i]);
                    if(!b) {
                        assert(sampleSegment != NULL);
                        stList_append(events, sampleSegment);
                        assert(referenceSegment != NULL);
                        stList_append(events, referenceSegment);
                        stList_append(events, stIntTuple_construct(3, i, (int32_t)sampleSeq[i], (int32_t)referenceSeq[i]));
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
        if(otherReferenceSeq != NULL) {
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
    if(otherReferenceEventString == NULL) {
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
                "totalSites=\"%i\" "
                "totalCorrect=\"%f\" "
                "totalErrors=\"%i\" "
                "totalCalls=\"%i\">", sampleEventString, referenceEventString, otherReferenceEventString,
                    totalSites, totalCorrect, totalErrors, totalCalls);

            for(int32_t i=0; i<stList_length(events); i+=3) {
                Segment *sampleSegment = stList_get(events, i);
                assert(segment_getSequence(sampleSegment) != NULL);
                Segment *referenceSegment = stList_get(events, i+1);
                assert(segment_getSequence(referenceSegment) != NULL);
                int32_t coordinate = stIntTuple_getPosition(stList_get(events, i+2), 0);
                char base1 = stIntTuple_getPosition(stList_get(events, i+2), 1);
                char base2 = stIntTuple_getPosition(stList_get(events, i+2), 2);
                fprintf(fileHandle, "%s %i %c %s %i %c\n",
                        sequence_getHeader(segment_getSequence(sampleSegment)), segment_getStart(sampleSegment) + (segment_getStrand(sampleSegment) ? coordinate : -coordinate) - sequence_getStart(segment_getSequence(sampleSegment)), base1,
                        sequence_getHeader(segment_getSequence(referenceSegment)), segment_getStart(referenceSegment) + (segment_getStrand(referenceSegment) ? coordinate : -coordinate) - sequence_getStart(segment_getSequence(referenceSegment)), base2);
                stIntTuple_destruct(stList_get(events, i+2));
            }
            stList_destruct(events);
            fprintf(fileHandle, "</statsForSample>\n");
        }
    }
    eventTree_destructIterator(eventIt);
    fprintf(fileHandle, "</substitutionStats>\n");

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    flower = NULL;
    cactusDisk_destruct(cactusDisk);

    //while(1);

    return 0;
}
