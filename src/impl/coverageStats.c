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

static void getMAFBlock2(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        stSortedSet *otherSampleEvents =
                    stSortedSet_construct3((int (*)(const void *, const void *))strcmp, NULL);
        Segment *segment;
        Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
        int32_t sampleNumber = 0;
        while ((segment = block_getNext(instanceIt)) != NULL) {
            const char *segmentEvent = event_getHeader(
                    segment_getEvent(segment));
            if (strcmp(segmentEvent, sampleEventString) == 0) {
                sampleNumber++;
            }
            if (strcmp(segmentEvent, referenceEventString) != 0) {
                stSortedSet_insert(otherSampleEvents, (void *)segmentEvent);
            }
        }
        block_destructInstanceIterator(instanceIt);
        baseCoverages[stSortedSet_size(otherSampleEvents)] += block_getLength(block) * sampleNumber;
        stSortedSet_destruct(otherSampleEvents);
    }
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "coverageStats");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<coverageStats>\n");
    EventTree_Iterator *eventIt = eventTree_getIterator(
            flower_getEventTree(flower));
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        sampleEventString = event_getHeader(event);
        if (sampleEventString != NULL && strcmp(sampleEventString, referenceEventString)
                != 0) {

            int32_t eventNumber = eventTree_getEventNumber(flower_getEventTree(flower));
            baseCoverages = st_calloc(sizeof(int32_t), eventNumber+1);

            getMAFsReferenceOrdered(flower, fileHandle, getMAFBlock2);

            ///////////////////////////////////////////////////////////////////////////
            // Print outputs
            ///////////////////////////////////////////////////////////////////////////

            fprintf(fileHandle, "<statsForSample "
                    "sampleName=\"%s\" "
                    "referenceName=\"%s\" "
                    "baseCoverages=\"",
                sampleEventString, referenceEventString);
            for(int32_t i=1; i<eventNumber; i++) {
                fprintf(fileHandle, "%i ", baseCoverages[i]);
            }
            fprintf(fileHandle, "\"/>\n");
            free(baseCoverages);
        }
    }
    eventTree_destructIterator(eventIt);
    fprintf(fileHandle, "</coverageStats>\n");

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}

