/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "sonLib.h"
#include "cactus.h"
#include "adjacencyClassification.h"
#include "referenceCommon.h"
#include "linkage.h"

/*
 * Function used to account for sequence within adjacencies.
 */

int32_t getTotalLengthOfAdjacencies(Flower *flower, const char *eventName) {
    stList *eventStrings = stList_construct();
    stList_append(eventStrings, (void *) eventName);

    int64_t totalLength = 0;

    //First sum the length of the sequences
    stSortedSet *metaSequences = getMetaSequencesForEvents(flower, eventStrings);
    stSortedSetIterator *metaSequenceIterator = stSortedSet_getIterator(metaSequences);
    MetaSequence *metaSequence;
    while ((metaSequence = stSortedSet_getNext(metaSequenceIterator)) != NULL) {
        totalLength += metaSequence_getLength(metaSequence);
    }
    stSortedSet_destructIterator(metaSequenceIterator);
    stSortedSet_destruct(metaSequences);
    //Now subtract of the length of the sequences
    stSortedSet *orderedSegments = getOrderedSegments(flower);
    stSortedSetIterator *segmentIt = stSortedSet_getIterator(orderedSegments);
    Segment *segment;
    while ((segment = stSortedSet_getNext(segmentIt)) != NULL) {
        if (strcmp(event_getHeader(segment_getEvent(segment)), eventName) == 0) {
            totalLength -= segment_getLength(segment);
        }
    }
    stSortedSet_destructIterator(segmentIt);
    stSortedSet_destruct(orderedSegments);

    stList_destruct(eventStrings);
    assert(totalLength >= 0);
    return totalLength;
}

/*
 * Global parameters.
 */

CapCodeParameters *capCodeParameters = NULL;
char *outputFile = NULL;
Flower *flower = NULL;
CactusDisk *cactusDisk = NULL;

char *referenceEventString = NULL;
char *otherReferenceEventString = NULL;

/*
 * Optional parameter used by copy number and substitution scripts.
 */
int32_t minimumBlockLength = 0;

/*
 * Parameters for the substitution script.
 */
int32_t ignoreFirstNBasesOfBlock = 0;
int32_t minimumIndentity = 0;
bool printIndelPositions = 0;

/*
 * For the linkage script.
 */
int32_t bucketNumber = 200;
int32_t upperLinkageBound = 200000000;
int32_t sampleNumber = 1000000;
bool doNotSampleDuplicatedPositions = 0;

void basicUsage(const char *programName) {
    fprintf(stderr, "%s\n", programName);
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --cactusDisk : The location of the cactus disk\n");
    fprintf(stderr, "-e --outputFile : The file to write the output in.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr,
            "-m --minimumNsForScaffoldGap : Minimum number of Ns in assembly sequence to denote a scaffold gap\n");
    fprintf(stderr, "-n --maximumDeletionLength : Maximum length of deletion\n");
    fprintf(stderr, "-o --maximumInsertionLength : Maximum length of insertion\n");

    fprintf(stderr, "-p --referenceEventString : The reference event string\n");

    //Parameters not used by everyone..
    fprintf(stderr, "The following are not used by all scripts\n");
    fprintf(stderr, "-t --minimumBlockLength : Minimum block length\n");
    fprintf(stderr, "-u --ignoreFirstNBasesOfBlock : Minimum block length\n");
    fprintf(stderr, "-v --minimumIdentity : Minimum identity of the block\n");
    fprintf(stderr, "-w --printIndelPositions : Print out valid columns containing only one haplotype\n");
    fprintf(stderr, "-x --bucketNumber : Number of buckets\n");
    fprintf(stderr, "-y --upperLinkageBound : Upper linkage bound\n");
    fprintf(stderr, "-z --sampleNumber : Number of samples\n");
    fprintf(stderr, "-A --doNotSampleDuplicatedPositions : Do not sample positions that contain duplication\n");
    fprintf(
            stderr,
            "-B --otherReferenceEventString : The other reference event string (only look at sites that also contain this reference)\n");
}

int parseBasicArguments(int argc, char *argv[], const char *programName) {
    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    int32_t k;
    capCodeParameters = capCodeParameters_construct(25, INT32_MAX, 100000);
    referenceEventString = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "cactusDisk",
                required_argument, 0, 'c' }, { "outputFile", required_argument, 0, 'e' },
                { "help", no_argument, 0, 'h' }, { "minimumNsForScaffoldGap", required_argument, 0, 'm' }, {
                        "maximumDeletionLength", required_argument, 0, 'n' }, { "maximumInsertionLength",
                        required_argument, 0, 'o' }, { "referenceEventString", required_argument, 0, 'p' }, {
                        "minimumBlockLength", required_argument, 0, 't' }, { "ignoreFirstNBasesOfBlock",
                        required_argument, 0, 'u' }, { "minimumIdentity", required_argument, 0, 'v' }, {
                        "printIndelPositions", no_argument, 0, 'w' }, { "bucketNumber", required_argument, 0, 'x' }, {
                        "upperLinkageBound", required_argument, 0, 'y' },
                { "sampleNumber", required_argument, 0, 'z' },
                { "doNotSampleDuplicatedPositions", no_argument, 0, 'A' }, { "otherReferenceEventString",
                        required_argument, 0, 'B' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:e:hm:n:o:p:t:u:v:wx:y:z:AB:", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'e':
                outputFile = stString_copy(optarg);
                break;
            case 'h':
                basicUsage(programName);
                return 0;
            case 'm':
                k = sscanf(optarg, "%i", &capCodeParameters->minimumNCount);
                assert(k == 1);
                break;
            case 'n':
                k = sscanf(optarg, "%i", &capCodeParameters->maxDeletionLength);
                assert(k == 1);
                break;
            case 'o':
                k = sscanf(optarg, "%i", &capCodeParameters->maxInsertionLength);
                assert(k == 1);
                break;
            case 'p':
                referenceEventString = stString_copy(optarg);
                break;
            case 't':
                k = sscanf(optarg, "%i", &minimumBlockLength);
                assert(k == 1);
                break;
            case 'u':
                k = sscanf(optarg, "%i", &ignoreFirstNBasesOfBlock);
                assert(k == 1);
                break;
            case 'v':
                k = sscanf(optarg, "%i", &minimumIndentity);
                assert(k == 1);
                if (minimumIndentity > 100 || minimumIndentity < 0) {
                    st_errAbort("The minimum identity was not in the range [0, 100]: %i", minimumIndentity);
                }
                break;
            case 'w':
                printIndelPositions = 1;
                break;
            case 'x':
                k = sscanf(optarg, "%i", &bucketNumber);
                assert(k == 1);
                if (bucketNumber < 1) {
                    st_errAbort("The number of buckets can not be less than 1: %i", bucketNumber);
                }
                break;
            case 'y':
                k = sscanf(optarg, "%i", &upperLinkageBound);
                assert(k == 1);
                if (upperLinkageBound < 1) {
                    st_errAbort("The total interval of linkage is too small: %i", upperLinkageBound);
                }
                break;
            case 'z':
                k = sscanf(optarg, "%i", &sampleNumber);
                assert(k == 1);
                if (sampleNumber < 0) {
                    st_errAbort("The number of samples can not be less than 0: %i", sampleNumber);
                }
                break;
            case 'A':
                doNotSampleDuplicatedPositions = 1;
                break;
            case 'B':
                otherReferenceEventString = stString_copy(optarg);
                break;
            default:
                st_errAbort("Unrecognised option %s", optarg);
                break;
        }
    }

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    if (outputFile == NULL) {
        st_errAbort("The output file was not specified");
    }
    if (cactusDiskDatabaseString == NULL) {
        st_errAbort("The cactus disk string was not specified");
    }
    if (referenceEventString == NULL) {
        st_errAbort("The reference event string was not specified");
    }

    assert(outputFile != NULL);
    assert(cactusDiskDatabaseString != NULL);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Output graph file : %s\n", outputFile);
    st_logInfo("The cactus disk string : %s\n", cactusDiskDatabaseString);
    st_logInfo("The reference event string : %s\n", referenceEventString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the cactus disk\n");

    //////////////////////////////////////////////
    //Cleanup
    //////////////////////////////////////////////

    free(cactusDiskDatabaseString);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    flower = cactusDisk_getFlower(cactusDisk, 0);
    assert(flower != NULL);
    st_logInfo("Parsed the top level flower of the cactus tree\n");

    return 0;
}

