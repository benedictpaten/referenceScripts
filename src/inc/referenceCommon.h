/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef COMMON_H_
#define COMMON_H_

#include "cactus.h"
#include "sonLib.h"
#include "adjacencyClassification.h"

/*
 * Returns the sum of bases contained within adjacencies for the sequences
 * with the given event names.
 */
int64_t getTotalLengthOfAdjacencies(Flower *flower, const char *eventName);

/*
 * Global parameters shared by all the scripts.
 */
extern CapCodeParameters *capCodeParameters;
extern char *outputFile;
extern Flower *flower;
extern CactusDisk *cactusDisk;

extern char *referenceEventString;
extern char *otherReferenceEventString;
extern char *outgroupEventString;

/*
 * Optional parameter used by copy number and substitution scripts.
 */
extern int64_t minimumBlockLength;

/*
 * Parameters for the substitution script.
 */
extern int64_t ignoreFirstNBasesOfBlock;
extern int64_t minimumIndentity;
extern bool printIndelPositions;

/*
 * Parameters for the path stats script.
 */
extern bool ignoreAdjacencies;

/*
 * For the linkage script.
 */
extern int64_t bucketNumber;
extern int64_t upperLinkageBound;
extern int64_t sampleNumber;
extern bool doNotSampleDuplicatedPositions;
extern bool ignoreSitesWithOtherReferencePresent;

/*
 * For the snps script
 */
extern int64_t minimumRecurrence;

stList *getEventStrings(const char *hapA1EventString, const char *hapA2EventString);

void basicUsage(const char *programName);

int parseBasicArguments(int argc, char *argv[], const char *programName);

#endif /* COMMON_H_ */
