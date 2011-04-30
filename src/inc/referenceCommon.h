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
 * Global parameters shared by all the scripts.
 */
extern CapCodeParameters *capCodeParameters;
extern char *outputFile;
extern Flower *flower;
extern CactusDisk *cactusDisk;

extern char *referenceEventString;

/*
 * Optional parameter used by copy number and substitution scripts.
 */
extern int32_t minimumBlockLength;

/*
 * Parameters for the substitution script.
 */
extern int32_t ignoreFirstNBasesOfBlock;
extern int32_t minimumIndentity;
extern bool printIndelPositions;

/*
 * For the linkage script.
 */
extern int32_t bucketNumber;
extern int32_t upperLinkageBound;
extern int32_t sampleNumber;

stList *getEventStrings(const char *hapA1EventString, const char *hapA2EventString);

void basicUsage(const char *programName);

int parseBasicArguments(int argc, char *argv[], const char *programName);

#endif /* COMMON_H_ */