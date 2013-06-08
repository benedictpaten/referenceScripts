/*
 *
 */

#include "assert.h"
#include "getopt.h"
#include "stdio.h"
#include "stdlib.h"

#include "sonLib.h"
#include "cactus.h"
#include "referenceCommon.h"

/*
 *Global parameters
 */
char *constraintFile = NULL;

/************************
 ******* CIGAR **********
 ************************/
struct _cigar{
    //cigar: hg19SEQ1 896 946 + 1kIndelContigsSEQ1 0 50 + 50 M 50
    char * header1;
    int64_t start1;
    int64_t end1;
    char strand1;
    char *header2;
    int64_t start2;
    int64_t end2;
    char strand2;
    int64_t score;
    char *cigarStr;
    stList *starts1; //start of matched blocks
    stList *ends1;
    stList *starts2;
    stList *ends2;
    int64_t correctBases;
    int numVisited;
};
typedef struct _cigar Cigar;

void cigar_setMatchedBlocks(Cigar *cigar){
    stList *items = stString_split(cigar->cigarStr);
    int64_t currstart1 = cigar->start1;
    int64_t currstart2 = cigar->start2;
    
    cigar->starts1 = stList_construct();
    cigar->ends1 = stList_construct();
    cigar->starts2 = stList_construct();
    cigar->ends2 = stList_construct();

    assert( stList_length(items)%2 == 0 );
    for(int64_t i=0; i < stList_length(items) - 1; i += 2){
        char * letter = stList_get(items, i);
        int64_t len = (int64_t) strtol( stList_get(items, i+1), (char **)NULL, 10 );
        if (strcmp(letter, "M") == 0){ //match
            int64_t s1 = currstart1;
            stList_append(cigar->starts1, s1);
            currstart1 += len;
            int64_t e1 = currstart1;
            stList_append(cigar->ends1, e1);
            
            int64_t s2 = currstart2;
            stList_append(cigar->starts2, s2);
            currstart2 += len;
            int64_t e2 = currstart2;
            stList_append(cigar->ends2, e2);
        }else if( strcmp(letter, "I") == 0 ){ //insertion
            currstart1 += len;
        }else if( strcmp(letter, "D") == 0 ) {
            currstart2 += len;
        }else {
            st_errAbort("Cigar string with letter that is different from [M, I, D]");
        }
    }
}

Cigar *cigar_construct(char *header1, int64_t start1, int64_t end1, char strand1,
                       char *header2, int64_t start2, int64_t end2, char strand2,
                       int64_t score, char *cigarStr)
{
    Cigar *cigar;
    cigar = st_malloc(sizeof(Cigar));
    cigar->header1 = stString_copy(header1);
    cigar->start1 = start1; 
    cigar->end1 = end1;
    cigar->strand1 = strand1;
    cigar->header2 = stString_copy(header2);
    cigar->start2 = start2;
    cigar->end2 = end2;
    cigar->strand2 = strand2;
    cigar->score = score;
    cigar->cigarStr = stString_copy(cigarStr);
    cigar_setMatchedBlocks(cigar);
    cigar->correctBases = 0;
    cigar->numVisited = 0;
    return cigar;
}

void cigar_destruct(Cigar *cigar){
    free(cigar);
}

/////////////////////////////////////

int cigar_cmp(Cigar *cig1, Cigar *cig2){ 
    //currently constraints have to be on positive strands
    assert( cig1.strand1 == '+' && cig1.strand2 == '+' );
    assert( cig2.strand1 == '+' && cig2.strand2 == '+' );

    int cmpHeader1 = strcmp(cig1->header1, cig2->header1);
    if ( cmpHeader1 == 0 ){ //same query contig/chromosome 
        int cmpStart1 = cig1->start1 - cig2->start1;
        if ( cmpStart1 == 0 ) { //same query start
            int cmpEnd1 = cig1->end1 - cig2->end1;
            if ( cmpEnd1 == 0 ){ //same query end
                int cmpHeader2 = strcmp(cig1->header2, cig2->header2);
                if ( cmpHeader2 == 0 ){
                    int cmpStart2 = cig1->start2 - cig2->start2;
                    if ( cmpStart2 == 0 ){
                        return cig1->end2 - cig2->end2;
                    } else { return cmpStart2; }
                } else { return cmpHeader2; }
            } else { return cmpEnd1; }
        } else { return cmpStart1; }
    } else { return cmpHeader1; }
}

stList * splitCigarsByQueryName( stList *cigars ){
    //Split cigars into sublists of cigars (1 list / queryName)
    stList * cigarsList = stList_construct();

    //Sort the cigars:
    stList_sort(cigars, cigar_cmp);
    st_logInfo("Done sorting cigars. Length: %" PRId64 "\n", stList_length(cigars));

    char *query = NULL;
    stList * currSublist = NULL; //list of cigar of current query
    for ( int64_t i=0; i < stList_length(cigars); i++ ){
        Cigar *cig = stList_get(cigars, i);
        if( query == NULL || strcmp(query, cig->header1) != 0 ){//new query
            if( currSublist != NULL ){//add list of old query to the big list
                stList_append(cigarsList, currSublist);
            }
            query = cig->header1;
            currSublist = stList_construct();
        }
        stList_append(currSublist, cig);
    }
    if( currSublist != NULL && stList_length(currSublist) > 0 ){//add list of old query to the big list
        stList_append(cigarsList, currSublist);
    }
    return cigarsList;
}

stList * readCigar(char *file){
    stList *cigars = stList_construct();
    FILE *fh = fopen(file, "r");
    if (fh == NULL) {
        fprintf(stderr, "Cannot open the cigar file %s.\n", file);
    }

    char header1[50], header2[50], cigarStr[100];
    int64_t start1, end1, start2, end2, score;
    char strand1, strand2;
    while ( fscanf(fh, "cigar: %s %" PRId64 " %" PRId64 " %c %s %" PRId64 " %" PRId64 
                       " %c %" PRId64 " %100[^\n]", header1, &start1, &end1, &strand1, 
                       header2, &start2, &end2, &strand2, &score, cigarStr) == 10 ){
        fscanf(fh, "\n");
        Cigar *cig;
        cig = cigar_construct(header1, start1, end1, strand1, header2, start2, end2, strand2, score, cigarStr);
        stList_append(cigars, cig);
    }
    fclose(fh);
    return cigars;
}

void writeCigar(Cigar *cig, FILE *fh){
    fprintf(fh, "cigar: %s %" PRId64 " %" PRId64 " %c %s %" PRId64 " %" PRId64 " %c %" PRId64 " %s\n", cig->header1, cig->start1, cig->end1, cig->strand1, cig->header2, cig->start2, cig->end2, cig->strand2, cig->score, cig->cigarStr);
}

void writeCigar_long(Cigar *cig, FILE *fh){
    writeCigar(cig, fh); //the original cigar
    fprintf(fh, "\tCorrectBases %" PRId64 ", Total Length: %" PRId64 ", pcCorrect: %f.\n", cig->correctBases, cig->end1 - cig->start1, (float)cig->correctBases/(cig->end1 - cig->start1));
    //write the left-over matched blocks
    for(int64_t i=0; i < stList_length(cig->starts1); i++){
        int64_t s1 = stList_get(cig->starts1, i);
        int64_t e1 = stList_get(cig->ends1, i);
        int64_t s2 = stList_get(cig->starts2, i);
        int64_t e2 = stList_get(cig->ends2, i);
        fprintf(fh, "\t[%" PRId64 ", %" PRId64 ")\t[%" PRId64 ", %" PRId64 ")\n", s1, e1, s2, e2);
    }
}

void writeCigars(stList *cigars, char *outfile){
    FILE *fh = fopen(outfile, "w");
    for (int64_t i=0; i < stList_length(cigars); i++) {
        Cigar *cig = stList_get(cigars, i);
        fprintf(fh, "cigar: %s %" PRId64 " %" PRId64 " %c %s %" PRId64 " %" PRId64 " %c %" PRId64 " %s\n", cig->header1, cig->start1, cig->end1, cig->strand1, cig->header2, cig->start2, cig->end2, cig->strand2, cig->score, cig->cigarStr);
    }
    fclose(fh);
}

/////////////////////////////////////////
///////// CALCULATE TRUE POSITIVES //////
/////////////////////////////////////////
int64_t segment_start(Segment *segment){
    assert(segment != NULL);
    Sequence * sequence = segment_getSequence(segment);
    int64_t start = segment_getStart(segment) - sequence_getStart(sequence);
    if (!segment_getStrand(segment)){//negative strand
        start = (sequence_getStart(sequence) + sequence_getLength(sequence) - 1) - segment_getStart(segment);
    }
    return start;
}

int64_t segment_end(Segment *segment){
    int64_t start = segment_start(segment);
    return start + segment_getLength(segment);
}

bool checkSegmentOverlap(Segment *segment, int64_t start, int64_t end){
    return segment_start(segment) < end  && start < segment_end(segment) ;
}

stList *block_getSegmentByHeader(Block *block, char *header){
    // Return all the segments in the block belong to the sequence matches with the input "header"
    stList *segments = stList_construct();
    Segment *seg = NULL;
    Block_InstanceIterator *blockIt = block_getInstanceIterator(block);
    while ( (seg = block_getNext(blockIt)) != NULL ){
        Sequence *sequence = segment_getSequence(seg);
        if ( strcmp(header, sequence_getHeader(sequence)) == 0 ){
            stList_append(segments, seg);
        }
    }
    block_destructInstanceIterator(blockIt);
    return segments;
}

Segment *checkStrandQuerySegment(Segment *segment, Cigar *cig){
    if (segment != NULL){
        if ( (segment_getStrand(segment) && cig->strand1 == '-') || 
             (!segment_getStrand(segment) && cig->strand1 == '+') ){ //diff strand
             return segment_getReverse(segment);
        }
    }
    return segment;
}

Segment *checkStrandQuerySegment2(Segment *segment, Cigar *cig){
    if (segment != NULL){
        if ( (segment_getStrand(segment) && cig->strand2 == '-') || 
             (!segment_getStrand(segment) && cig->strand2 == '+') ){ //diff strand
             return segment_getReverse(segment);
        }
    }
    return segment;
}

stList * blockGetConstraints(Block *block, stList *constraintsList){
    //Get a list of constraints whose queries ovelap with the block's threads
    stList * blockConstraints = stList_construct();
    for (int64_t i = 0; i< stList_length(constraintsList); i++){
        stList *cigars = stList_get(constraintsList, i);//current sublist of constraint
        if(stList_length(cigars) == 0){
            continue;
        }
        Cigar *firstCig = stList_get(cigars, 0);
        stList *refSegments = block_getSegmentByHeader(block, firstCig->header1);
        
        if ( stList_length(refSegments) == 0 ){ //block does not contain the current query
            continue;
        }
        
        for (int64_t k = 0; k < stList_length(refSegments); k++){//each thread of the block that matches the query
            Segment *refSegment = stList_get(refSegments, k);
            for (int64_t j = 0; j < stList_length(cigars); j++){//find all the constraints that overlap with that thread
                Cigar *cig = stList_get(cigars, j);
                refSegment = checkStrandQuerySegment(refSegment, cig);//make sure Thread has the same strand with query of the current constraint
                if( cig->start1 >= segment_end(refSegment) ){//constraints pass current block, stop
                    break;
                }
                if( checkSegmentOverlap(refSegment, cig->start1, cig->end1) ){
                    if( stList_contains(blockConstraints, cig) == 0 ){//not already added
                        cig->numVisited += 1;
                        stList_append(blockConstraints, cig);
                        /*//DEBUG
                        st_logInfo("\nVisit %d: \n", cig->numVisited);
                        for(int64_t s=0; s < stList_length(cig->starts1); s++){
                            int64_t s1 = stList_get(cig->starts1, s);
                            int64_t e1 = stList_get(cig->ends1, s);
                            int64_t s2 = stList_get(cig->starts2, s);
                            int64_t e2 = stList_get(cig->ends2, s);
                            st_logInfo("\t(%" PRId64 ", %" PRId64 ")\t(%" PRId64 ", %" PRId64 ")\n", s1, e1, s2, e2);
                        }
                        //END DEBUG*/
                    }
                }
            }
            
            /*//DEBUGGING: dbbSEQ10 1 158
            //st_logInfo("Get the set of constraints within current block. Total current constraints: %" PRId64 "\n", stList_length(constraints));
            Sequence *refSeq = segment_getSequence(refSegment);
            char *refheader = sequence_getHeader(refSeq);
            int64_t refstart = -1;
            if ( segment_getStrand(refSegment) ){refstart = segment_start(refSegment);}else{refstart = segment_start(segment_getReverse(refSegment));}
            int64_t refend = refstart + segment_getLength(refSegment);
            if ( strcmp(refheader, "dbbSEQ10") == 0 && refstart < 158 && 1 < refend ){
                st_logInfo("%s\t%" PRId64 "\t%" PRId64 "\n", firstCig->header1, refstart, refend);
                Segment *other = NULL;
                Block_InstanceIterator *blockIt = block_getInstanceIterator(block);
                while ( (other = block_getNext(blockIt)) != NULL ){
                    Sequence *sequence = segment_getSequence(other);
                    char *header = sequence_getHeader(sequence);
                    int64_t ostart = -1;
                    if ( segment_getStrand(other) ){ostart = segment_start(other);}else{ostart = segment_start(segment_getReverse(other));}
                    int64_t oend = ostart + segment_getLength(other);
                    st_logInfo("\t%s\t%" PRId64 "\t%" PRId64 "\n", header, ostart, oend);
                }
                block_destructInstanceIterator(blockIt);
            }
            //END DEBUGGING */
        }
    }
    return blockConstraints;
}

stList * getNonOverlap(int64_t s1, int64_t e1, int64_t s2, int64_t e2){
    stList *positions = stList_construct();
    assert (s1 < e2 && s2 < e1);
    if ( s1 < s2 ){ 
        stList_append(positions, s1);
        stList_append(positions, s2);
    }
    if ( e2 < e1 ){
        stList_append(positions, e2);
        stList_append(positions, e1);
    }
    return positions;
}

stList * getNonOverlap2(int64_t offset, stList *positions1){
    stList *positions2 = stList_construct();
    for(int64_t i=0; i < stList_length(positions1); i++){
        int64_t pos1 = (int64_t) stList_get(positions1, i);
        stList_append(positions2, pos1 + offset);
    }
    return positions2;
}

int64_t getOverlapLength(int64_t s1, int64_t e1, int64_t s2,int64_t e2){
    int64_t start = s1;
    if (s2 > s1){
        start = s2;
    }
    int64_t end = e1;
    if (e2 < e1){
        end = e2;
    }
    return end - start;
}

int cmpCigarAndBlock2(Cigar *cig, Segment *query, Segment *target, int64_t *truePosPointer, FILE *fh){
    //list of part of matched blocks in the cigar that are not overlapped with the alignment
    //return 1 = remove the constraint from the constraint list. Else return 0
    stList *starts1 = stList_construct(); 
    stList *ends1 = stList_construct();
    stList *starts2 = stList_construct();
    stList *ends2 = stList_construct();
    int disagree = 0;
    float minCorrect = 0.98;
    //float minCorrect = 1.0;

    int64_t qstart = segment_start(query);
    int64_t qend = qstart + segment_getLength(query);
    int64_t tstart = segment_start(target);
    int64_t tend = tstart + segment_getLength(target);
    
    while( stList_length(cig->starts1) > 0 ){ //each matched block of the constraint
        int64_t s1 = (int64_t) stList_removeFirst( cig->starts1 );
        int64_t e1 = (int64_t) stList_removeFirst( cig->ends1 );
        int64_t s2 = (int64_t) stList_removeFirst( cig->starts2 );
        int64_t e2 = (int64_t) stList_removeFirst( cig->ends2 );
        //check to see if the alignment agrees with the alignment in the constraint cigar
        if ((s1 < qend && qstart < e1) && (s2 < tend && tstart < e2) &&//both q&t have overlap with block threads
            (s1 - qstart == s2 - tstart)){ //alignment agrees with the constraint cigar
            cig->correctBases += getOverlapLength(s1, e1, qstart, qend);
            stList *positions1 = getNonOverlap(s1, e1, qstart, qend);
            stList *positions2 = getNonOverlap(s2, e2, tstart, tend);
            assert( stList_length(positions1) == stList_length(positions2) );
            for (int i=0; i < stList_length(positions1) - 1; i += 2){
                stList_append(starts1, stList_get(positions1, i));
                stList_append(ends1, stList_get(positions1, i+1));
                stList_append(starts2, stList_get(positions2, i));
                stList_append(ends2, stList_get(positions2, i+1));
            }
        }else{
            stList_append(starts1, s1);
            stList_append(ends1, e1);
            stList_append(starts2, s2);
            stList_append(ends2, e2);
        }
    }

    float correctCoverage = (float)cig->correctBases/(cig->end1 - cig->start1);
    if (correctCoverage >= minCorrect) { //majority of the constraint was kept by cactus
        disagree = 1; //remove constraint
        * truePosPointer += 1; //count as True Positive
    }else {
        if ( stList_length(starts1) > 0 ){
            cig->starts1 = starts1;
            cig->ends1 = ends1;
            cig->starts2 = starts2;
            cig->ends2 = ends2;
        }else{ //fasle positive, constraint was not kept by cactus
            writeCigar_long(cig, fh);
            disagree = 1; //remove constraint
        }
    }
    return disagree;
}

int removeDisagreement(Cigar *cig, stList *segments, FILE *fh){
    for (int64_t j=0; j< stList_length(segments); j++){
        Segment *segment = stList_get(segments, j);
        int64_t qstart = segment_start(segment);
        int64_t qend = qstart + segment_getLength(segment);
        
        stList *starts1 = stList_construct(); 
        stList *ends1 = stList_construct();
        stList *starts2 = stList_construct();
        stList *ends2 = stList_construct();
            
        while( stList_length(cig->starts1) > 0 ){ //each matched block of the constraint
            int64_t s1 = (int64_t) stList_removeFirst( cig->starts1 );
            int64_t e1 = (int64_t) stList_removeFirst( cig->ends1 );
            int64_t s2 = (int64_t) stList_removeFirst( cig->starts2 );
            int64_t e2 = (int64_t) stList_removeFirst( cig->ends2 );
            
            if (s1 < qend && qstart < e1){
                stList *positions1 = getNonOverlap(s1, e1, qstart, qend);
                stList *positions2 = getNonOverlap2(s2 - s1, positions1);
                assert( stList_length(positions1) == stList_length(positions2) );
                for (int64_t i=0; i < stList_length(positions1) - 1; i += 2){
                    stList_append(starts1, stList_get(positions1, i));
                    stList_append(ends1, stList_get(positions1, i+1));
                    stList_append(starts2, stList_get(positions2, i));
                    stList_append(ends2, stList_get(positions2, i+1));
                }
            }else{
                stList_append(starts1, s1);
                stList_append(ends1, e1);
                stList_append(starts2, s2);
                stList_append(ends2, e2);
            }
        }
        if ( stList_length(starts1) > 0 ){
            cig->starts1 = starts1;
            cig->ends1 = ends1;
            cig->starts2 = starts2;
            cig->ends2 = ends2;
        }else{ //all the matched blocks of this constraint were removed
            writeCigar_long(cig, fh);
            return 1;
        }
    }
    
    //Count how many left-over bases
    /*float minCorrect = 0.98;
    int64_t leftoverBases = 0;
    for (int64_t i=0; i < stList_length(cig->starts1); i ++){
        int64_t s1 = stList_get(cig->starts1, i);
        int64_t e1 = stList_get(cig->ends1, i);
        leftoverBases += e1 - s1;
    }
    float maxPossibleCorrect = (float)(leftoverBases + cig->correctBases)/(cig->end1 - cig->start1);
    if( maxPossibleCorrect <  minCorrect ){
        return 1;
    }*/
    return 0;
}

int cmpCigarAndBlock(Cigar *cig, stList *querySegments, stList *targetSegments, int64_t *truePosPointer, FILE *fh){
    int disagree = 0;
    //checking for parts of the constraint that are covered by the alignments in the block
    for(int64_t i=0; i < stList_length(querySegments); i++){
        Segment *query = stList_get(querySegments, i);
        for(int64_t j=0; j < stList_length(targetSegments); j++){
            Segment *target = stList_get(targetSegments, j);
            int status = cmpCigarAndBlock2(cig, query, target, truePosPointer, fh);
            if (status == 1){ //constraint is covered by the alignment, stop checking
                return 1;
            }
        }
    }
    //remove parts of the query contig that were covered in the block
    //but the constraint alignments were not
    disagree = removeDisagreement(cig, querySegments, fh);
    return disagree;
}

void removeConstraint(stList *constraintsList, char *query, Cigar *cigar){
    for (int64_t i =0; i < stList_length(constraintsList); i++){
        stList *constraints = stList_get(constraintsList, i);
        if ( stList_length(constraints) == 0 ){
            continue;
        }else {
            Cigar *firstCig = stList_get(constraints, 0);
            if( strcmp(firstCig->header1, query) == 0 ){
                stList_removeItem(constraints, cigar);
                return;
            }
        }
    }
}

int64_t blockCheckConstraints(Block *block, stList *constraintsList, FILE *fh){
    //returns the number of constraints that the alignment keeps
    int64_t correct = 0;

    //Find list of constraints that overlap with block coordinate range
    stList * blockConstraints = blockGetConstraints( block, constraintsList );
     
    //check the alignment in block to see if it contains the constraints
    for (int64_t i=0; i < stList_length(blockConstraints); i++){
        Cigar * cig = stList_get(blockConstraints, i);
        //First check to see if the block contains both the query and target of the cigar constraint
        stList *segments1 = block_getSegmentByHeader(block, cig->header1);
        stList *segments2 = block_getSegmentByHeader(block, cig->header2);
        
        //make sure segments1 and segments2 are on the same strand with the query and target of the constraint
        stList *newSegments1 = stList_construct();
        for(int64_t j=0; j < stList_length(segments1); j++){
            Segment *segment1 = stList_get(segments1, j);
            stList_append(newSegments1, checkStrandQuerySegment(segment1, cig));
        }

        stList *newSegments2 = stList_construct();
        for(int64_t j=0; j < stList_length(segments2); j++){
            Segment *segment2 = stList_get(segments2, j);
            stList_append(newSegments2, checkStrandQuerySegment2(segment2, cig));
        }

        if ( stList_length(newSegments2) == 0 ){//block contains query but Not target
            if( removeDisagreement(cig, newSegments1, fh) == 1 ){
                removeConstraint(constraintsList, cig->header1, cig);
            }
            continue;
        }

        //Block contains both query & target, now check matched blocks
        if ( cmpCigarAndBlock(cig, newSegments1, newSegments2, &correct, fh) == 1 ){
            removeConstraint(constraintsList, cig->header1, cig);
        }
    }
    return correct; 
}

int64_t checkConstraintTruePos(Flower *flower, stList *constraintsList, FILE *fh){
    //Iterate through each block and check all the constraints overlap with the block
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    int64_t truePos = 0;
    while ( (block = flower_getNextBlock(blockIterator)) != NULL ) {
        truePos += blockCheckConstraints(block, constraintsList, fh);
    }
    flower_destructBlockIterator(blockIterator);

    //Call child flowers recursively
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ( (group = flower_getNextGroup(groupIterator)) != NULL ) {
        if (!group_isLeaf(group)) {
            truePos += checkConstraintTruePos( group_getNestedFlower(group), constraintsList, fh );
        }
    }
    flower_destructGroupIterator(groupIterator);
    return truePos;
}

//////////////////////////////////////////
//////////////// USAGE ///////////////////
//////////////////////////////////////////
void usage(){
    fprintf(stderr, "checkConstraints\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --constraintFile : The constraint cigar file\n");
    fprintf(stderr, "-c --cactusDisk : The location of the cactus disk\n");
    fprintf(stderr, "-e --outputFile : The file to write the output in.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int parseArguments(int64_t argc, char *argv[]){
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    while (1) {
        static struct option long_options[] = { 
            {"logLevel", required_argument, 0, 'a'},
            {"constraintFile", required_argument, 0, 'b'},
            {"cactusDisk", required_argument, 0, 'c'},
            {"outputFile", required_argument, 0, 'e'},
            {0, 0, 0, 0}
        };
        
        int option_index = 0;
        int key = getopt_long(argc, argv, "a:b:c:e:h", long_options, &option_index);
        
        if (key == -1){
            break;
        }
        switch(key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                constraintFile = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'e':
                outputFile = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 1;
        }
    }
    
    //Set up logging
    st_setLogLevelFromString(logLevelString);
    if (outputFile == NULL) {
        st_errAbort("The output file was not specified");
    }
    if (cactusDiskDatabaseString == NULL) {
        st_errAbort("The cactus disk string was not specified");
    }
    if (constraintFile == NULL) {
        st_errAbort("The constraint cigar file was not specified");
    }

    assert(outputFile != NULL);
    assert(cactusDiskDatabaseString != NULL);
    assert(constraintFile != NULL);

    //Log the inputs
    st_logInfo("The cactus disk string: %s\n", cactusDiskDatabaseString);
    st_logInfo("The constraint file: %s\n", constraintFile);
    st_logInfo("The output file: %s\n", outputFile);

    //Load the database
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the cactus disk\n");

    //Clean up
    free(cactusDiskDatabaseString);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    
    //Parse the basic flower
    flower = cactusDisk_getFlower(cactusDisk, 0);
    assert(flower != NULL);
    st_logInfo("Parsed the top level flower of the cactus tree\n");

    return 0;
}

int main(int argc, char *argv[]){
    if ( parseArguments(argc, argv) ){
        return 0;
    }

    //Read in the constraints
    stList *constraints = readCigar(constraintFile);
    int64_t numConstraints = stList_length(constraints);
    if( numConstraints == 0 ){
        return 1;
    }

    st_logInfo("Done reading in the constraints: %" PRId64 " cigars total.\n", numConstraints);
    stList *constraintsList = splitCigarsByQueryName(constraints);
    st_logInfo("Number of constraint lists: %" PRId64 ".\n", stList_length(constraintsList));

    //Check how many constraints the cactus alignment abide by
    FILE *fh = fopen(outputFile, "w");
    int64_t truePos = checkConstraintTruePos(flower, constraintsList, fh); 
    st_logInfo("Done counting the true positive constraints\n");

    //DEBUG
    //fprintf(stdout, "\n****** LEFT-OVER CONSTRAINTS ******\n");
    for(int64_t i = 0; i< stList_length(constraintsList); i++){
        stList *cigars = stList_get(constraintsList, i);
        for(int64_t j = 0; j< stList_length(cigars); j++){
            Cigar *cig = stList_get(cigars, j);
            writeCigar_long(cig, fh);
        }
    }
    //END DEBUG */
    
    float pcTP = (float) truePos *100.0/ (float) numConstraints;
    fprintf(fh, "\n\nNumberOfConstraints\tTruePositives\tPercentTruePositives\n");
    fprintf(fh, "%" PRId64 "\t%" PRId64 "\t%.3f\n", numConstraints, truePos, pcTP);
    fclose(fh);
}

