rootPath = ../../
binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}data
outputPath = ${rootPath}output

newickTree=((HUMAN:0.000680, HUMAN2:0.000680):0.0000001,CHIMP:0.000520);
dataDir=${dataPath}/test
sequences= ${dataDir}/HUMAN ${dataDir}/HUMAN ${dataDir}/CHIMP
#CHIMP
requiredSequences=HUMAN HUMAN2
outputDir=${outputPath}/tests/little
referenceSpecies=reference HUMAN

minimumBlockDegreeRange=1 2
#2 3
referenceAlgorithms = greedy 
#maxCardinality maxWeight blossom5

include ../../include.mk

all : basic

