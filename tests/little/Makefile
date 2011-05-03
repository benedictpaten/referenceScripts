rootPath = ../../
binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}data
outputPath = ${rootPath}output

newickTree=(CHIMP:0.007060,(HUMAN:0.002000,HUMAN2:0.002000,HUMAN3:0.002000):0.005850);
dataDir=${dataPath}/test
sequences= ${dataDir}/CHIMP ${dataDir}/HUMAN ${dataDir}/HUMAN2 ${dataDir}/HUMAN3
requiredSequences=HUMAN HUMAN2 HUMAN3
outputDir=${outputPath}/tests/little
referenceSpecies=reference HUMAN

minimumBlockDegreeRange=1 2
#2 3
referenceAlgorithms = greedy maxCardinality
#maxCardinality maxWeight blossom5

include ../../include.mk

all : basic

