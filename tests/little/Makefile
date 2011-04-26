rootPath = ../../
binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}data
outputPath = ${rootPath}output

newickTree=(qbl:0.000680,mann:0.000520);
dataDir=${dataPath}/mhcHumanVariants
sequences= ${dataDir}/qbl ${dataDir}/mann
requiredSequences=qbl mann
outputDir=${outputPath}/tests/little

minimumBlockDegreeRange=1 2
#2 3
referenceAlgorithms = greedy 
#maxCardinality maxWeight blossom5

include ../../include.mk

all : basic

