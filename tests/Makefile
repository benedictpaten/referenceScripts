rootPath = ../
binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}dataDir
outputPath = ${rootPath}output

#newickTree=(CHIMP:0.007060,(HUMAN:0.002000,HUMAN2:0.002000,HUMAN3:0.002000,HUMAN4:0.002000):0.005850);
newickTree=(CHIMP:0.01291,HUMAN:0.002000,HUMAN2:0.002000,HUMAN3:0.002000,HUMAN4:0.002000)reference;
dataDir=${dataPath}/test
sequences= ${dataDir}/CHIMP ${dataDir}/HUMAN ${dataDir}/HUMAN2 ${dataDir}/HUMAN3 ${dataDir}/HUMAN4
outputDir=${outputPath}/tests
referenceSpecies=reference HUMAN
minimumBlockDegreeRange=2
referenceAlgorithms = maxCardinality
baseLevel = 1
maxNumberOfChains = 0
#blastAlignmentStrings = --nogapped % --step=10 --seed=match12 --notransition --mismatch=2,100 --match=1,5 --nogapped
blastAlignmentStrings = 
#% --step=10 --seed=match12 --notransition --mismatch=2,100 --match=1,5 --nogapped
theta = 0.001
useSimulatedAnnealing = 0
permutations = 500
#heldOutSequences=CHIMP
outgroupEvent=HUMAN3
gapGamma = 0.0 0.2
constraints=${dataDir}/constraints.cig
singleCopyIngroup=1 0

#maxCardinality maxWeight blossom5

include ../include.mk

all : basic

