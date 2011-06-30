rootPath = ../../
binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}data
outputPath = ${rootPath}output

newickTree=(panTro2:0.007060,(hg19:0.002000,mann:0.002000,NA12878:0.002000):0.005850);
dataDir=${dataPath}/mhcHumanVariantsNsRemoved
sequences= ${dataDir}/panTro2 ${dataDir}/hg19 ${dataDir}/mann ${dataDir}/NA12878.fa
requiredSequences=hg19 mann NA12878
outputDir=${outputPath}/tests/big
referenceSpecies=reference hg19
minimumBlockDegreeRange=1 2 3
#algorithms: greedy maxCardinality maxWeight blossom5
referenceAlgorithms = maxCardinality maxWeight
baseLevel = 1
maxNumberOfChains = 0
blastAlignmentStrings = --nogapped %  % --step=10 --seed=match12 --notransition --mismatch=2,100 --match=1,5 --nogapped
maxNumberOfChainsToSolvePerRound = 10 1000

include ../../include.mk

all : basic

