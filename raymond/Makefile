rootPath = ../
experimentName=raymondBasic
binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}dataDir
outputPath = ${rootPath}output

#Tree for the nine haplotypes
#(panTro2:0.007060,(cox:0.000890,((qbl:0.000680,mann:0.000520):0.000270,(((mcf:0.001010,dbb:0.000790):0.000210,apd:0.000740):0.000040,(ssto:0.001000,hg19:0.000800):0.000260):0.000040):0.000100):0.005850);
#Base tree (panTro2:0.007060,x:0.005850);

newickTree=(chimp:0.007060,(00576-1:0.002,01018-2:0.002,03715-1:0.002,04535-2:0.002,10923-1:0.002,14660-2:0.002,14663-1:0.002,mann:0.002,qbl:0.002,00576-2:0.002,01960-1:0.002,03715-2:0.002,10540-1:0.002,10923-2:0.002,14661-1:0.002,14663-2:0.002,cox:0.002,mcf:0.002,01018-1:0.002,01960-2:0.002,04535-1:0.002,10540-2:0.002,14660-1:0.002,14661-2:0.002,dbb:0.002,pgf:0.002,ssto:0.002):0.005850);
dataDir=${dataPath}/raymond
sequences=${dataDir}/chimp ${dataDir}/00576-1 ${dataDir}/01018-2 ${dataDir}/03715-1 ${dataDir}/04535-2 ${dataDir}/10923-1 ${dataDir}/14660-2 ${dataDir}/14663-1 ${dataDir}/mann ${dataDir}/qbl ${dataDir}/00576-2 ${dataDir}/01960-1 ${dataDir}/03715-2 ${dataDir}/10540-1 ${dataDir}/10923-2 ${dataDir}/14661-1 ${dataDir}/14663-2 ${dataDir}/cox ${dataDir}/mcf ${dataDir}/01018-1 ${dataDir}/01960-2 ${dataDir}/04535-1 ${dataDir}/10540-2 ${dataDir}/14660-1 ${dataDir}/14661-2 ${dataDir}/dbb ${dataDir}/pgf ${dataDir}/ssto
singleCopySpecies=00576-1 01018-2 03715-1 04535-2 10923-1 14660-2 14663-1 mann qbl 00576-2 01960-1 03715-2 10540-1 10923-2 14661-1 14663-2 cox mcf 01018-1 01960-2 04535-1 10540-2 14660-1 14661-2 dbb pgf ssto
outputDir=${outputPath}/raymond/${experimentName}
referenceSpecies=reference 00576-1
minimumBlockDegreeRange=1 2
#algorithms: greedy maxCardinality maxWeight blossom5
referenceAlgorithms = maxWeight
baseLevel = 1
maxNumberOfChains = 10000
blastAlignmentStrings = 
#--nogapped %  % --step=10 --seed=match12 --notransition --mismatch=2,100 --match=1,5 --nogapped
theta = 0.00001
useSimulatedAnnealing = 0
permutations = 100
outgroupEvent=chimp
gapGamma = 0.2
constraints=${dataDir}/constraints.cig

include ${rootPath}/include.mk
sampleNumber=100000000

all : basic
