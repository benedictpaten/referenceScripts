rootPath = ../
experimentName=raymondBasic
binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}dataDir
outputPath = ${rootPath}output

#Tree for the nine haplotypes
#(panTro2:0.007060,(cox:0.000890,((qbl:0.000680,mann:0.000520):0.000270,(((mcf:0.001010,dbb:0.000790):0.000210,apd:0.000740):0.000040,(ssto:0.001000,hg19:0.000800):0.000260):0.000040):0.000100):0.005850);
#Base tree (panTro2:0.007060,x:0.005850);

newickTree=(AY663401.fa:0.007060,(AY663393.fa:0.002,AY663396.fa:0.002,AY663399.fa:0.002,AY663404.fa:0.002,AY663407.fa:0.002,AY663410.fa:0.002,AY663413.fa:0.002,AY663394.fa:0.002,AY663397.fa:0.002,AY663400.fa:0.002,AY663405.fa:0.002,AY663408.fa:0.002,AY663411.fa:0.002,AY663414.fa:0.002,AY663395.fa:0.002,AY663398.fa:0.002,AY663406.fa:0.002,AY663409.fa:0.002,AY663412.fa:0.002,AY663415.fa:0.002):0.005850);
dataDir=${dataPath}/raymond
sequences=${dataDir}/AY663401.fa ${dataDir}/AY663393.fa	${dataDir}/AY663396.fa ${dataDir}/AY663399.fa ${dataDir}/AY663404.fa ${dataDir}/AY663407.fa ${dataDir}/AY663410.fa ${dataDir}/AY663413.fa ${dataDir}/AY663394.fa ${dataDir}/AY663397.fa ${dataDir}/AY663400.fa ${dataDir}/AY663405.fa ${dataDir}/AY663408.fa ${dataDir}/AY663411.fa ${dataDir}/AY663414.fa ${dataDir}/AY663395.fa ${dataDir}/AY663398.fa ${dataDir}/AY663406.fa ${dataDir}/AY663409.fa ${dataDir}/AY663412.fa ${dataDir}/AY663415.fa
singleCopySpecies=AY663393.fa AY663396.fa AY663399.fa AY663404.fa AY663407.fa AY663410.fa AY663413.fa AY663394.fa AY663397.fa AY663400.fa AY663405.fa AY663408.fa AY663411.fa AY663414.fa AY663395.fa AY663398.fa AY663406.fa AY663409.fa AY663412.fa AY663415.fa
outputDir=${outputPath}/raymond/${experimentName}
referenceSpecies=reference AY663393.fa
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
outgroupEvent=AY663401.fa
gapGamma = 0.2
constraints=${dataDir}/constraints.cig

include ${rootPath}/include.mk
sampleNumber=100000000

all : basic
