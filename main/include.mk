binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}dataDir
outputPath = ${rootPath}output

#Tree for the nine haplotypes
#(panTro2:0.007060,(cox:0.000890,((qbl:0.000680,mann:0.000520):0.000270,(((mcf:0.001010,dbb:0.000790):0.000210,apd:0.000740):0.000040,(ssto:0.001000,hg19:0.000800):0.000260):0.000040):0.000100):0.005850);
#Base tree (panTro2:0.007060,x:0.005850);

newickTree=(panTro3:0.01291,NA12878:0.002000,NA12892:0.002000,NA19239:0.002000,apd:0.002000,dbb:0.002000,mann:0.002000,nigerian:0.002000,qbl:0.002000,venter:0.002000,yanhuang:0.002000,NA19238:0.002000,NA19240:0.002000,cox:0.002000,hg19:0.002000,mcf:0.002000,ssto:0.002000,cgIndelContigs.fa:0.002,1kIndelContigs.fa:0.002,snp135IndelContigs.fa:0.002)reference;
dataDir=${dataPath}/${experimentName}
sequences=${dataDir}/panTro3 ${dataDir}/NA12878 ${dataDir}/NA12892 ${dataDir}/NA19239 ${dataDir}/apd ${dataDir}/dbb ${dataDir}/mann ${dataDir}/nigerian ${dataDir}/qbl ${dataDir}/venter ${dataDir}/yanhuang ${dataDir}/NA19238 ${dataDir}/NA19240 ${dataDir}/cox ${dataDir}/hg19 ${dataDir}/mcf ${dataDir}/ssto ${dataDir}/cgIndelContigs.fa ${dataDir}/1kIndelContigs.fa ${dataDir}/snp135IndelContigs.fa
outputDir=${outputPath}/main/${experimentName}
referenceSpecies=reference hg19
minimumBlockDegreeRange=2 
#3 4 5 6 7 8
#algorithms: greedy maxCardinality maxWeight blossom5
referenceAlgorithms = maxWeight
baseLevel = 1
maxNumberOfChains = 10000
blastAlignmentStrings = 
#--nogapped %  % --step=10 --seed=match12 --notransition --mismatch=2,100 --match=1,5 --nogapped
theta = 0.00001
#theta = 0.0 0.0000001 0.000001 0.00001 0.01
useSimulatedAnnealing = 0
permutations = 100
outgroupEvent=panTro3
gapGamma = 0.0
#gapGamma = 0.2
constraints=${dataDir}/constraints.cig

include ${rootPath}/include.mk

sampleNumber=100000000



