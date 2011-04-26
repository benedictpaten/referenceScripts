rootPath = ../../
binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}data
outputPath = ${rootPath}output

newickTree=(panTro2:0.007060,(cox:0.000890,((qbl:0.000680,mann:0.000520):0.000270,(((mcf:0.001010,dbb:0.000790):0.000210,apd:0.000740):0.000040,(ssto:0.001000,hg19:0.000800):0.000260):0.000040):0.000100):0.005850);
dataDir=${dataPath}/mhcHumanVariants
sequences= ${dataDir}/panTro2 ${dataDir}/cox ${dataDir}/qbl ${dataDir}/mann ${dataDir}/mcf ${dataDir}/dbb ${dataDir}/apd ${dataDir}/ssto ${dataDir}/hg19
requiredSequences=cox qbl mann mcf dbb apd ssto hg19
outputDir=${outputPath}/main/mhcHumanVariants

minimumBlockDegreeRange=1 2 3
referenceAlgorithms = greedy maxCardinality maxWeight blossom5

include ../../include.mk

all : basic

