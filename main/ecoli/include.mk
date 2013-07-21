binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}dataDir
outputPath = ${rootPath}output

#Tree for the nine haplotypes
#(panTro2:0.007060,(cox:0.000890,((qbl:0.000680,mann:0.000520):0.000270,(((mcf:0.001010,dbb:0.000790):0.000210,apd:0.000740):0.000040,(ssto:0.001000,hg19:0.000800):0.000260):0.000040):0.000100):0.005850);
#Base tree (panTro2:0.007060,x:0.005850);

newickTree=(EscherichiaColi042Uid161985:0.0001,EscherichiaColi536Uid58531:0.0001,EscherichiaColi55989Uid59383:0.0001,EscherichiaColiAbu83972Uid161975:0.0001,EscherichiaColiApecO1Uid58623:0.0001,EscherichiaColiAtcc8739Uid58783:0.0001,EscherichiaColiBRel606Uid58803:0.0001,EscherichiaColiBl21De3Uid161947:0.0001,EscherichiaColiBl21De3Uid161949:0.0001,EscherichiaColiBl21GoldDe3PlyssAgUid59245:0.0001,EscherichiaColiBw2952Uid59391:0.0001,EscherichiaColiCft073Uid57915:0.0001,EscherichiaColiCloneDI14Uid162049:0.0001,EscherichiaColiCloneDI2Uid162047:0.0001,EscherichiaColiDh1Uid161951:0.0001,EscherichiaColiDh1Uid162051:0.0001,EscherichiaColiE24377aUid58395:0.0001,EscherichiaColiEd1aUid59379:0.0001,EscherichiaColiEtecH10407Uid161993:0.0001,EscherichiaColiHsUid58393:0.0001,EscherichiaColiIai1Uid59377:0.0001,EscherichiaColiIai39Uid59381:0.0001,EscherichiaColiIhe3034Uid162007:0.0001,EscherichiaColiK12SubstrDh10bUid58979:0.0001,EscherichiaColiK12SubstrMg1655Uid57779:0.0001,EscherichiaColiK12SubstrW3110Uid161931:0.0001,EscherichiaColiKo11flUid162099:0.0001,EscherichiaColiKo11flUid52593:0.0001,EscherichiaColiLf82Uid161965:0.0001,EscherichiaColiNa114Uid162139:0.0001,EscherichiaColiO103H212009Uid41013:0.0001,EscherichiaColiO104H42009el2050Uid175905:0.0001,EscherichiaColiO104H42009el2071Uid176128:0.0001,EscherichiaColiO104H42011c3493Uid176127:0.0001,EscherichiaColiO111H11128Uid41023:0.0001,EscherichiaColiO127H6E234869Uid59343:0.0001,EscherichiaColiO157H7Ec4115Uid59091:0.0001,EscherichiaColiO157H7Edl933Uid57831:0.0001,EscherichiaColiO157H7SakaiUid57781:0.0001,EscherichiaColiO157H7Tw14359Uid59235:0.0001,EscherichiaColiO26H1111368Uid41021:0.0001,EscherichiaColiO55H7Cb9615Uid46655:0.0001,EscherichiaColiO55H7Rm12579Uid162153:0.0001,EscherichiaColiO7K1Ce10Uid162115:0.0001,EscherichiaColiO83H1Nrg857cUid161987:0.0001,EscherichiaColiP12bUid162061:0.0001,EscherichiaColiS88Uid62979:0.0001,EscherichiaColiSe11Uid59425:0.0001,EscherichiaColiSe15Uid161939:0.0001,EscherichiaColiSms35Uid58919:0.0001,EscherichiaColiUm146Uid162043:0.0001,EscherichiaColiUmn026Uid62981:0.0001,EscherichiaColiUmnk88Uid161991:0.0001,EscherichiaColiUti89Uid58541:0.0001,EscherichiaColiWUid162011:0.0001,EscherichiaColiWUid162101:0.0001,EscherichiaColiXuzhou21Uid163995:0.0001,ShigellaBoydiiCdc308394Uid58415:0.0001,ShigellaBoydiiSb227Uid58215:0.0001,ShigellaDysenteriaeSd197Uid58213:0.0001,ShigellaFlexneri2002017Uid159233:0.0001,ShigellaFlexneri2a2457tUid57991:0.0001,ShigellaFlexneri2a301Uid62907:0.0001,ShigellaFlexneri58401Uid58583:0.0001,ShigellaSonnei53gUid84383:0.0001,ShigellaSonneiSs046Uid58217:0.0001)reference;
dataDir=${dataPath}/${experimentName}
sequences=${dataDir}/EscherichiaColi042Uid161985 ${dataDir}/EscherichiaColi536Uid58531 ${dataDir}/EscherichiaColi55989Uid59383 ${dataDir}/EscherichiaColiAbu83972Uid161975 ${dataDir}/EscherichiaColiApecO1Uid58623 ${dataDir}/EscherichiaColiAtcc8739Uid58783 ${dataDir}/EscherichiaColiBRel606Uid58803 ${dataDir}/EscherichiaColiBl21De3Uid161947 ${dataDir}/EscherichiaColiBl21De3Uid161949 ${dataDir}/EscherichiaColiBl21GoldDe3PlyssAgUid59245 ${dataDir}/EscherichiaColiBw2952Uid59391 ${dataDir}/EscherichiaColiCft073Uid57915 ${dataDir}/EscherichiaColiCloneDI14Uid162049 ${dataDir}/EscherichiaColiCloneDI2Uid162047 ${dataDir}/EscherichiaColiDh1Uid161951 ${dataDir}/EscherichiaColiDh1Uid162051 ${dataDir}/EscherichiaColiE24377aUid58395 ${dataDir}/EscherichiaColiEd1aUid59379 ${dataDir}/EscherichiaColiEtecH10407Uid161993 ${dataDir}/EscherichiaColiHsUid58393 ${dataDir}/EscherichiaColiIai1Uid59377 ${dataDir}/EscherichiaColiIai39Uid59381 ${dataDir}/EscherichiaColiIhe3034Uid162007 ${dataDir}/EscherichiaColiK12SubstrDh10bUid58979 ${dataDir}/EscherichiaColiK12SubstrMg1655Uid57779 ${dataDir}/EscherichiaColiK12SubstrW3110Uid161931 ${dataDir}/EscherichiaColiKo11flUid162099 ${dataDir}/EscherichiaColiKo11flUid52593 ${dataDir}/EscherichiaColiLf82Uid161965 ${dataDir}/EscherichiaColiNa114Uid162139 ${dataDir}/EscherichiaColiO103H212009Uid41013 ${dataDir}/EscherichiaColiO104H42009el2050Uid175905 ${dataDir}/EscherichiaColiO104H42009el2071Uid176128 ${dataDir}/EscherichiaColiO104H42011c3493Uid176127 ${dataDir}/EscherichiaColiO111H11128Uid41023 ${dataDir}/EscherichiaColiO127H6E234869Uid59343 ${dataDir}/EscherichiaColiO157H7Ec4115Uid59091 ${dataDir}/EscherichiaColiO157H7Edl933Uid57831 ${dataDir}/EscherichiaColiO157H7SakaiUid57781 ${dataDir}/EscherichiaColiO157H7Tw14359Uid59235 ${dataDir}/EscherichiaColiO26H1111368Uid41021 ${dataDir}/EscherichiaColiO55H7Cb9615Uid46655 ${dataDir}/EscherichiaColiO55H7Rm12579Uid162153 ${dataDir}/EscherichiaColiO7K1Ce10Uid162115 ${dataDir}/EscherichiaColiO83H1Nrg857cUid161987 ${dataDir}/EscherichiaColiP12bUid162061 ${dataDir}/EscherichiaColiS88Uid62979 ${dataDir}/EscherichiaColiSe11Uid59425 ${dataDir}/EscherichiaColiSe15Uid161939 ${dataDir}/EscherichiaColiSms35Uid58919 ${dataDir}/EscherichiaColiUm146Uid162043 ${dataDir}/EscherichiaColiUmn026Uid62981 ${dataDir}/EscherichiaColiUmnk88Uid161991 ${dataDir}/EscherichiaColiUti89Uid58541 ${dataDir}/EscherichiaColiWUid162011 ${dataDir}/EscherichiaColiWUid162101 ${dataDir}/EscherichiaColiXuzhou21Uid163995 ${dataDir}/ShigellaBoydiiCdc308394Uid58415 ${dataDir}/ShigellaBoydiiSb227Uid58215 ${dataDir}/ShigellaDysenteriaeSd197Uid58213 ${dataDir}/ShigellaFlexneri2002017Uid159233 ${dataDir}/ShigellaFlexneri2a2457tUid57991 ${dataDir}/ShigellaFlexneri2a301Uid62907 ${dataDir}/ShigellaFlexneri58401Uid58583 ${dataDir}/ShigellaSonnei53gUid84383 ${dataDir}/ShigellaSonneiSs046Uid58217
outputDir=${outputPath}/main/${experimentName}
referenceSpecies=reference EscherichiacoiK12substrDH10Buid58979

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
outgroupEvent=ShigellaSonneiSs046Uid58217
gapGamma = 0.0
#gapGamma = 0.2
constraints=${dataDir}/constraints.cig
singleCopyIngroup=1 0

include ${rootPath}/include.mk

sampleNumber=1000000



