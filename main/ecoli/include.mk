binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}dataDir
outputPath = ${rootPath}output

#Tree for the nine haplotypes
#(panTro2:0.007060,(cox:0.000890,((qbl:0.000680,mann:0.000520):0.000270,(((mcf:0.001010,dbb:0.000790):0.000210,apd:0.000740):0.000040,(ssto:0.001000,hg19:0.000800):0.000260):0.000040):0.000100):0.005850);
#Base tree (panTro2:0.007060,x:0.005850);

newickTree=(Escherichiacoli042uid161985:0.0001,Escherichiacoli536uid58531:0.0001,Escherichiacoli55989uid59383:0.0001,EscherichiacoliABU83972uid161975:0.0001,EscherichiacoliAPECO1uid58623:0.0001,EscherichiacoliATCC8739uid58783:0.0001,EscherichiacoliBL21DE3uid161947:0.0001,EscherichiacoliBL21DE3uid161949:0.0001,EscherichiacoliBW2952uid59391:0.0001,EscherichiacoliBREL606uid58803:0.0001,EscherichiacoliCFT073uid57915:0.0001,EscherichiacoliDH1uid161951:0.0001,EscherichiacoliDH1uid162051:0.0001,EscherichiacoliE24377Auid58395:0.0001,EscherichiacoliED1auid59379:0.0001,EscherichiacoliETECH10407uid161993:0.0001,EscherichiacoliHSuid58393:0.0001,EscherichiacoliIAI1uid59377:0.0001,EscherichiacoliIAI39uid59381:0.0001,EscherichiacoliIHE3034uid162007:0.0001,EscherichiacoliKO11FLuid162099:0.0001,EscherichiacoliKO11FLuid52593:0.0001,EscherichiacoliK12substrDH10Buid58979:0.0001,EscherichiacoliK12substrMG1655uid57779:0.0001,EscherichiacoliK12substrW3110uid161931:0.0001,EscherichiacoliLF82uid161965:0.0001,EscherichiacoliNA114uid162139:0.0001,EscherichiacoliO103H212009uid41013:0.0001,EscherichiacoliO104H42009EL2050uid175905:0.0001,EscherichiacoliO104H42009EL2071uid176128:0.0001,EscherichiacoliO104H42011C3493uid176127:0.0001,EscherichiacoliO111H11128uid41023:0.0001,EscherichiacoliO127H6E234869uid59343:0.0001,EscherichiacoliO157H7EC4115uid59091:0.0001,EscherichiacoliO157H7EDL933uid57831:0.0001,EscherichiacoliO157H7Sakaiuid57781:0.0001,EscherichiacoliO157H7TW14359uid59235:0.0001,EscherichiacoliO26H1111368uid41021:0.0001,EscherichiacoliO55H7CB9615uid46655:0.0001,EscherichiacoliO55H7RM12579uid162153:0.0001,EscherichiacoliO7K1CE10uid162115:0.0001,EscherichiacoliO83H1NRG857Cuid161987:0.0001,EscherichiacoliP12buid162061:0.0001,EscherichiacoliS88uid62979:0.0001,EscherichiacoliSE11uid59425:0.0001,EscherichiacoliSE15uid161939:0.0001,EscherichiacoliSMS35uid58919:0.0001,EscherichiacoliUM146uid162043:0.0001,EscherichiacoliUMN026uid62981:0.0001,EscherichiacoliUMNK88uid161991:0.0001,EscherichiacoliUTI89uid58541:0.0001,EscherichiacoliWuid162011:0.0001,EscherichiacoliWuid162101:0.0001,EscherichiacoliXuzhou21uid163995:0.0001,EscherichiacoliBL21GoldDE3pLysSAGuid59245:0.0001,EscherichiacolicloneDi14uid162049:0.0001,EscherichiacolicloneDi2uid162047:0.0001,ShigellaboydiiCDC308394uid58415:0.0001,ShigellaboydiiSb227uid58215:0.0001,ShigelladysenteriaeSd197uid58213:0.0001,Shigellaflexneri2002017uid159233:0.0001,Shigellaflexneri2a2457Tuid57991:0.0001,Shigellaflexneri2a301uid62907:0.0001,Shigellaflexneri58401uid58583:0.0001,Shigellasonnei53Guid84383:0.0001,ShigellasonneiSs046uid58217:0.0001)reference;
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
outgroupEvent=ShigellasonneiSs046uid58217
gapGamma = 0.0
#gapGamma = 0.2
constraints=${dataDir}/constraints.cig
singleCopyIngroup=1 0

include ${rootPath}/include.mk

sampleNumber=1000000



