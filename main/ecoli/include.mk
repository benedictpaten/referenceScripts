binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}dataDir
outputPath = ${rootPath}output

#Tree for the nine haplotypes
#(panTro2:0.007060,(cox:0.000890,((qbl:0.000680,mann:0.000520):0.000270,(((mcf:0.001010,dbb:0.000790):0.000210,apd:0.000740):0.000040,(ssto:0.001000,hg19:0.000800):0.000260):0.000040):0.000100):0.005850);
#Base tree (panTro2:0.007060,x:0.005850);

newickTree=(Escherichia_coli_042_uid161985:0.0001,Escherichia_coli_536_uid58531:0.0001,Escherichia_coli_55989_uid59383:0.0001,Escherichia_coli_ABU_83972_uid161975:0.0001,Escherichia_coli_APEC_O1_uid58623:0.0001,Escherichia_coli_ATCC_8739_uid58783:0.0001,Escherichia_coli_BL21_DE3__uid161947:0.0001,Escherichia_coli_BL21_DE3__uid161949:0.0001,Escherichia_coli_BW2952_uid59391:0.0001,Escherichia_coli_B_REL606_uid58803:0.0001,Escherichia_coli_CFT073_uid57915:0.0001,Escherichia_coli_DH1_uid161951:0.0001,Escherichia_coli_DH1_uid162051:0.0001,Escherichia_coli_E24377A_uid58395:0.0001,Escherichia_coli_ED1a_uid59379:0.0001,Escherichia_coli_ETEC_H10407_uid161993:0.0001,Escherichia_coli_HS_uid58393:0.0001,Escherichia_coli_IAI1_uid59377:0.0001,Escherichia_coli_IAI39_uid59381:0.0001,Escherichia_coli_IHE3034_uid162007:0.0001,Escherichia_coli_KO11FL_uid162099:0.0001,Escherichia_coli_KO11FL_uid52593:0.0001,Escherichia_coli_K_12_substr__DH10B_uid58979:0.0001,Escherichia_coli_K_12_substr__MG1655_uid57779:0.0001,Escherichia_coli_K_12_substr__W3110_uid161931:0.0001,Escherichia_coli_LF82_uid161965:0.0001,Escherichia_coli_NA114_uid162139:0.0001,Escherichia_coli_O103_H2_12009_uid41013:0.0001,Escherichia_coli_O104_H4_2009EL_2050_uid175905:0.0001,Escherichia_coli_O104_H4_2009EL_2071_uid176128:0.0001,Escherichia_coli_O104_H4_2011C_3493_uid176127:0.0001,Escherichia_coli_O111_H__11128_uid41023:0.0001,Escherichia_coli_O127_H6_E2348_69_uid59343:0.0001,Escherichia_coli_O157_H7_EC4115_uid59091:0.0001,Escherichia_coli_O157_H7_EDL933_uid57831:0.0001,Escherichia_coli_O157_H7_Sakai_uid57781:0.0001,Escherichia_coli_O157_H7_TW14359_uid59235:0.0001,Escherichia_coli_O26_H11_11368_uid41021:0.0001,Escherichia_coli_O55_H7_CB9615_uid46655:0.0001,Escherichia_coli_O55_H7_RM12579_uid162153:0.0001,Escherichia_coli_O7_K1_CE10_uid162115:0.0001,Escherichia_coli_O83_H1_NRG_857C_uid161987:0.0001,Escherichia_coli_P12b_uid162061:0.0001,Escherichia_coli_S88_uid62979:0.0001,Escherichia_coli_SE11_uid59425:0.0001,Escherichia_coli_SE15_uid161939:0.0001,Escherichia_coli_SMS_3_5_uid58919:0.0001,Escherichia_coli_UM146_uid162043:0.0001,Escherichia_coli_UMN026_uid62981:0.0001,Escherichia_coli_UMNK88_uid161991:0.0001,Escherichia_coli_UTI89_uid58541:0.0001,Escherichia_coli_W_uid162011:0.0001,Escherichia_coli_W_uid162101:0.0001,Escherichia_coli_Xuzhou21_uid163995:0.0001,Escherichia_coli__BL21_Gold_DE3_pLysS_AG__uid59245:0.0001,Escherichia_coli__clone_D_i14__uid162049:0.0001,Escherichia_coli__clone_D_i2__uid162047:0.0001,Shigella_boydii_CDC_3083_94_uid58415:0.0001,Shigella_boydii_Sb227_uid58215:0.0001,Shigella_dysenteriae_Sd197_uid58213:0.0001,Shigella_flexneri_2002017_uid159233:0.0001,Shigella_flexneri_2a_2457T_uid57991:0.0001,Shigella_flexneri_2a_301_uid62907:0.0001,Shigella_flexneri_5_8401_uid58583:0.0001,Shigella_sonnei_53G_uid84383:0.0001,Shigella_sonnei_Ss046_uid58217)reference;
dataDir=${dataPath}/${experimentName}
sequences=${dataDir}/Escherichia_coli_042_uid161985 ${dataDir}/Escherichia_coli_536_uid58531 ${dataDir}/Escherichia_coli_55989_uid59383 ${dataDir}/Escherichia_coli_ABU_83972_uid161975 ${dataDir}/Escherichia_coli_APEC_O1_uid58623 ${dataDir}/Escherichia_coli_ATCC_8739_uid58783 ${dataDir}/Escherichia_coli_BL21_DE3__uid161947 ${dataDir}/Escherichia_coli_BL21_DE3__uid161949 ${dataDir}/Escherichia_coli_BW2952_uid59391 ${dataDir}/Escherichia_coli_B_REL606_uid58803 ${dataDir}/Escherichia_coli_CFT073_uid57915 ${dataDir}/Escherichia_coli_DH1_uid161951 ${dataDir}/Escherichia_coli_DH1_uid162051 ${dataDir}/Escherichia_coli_E24377A_uid58395 ${dataDir}/Escherichia_coli_ED1a_uid59379 ${dataDir}/Escherichia_coli_ETEC_H10407_uid161993 ${dataDir}/Escherichia_coli_HS_uid58393 ${dataDir}/Escherichia_coli_IAI1_uid59377 ${dataDir}/Escherichia_coli_IAI39_uid59381 ${dataDir}/Escherichia_coli_IHE3034_uid162007 ${dataDir}/Escherichia_coli_KO11FL_uid162099 ${dataDir}/Escherichia_coli_KO11FL_uid52593 ${dataDir}/Escherichia_coli_K_12_substr__DH10B_uid58979 ${dataDir}/Escherichia_coli_K_12_substr__MG1655_uid57779 ${dataDir}/Escherichia_coli_K_12_substr__W3110_uid161931 ${dataDir}/Escherichia_coli_LF82_uid161965 ${dataDir}/Escherichia_coli_NA114_uid162139 ${dataDir}/Escherichia_coli_O103_H2_12009_uid41013 ${dataDir}/Escherichia_coli_O104_H4_2009EL_2050_uid175905 ${dataDir}/Escherichia_coli_O104_H4_2009EL_2071_uid176128 ${dataDir}/Escherichia_coli_O104_H4_2011C_3493_uid176127 ${dataDir}/Escherichia_coli_O111_H__11128_uid41023 ${dataDir}/Escherichia_coli_O127_H6_E2348_69_uid59343 ${dataDir}/Escherichia_coli_O157_H7_EC4115_uid59091 ${dataDir}/Escherichia_coli_O157_H7_EDL933_uid57831 ${dataDir}/Escherichia_coli_O157_H7_Sakai_uid57781 ${dataDir}/Escherichia_coli_O157_H7_TW14359_uid59235 ${dataDir}/Escherichia_coli_O26_H11_11368_uid41021 ${dataDir}/Escherichia_coli_O55_H7_CB9615_uid46655 ${dataDir}/Escherichia_coli_O55_H7_RM12579_uid162153 ${dataDir}/Escherichia_coli_O7_K1_CE10_uid162115 ${dataDir}/Escherichia_coli_O83_H1_NRG_857C_uid161987 ${dataDir}/Escherichia_coli_P12b_uid162061 ${dataDir}/Escherichia_coli_S88_uid62979 ${dataDir}/Escherichia_coli_SE11_uid59425 ${dataDir}/Escherichia_coli_SE15_uid161939 ${dataDir}/Escherichia_coli_SMS_3_5_uid58919 ${dataDir}/Escherichia_coli_UM146_uid162043 ${dataDir}/Escherichia_coli_UMN026_uid62981 ${dataDir}/Escherichia_coli_UMNK88_uid161991 ${dataDir}/Escherichia_coli_UTI89_uid58541 ${dataDir}/Escherichia_coli_W_uid162011 ${dataDir}/Escherichia_coli_W_uid162101 ${dataDir}/Escherichia_coli_Xuzhou21_uid163995 ${dataDir}/Escherichia_coli__BL21_Gold_DE3_pLysS_AG__uid59245 ${dataDir}/Escherichia_coli__clone_D_i14__uid162049 ${dataDir}/Escherichia_coli__clone_D_i2__uid162047 ${dataDir}/Shigella_boydii_CDC_3083_94_uid58415 ${dataDir}/Shigella_boydii_Sb227_uid58215 ${dataDir}/Shigella_dysenteriae_Sd197_uid58213 ${dataDir}/Shigella_flexneri_2002017_uid159233 ${dataDir}/Shigella_flexneri_2a_2457T_uid57991 ${dataDir}/Shigella_flexneri_2a_301_uid62907 ${dataDir}/Shigella_flexneri_5_8401_uid58583 ${dataDir}/Shigella_sonnei_53G_uid84383 ${dataDir}/Shigella_sonnei_Ss046_uid58217
outputDir=${outputPath}/main/${experimentName}
referenceSpecies=reference Escherichia_coli_K_12_substr__DH10B_uid58979

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
outgroupEvent=
gapGamma = 0.0
#gapGamma = 0.2
constraints=${dataDir}/constraints.cig

include ${rootPath}/include.mk

sampleNumber=100000000



