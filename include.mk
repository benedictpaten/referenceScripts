jobTreeFlags = --batchSystem parasol --logDebug --retryCount 0 --maxThreads 4 --defaultMemory 4294967296
#jobTreeFlags = --batchSystem singleMachine --maxThreads 8 --logDebug --retryCount 0
configFile=${libPath}/config_fast.xml
minimumNsForScaffoldGap=15

#########
#Build basic cactus alignment
#########
		
pipeline :
	rm -rf ./jobTree2
	#Running pipeline to build comparisons
	#python ${binPath}/pipeline.py --maxNumberOfChainsToSolvePerRound '${maxNumberOfChainsToSolvePerRound}' --blastAlignmentStrings '${blastAlignmentStrings}' --baseLevel '${baseLevel}' --maxNumberOfChains '${maxNumberOfChains}' --referenceSpecies '${referenceSpecies}' --haplotypeSequences '${sequences}' --newickTree '${newickTree}' --outputDir ${outputDir} --configFile ${configFile} --referenceAlgorithms '${referenceAlgorithms}' --requiredSpecies '${requiredSequences}' --minimumNsForScaffoldGap ${minimumNsForScaffoldGap} --rangeOfMinimumBlockDegrees '${minimumBlockDegreeRange}' --jobTree ./jobTree ${jobTreeFlags}
	python ${binPath}/../pipeline2.py --maxNumberOfChainsToSolvePerRound '${maxNumberOfChainsToSolvePerRound}' --blastAlignmentStrings '${blastAlignmentStrings}' --baseLevel '${baseLevel}' --maxNumberOfChains '${maxNumberOfChains}' --referenceSpecies '${referenceSpecies}' --haplotypeSequences '${sequences}' --newickTree '${newickTree}' --outputDir ${outputDir} --configFile ${configFile} --referenceAlgorithms '${referenceAlgorithms}' --requiredSpecies '${requiredSequences}' --minimumNsForScaffoldGap ${minimumNsForScaffoldGap} --rangeOfMinimumBlockDegrees '${minimumBlockDegreeRange}' --jobTree ./jobTree2 ${jobTreeFlags}
	jobTreeStatus --jobTree2 ./jobTree2 --failIfNotComplete
	rm -rf ./jobTree2
	
basic : pipeline

#############
#Cleaning stuff
#############

clean :
	rm -rf ${outputDir}/*
