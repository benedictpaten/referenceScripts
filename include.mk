jobTreeFlags = --batchSystem parasol --logDebug --retryCount 0 --maxThreads 4 --defaultMemory 4294967296
#jobTreeFlags = --batchSystem singleMachine --maxThreads 2 --logDebug --retryCount 0
configFile=${libPath}/config_fast.xml
minimumNsForScaffoldGap=15
sampleNumber=1000000

#########
#Build basic cactus alignment
#########
		
pipeline :
	rm -rf ./jobTree
	#Running pipeline to build comparisons
	python ${binPath}/pipeline.py --heldOutSequences '${heldOutSequences}' --sampleNumber ${sampleNumber} --permutations '${permutations}' --useSimulatedAnnealing '${useSimulatedAnnealing}' --theta '${theta}' --blastAlignmentStrings '${blastAlignmentStrings}' --baseLevel '${baseLevel}' --maxNumberOfChains '${maxNumberOfChains}' --referenceSpecies '${referenceSpecies}' --singleCopySpecies '${singleCopySpecies}' --haplotypeSequences '${sequences}' --newickTree '${newickTree}' --outputDir ${outputDir} --configFile ${configFile} --referenceAlgorithms '${referenceAlgorithms}' --requiredSpecies '${requiredSequences}' --minimumNsForScaffoldGap ${minimumNsForScaffoldGap} --rangeOfMinimumBlockDegrees '${minimumBlockDegreeRange}' --jobTree ./jobTree ${jobTreeFlags}
	jobTreeStatus --jobTree ./jobTree --failIfNotComplete
	rm -rf ./jobTree
	
basic : pipeline

#############
#Cleaning stuff
#############

clean :
	rm -rf ${outputDir}/*
