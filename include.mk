#jobTreeFlags = --batchSystem parasol --logDebug --retryCount 0 --maxThreads 4
jobTreeFlags = --batchSystem singleMachine --maxThreads 1 --logDebug --retryCount 0
configFile=${libPath}/config_fast.xml
minimumNsForScaffoldGap=15

#########
#Build basic cactus alignment
#########
		
pipeline :
	rm -rf ./jobTree
	#Running pipeline to build comparisons
	python ${binPath}/pipeline.py --referenceSpecies '${referenceSpecies}' --haplotypeSequences '${sequences}' --newickTree '${newickTree}' --outputDir ${outputDir} --configFile ${configFile} --referenceAlgorithms '${referenceAlgorithms}' --requiredSpecies '${requiredSequences}' --minimumNsForScaffoldGap ${minimumNsForScaffoldGap} --rangeOfMinimumBlockDegrees '${minimumBlockDegreeRange}' --jobTree ./jobTree ${jobTreeFlags}
	jobTreeStatus --jobTree ./jobTree --failIfNotComplete
	rm -rf ./jobTree
	
basic : pipeline

#############
#Cleaning stuff
#############

clean :
	rm -rf ${outputDir}/*