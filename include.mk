#jobTreeFlags = --batchSystem parasol --logDebug --retryCount 0 --maxThreads 4 --defaultMemory 4294967296
jobTreeFlags = --batchSystem singleMachine --maxThreads 20 --logDebug --retryCount 0 --maxCpus 20
configFile=${libPath}/config_fast.xml
minimumNsForScaffoldGap=15
sampleNumber=1000000

batchSystemForAlignments = singleMachine
databaseHost = localhost
parasolCommandForAlignment = parasol
extraJobTreeArgs = --maxMemory 746870912000 --maxCpus 25 --maxThreads 25

#batchSystemForAlignments = parasol
#databaseHost = kolossus-10
#parasolCommandForAlignment = /hive/users/benedict/parasol -host=swarm-10
#batchSystemForAlignments = parasol
#extraJobTreeArgs = --bigBatchSystem singleMachine --bigMemoryThreshold 2147483648 --bigMaxMemory 746870912000 --bigCpuThreshold 4 --bigMaxCpus 25 --maxThreads 25

#########
#Build basic cactus alignment
#########
		
pipeline :
	rm -rf ./jobTree
	#Running pipeline to build comparisons
	python ${binPath}/pipeline.py --defaultMemory 8589934593 --databaseHost='${databaseHost}' --extraJobTreeArgs '${extraJobTreeArgs}' --batchSystemForAlignments='${batchSystemForAlignments}' --parasolCommandForAlignment='${parasolCommandForAlignment}' --singleCopyIngroup '${singleCopyIngroup}' --gapGamma '${gapGamma}' --constraints ${constraints} --outgroupEvent ${outgroupEvent} --heldOutSequences '${heldOutSequences}' --sampleNumber ${sampleNumber} --permutations '${permutations}' --useSimulatedAnnealing '${useSimulatedAnnealing}' --theta '${theta}' --blastAlignmentStrings '${blastAlignmentStrings}' --baseLevel '${baseLevel}' --maxNumberOfChains '${maxNumberOfChains}' --referenceSpecies '${referenceSpecies}' --haplotypeSequences '${sequences}' --newickTree '${newickTree}' --outputDir ${outputDir} --configFile ${configFile} --referenceAlgorithms '${referenceAlgorithms}' --minimumNsForScaffoldGap ${minimumNsForScaffoldGap} --rangeOfMinimumBlockDegrees '${minimumBlockDegreeRange}' --jobTree ./jobTree ${jobTreeFlags}
	jobTreeStatus --jobTree ./jobTree --failIfNotComplete
	rm -rf ./jobTree
	
basic : pipeline

#############
#Cleaning stuff
#############

clean :
	rm -rf ${outputDir}/*
