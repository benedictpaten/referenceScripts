#!/usr/bin/env python

"""Script runs cactus to build a bunch of reference genome assemblies for a given locus
"""

import os
import xml.etree.ElementTree as ET
import xml
import sys
from optparse import OptionParser

from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack

from jobTree.src.bioio import logger
from jobTree.src.bioio import setLoggingFromOptions

from cactus.shared.config import CactusWorkflowExperiment
from cactus.shared.common import runCactusWorkflow

from sonLib.bioio import getTempFile, getTempDirectory
from sonLib.bioio import fastaRead, fastaWrite
from sonLib.bioio import system

from jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

def getRootPathString():
    """
    function for finding external location
    """
    import os
    import referencePaper.bin.pipeline
    i = os.path.abspath(referencePaper.bin.pipeline.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def getCactusDiskString(alignmentFile):
    return "<st_kv_database_conf type=\"tokyo_cabinet\"><tokyo_cabinet database_dir=\"%s\"/></st_kv_database_conf>" % alignmentFile

class MakeAlignment(Target):
    """Target runs the alignment.
    """
    def __init__(self, newickTree, haplotypeSequences, requiredSpecies,
                 outputDir, referenceAlgorithm, minimumBlockDegree):
        Target.__init__(self, cpu=1, memory=8000000000)
        self.newickTree = newickTree
        self.haplotypeSequences = haplotypeSequences
        self.requiredSpecies = requiredSpecies
        self.outputDir = outputDir
        self.referenceAlgorithm = referenceAlgorithm
        self.minimumBlockDegree = minimumBlockDegree
    
    def run(self):
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)
        outputFile = os.path.join(self.outputDir, "cactusAlignment")
        if not os.path.exists(outputFile):
            config = ET.parse(os.path.join(getRootPathString(), "lib", "cactus_workflow_config.xml")).getroot()
            
            #Set the reference algorithm
            config.find("reference").attrib["matching_algorithm"] = self.referenceAlgorithm
            
            #Do the minimum block degree configuration
            minimumBlastBlockDegree = self.minimumBlockDegree
            if minimumBlastBlockDegree <= 1:
                minimumBlastBlockDegree = 2
            config.find("alignment").find("iterations").findall("iteration")[0].find("core").attrib["minimumBlockDegree"] = str(minimumBlastBlockDegree)
            config.find("alignment").find("iterations").findall("iteration")[1].attrib["minimumBlockDegree"] = str(self.minimumBlockDegree)
            
            #Write the config file
            tempConfigFile = os.path.join(self.getLocalTempDir(), "config.xml")
            fileHandle = open(tempConfigFile, 'w')
            tree = ET.ElementTree(config)
            tree.write(fileHandle)
            fileHandle.close()
         
            #Make the supporting temporary files
            tempExperimentFile = os.path.join(self.getLocalTempDir(), "experiment.xml")
            tempJobTreeDir = os.path.join(self.getLocalTempDir(), "jobTree")
            #Make the experiment file
            cactusWorkflowExperiment = CactusWorkflowExperiment(
                                                 sequences=self.haplotypeSequences, 
                                                 newickTreeString=self.newickTree, 
                                                 requiredSpecies=self.requiredSpecies,
                                                 databaseName="cactusAlignment",
                                                 outputDir=self.getLocalTempDir(),
                                                 configFile=tempConfigFile)
            cactusWorkflowExperiment.writeExperimentFile(tempExperimentFile)
            #Now run cactus workflow
            runCactusWorkflow(experimentFile=tempExperimentFile, jobTreeDir=tempJobTreeDir, 
                              setupAndBuildAlignments=True,
                              buildTrees=False, buildFaces=False, buildReference=True,
                              batchSystem="single_machine", maxThreads=1, jobTreeStats=True)
            logger.info("Ran the workflow")
            #Check if the jobtree completed sucessively.
            runJobTreeStatusAndFailIfNotComplete(tempJobTreeDir)
            logger.info("Checked the job tree dir")
            #Now copy the true assembly back to the output
            system("mv %s %s/experiment.xml" % (tempExperimentFile, self.outputDir))
            system("mv %s %s/config.xml" % (tempConfigFile, self.outputDir))
            #Compute the stats
            system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir, self.outputDir))
            #Copy across the final alignment
            system("mv %s %s" % (os.path.join(self.getLocalTempDir(), "cactusAlignment"), outputFile))
            #We're done!
        self.addChildTarget(MakeStats(outputFile, self.outputDir))
    
class MakeAlignments(Target):
    """Makes alignments using pipeline.
    """
    def __init__(self, newickTree, haplotypeSequences, requiredSpecies,
                 outputDir, configFile, referenceAlgorithms,
                 rangeOfMinimumBlockDegrees):
        Target.__init__(self)
        self.newickTree = newickTree
        self.haplotypeSequences = haplotypeSequences
        self.outputDir = outputDir
        self.configFile = configFile
        self.requiredSpecies = requiredSpecies
        self.referenceAlgorithms = referenceAlgorithms
        self.rangeOfMinimumBlockDegrees = rangeOfMinimumBlockDegrees
    
    def run(self):
        statsFiles = []
        statsNames = []
        for requiredSpecies in (None, self.requiredSpecies):
            for referenceAlgorithm in self.referenceAlgorithms:
                for minimumBlockDegree in self.rangeOfMinimumBlockDegrees:
                    os.path.exists(self.outputDir)
                    def fn(i):
                        if i == None:
                            return "no-required-species"
                        return "required-species"
                    jobOutputDir = "%s-%s-%s" % (fn(requiredSpecies), referenceAlgorithm, minimumBlockDegree)
                    statsNames.append(jobOutputDir)
                    absJobOutputDir = os.path.join(self.outputDir, jobOutputDir)
                    statsFiles.append(os.path.join(absJobOutputDir, "treeStats.xml"))
                    self.addChildTarget(MakeAlignment(self.newickTree, self.haplotypeSequences, 
                                                      self.requiredSpecies, absJobOutputDir, 
                                                      referenceAlgorithm, minimumBlockDegree))
        self.setFollowOnTarget(MakeComparativeStats(statsFiles, statsNames, self.outputDir))

class MakeStats(Target):
    """Builds basic stats and the maf alignment.
    """
    def __init__(self, alignment, outputDir, cpu=4, memory=8000000000):
        Target.__init__(self, cpu=cpu, memory=memory)
        self.alignment = alignment
        self.outputDir = outputDir
        
    def run(self):
        binPath = os.path.join(getRootPathString(), "bin")
        
        tempOutputFile = os.path.join(self.getLocalTempDir(), "tempOutput")
        treeStatsOutputFile = os.path.join(self.outputDir, "treeStats.xml")
        if not os.path.exists(treeStatsOutputFile):
            system("cactus_treeStats --cactusDisk '%s' --flowerName 0 --outputFile %s --noPerColumnStats" % (getCactusDiskString(self.alignment), tempOutputFile))
            system("mv %s %s" % (tempOutputFile, treeStatsOutputFile))
        outputFile = os.path.join(self.outputDir, "output.maf")
        if not os.path.exists(outputFile):
            system("cactus_MAFGenerator --cactusDisk '%s' --flowerName 0 --outputFile %s --orderByReference" % (getCactusDiskString(self.alignment), tempOutputFile))
            system("mv %s %s" % (tempOutputFile, outputFile))
        outputFile = os.path.join(self.outputDir, "adjacenciesVsPseudoAdjacenciesPerTerminalGroup.txt")
        if not os.path.exists(outputFile):
            system("python %s/scatterPlot.py %s reference.pseudo_adjacencies_per_terminal_group reference.true_pseudo_adjacencies_per_terminal_group %s" % \
                   (binPath, tempOutputFile, treeStatsOutputFile))
            system("mv %s %s" % (tempOutputFile, outputFile))
            
class MakeComparativeStats(Target):
    def __init__(self, statsFiles, statsNames, outputDir):
        Target.__init__(self)
        self.statsFiles = statsFiles
        self.statsNames = statsNames
        self.outputDir = outputDir
        
    def run(self):
        binPath = os.path.join(getRootPathString(), "bin")
        
        statsString = " ".join([ "%s %s" % (i, j) for (i, j) in zip(self.statsFiles, self.statsNames) ])
        
        logger.info("The stats string for making comparative stats is: %s" % statsString)
        
        system("cactus_treeStatsToLatexTables.py --outputFile %s/cactusStatsLatexTables.tex %s" % (self.outputDir, statsString))
                        
        statsFileString = " ".join(self.statsFiles)                                                                   
        system("python %s/tabulateFrequencies.py %s/contigLengths.txt reference2.contig_lengths_filtered %s" % (binPath, self.outputDir, statsString))
    
        system("python %s/tabulateFrequencies.py %s/blockLengths.txt blocks.lengths %s" % (binPath, self.outputDir, statsString))
    
        system("python %s/tabulateFrequencies.py %s/chainLengths.txt chains.base_block_lengths %s" % (binPath, self.outputDir, statsString))
    
        system("python %s/tabulateFrequencies.py %s/copyNumber.txt blocks.leaf_degrees %s" % (binPath, self.outputDir, statsString))
    
        system("python %s/tabulateFrequencies.py %s/flowerDepths.txt flowers.depths %s" % (binPath, self.outputDir, statsString))
    
        system("python %s/tabulateFrequencies.py %s/endsPerTerminalGroup.txt nets.total_end_numbers_per_terminal_group %s" % (binPath, self.outputDir, statsString))
    
        system("python %s/tabulateFrequencies.py %s/groupsPerNet.txt nets.total_groups_per_net %s" % (binPath, self.outputDir, statsString))
    
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
 
    parser.add_option("--haplotypeSequences", dest="haplotypeSequences")
    parser.add_option("--newickTree", dest="newickTree")
    parser.add_option("--outputDir", dest="outputDir")
    parser.add_option("--configFile", dest="configFile")
    parser.add_option("--referenceAlgorithms", dest="referenceAlgorithms")
    parser.add_option("--requiredSpecies", dest="requiredSpecies")
    parser.add_option("--rangeOfMinimumBlockDegrees", dest="rangeOfMinimumBlockDegrees")
    
    Stack.addJobTreeOptions(parser)
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))
    
    Stack(MakeAlignments(newickTree=options.newickTree, 
                         haplotypeSequences=options.haplotypeSequences.split(), 
                         requiredSpecies=options.requiredSpecies,
                         outputDir=options.outputDir,
                         configFile=options.configFile, 
                         referenceAlgorithms=options.referenceAlgorithms.split(),
                         rangeOfMinimumBlockDegrees=[ int(i) for i in options.rangeOfMinimumBlockDegrees.split() ])).startJobTree(options)
    logger.info("Done with job tree")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from referencePaper.bin.pipeline import *
    _test()
    main()

