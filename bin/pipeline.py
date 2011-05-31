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

from cactusTools.shared.common import runCactusTreeStats
from cactusTools.shared.common import runCactusMAFGenerator
from cactusTools.shared.common import runCactusAddReferenceSequence

from sonLib.bioio import getTempFile, getTempDirectory
from sonLib.bioio import fastaRead, fastaWrite
from sonLib.bioio import system

from jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

def getRootPathString():
    """
    function for finding external location
    """
    import os
    import referenceScripts.bin.pipeline
    i = os.path.abspath(referenceScripts.bin.pipeline.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def getCactusDiskString(alignmentFile):
    return "<st_kv_database_conf type=\"tokyo_cabinet\"><tokyo_cabinet database_dir=\"%s\"/></st_kv_database_conf>" % alignmentFile

class MakeAlignment(Target):
    """Target runs the alignment.
    """
    def __init__(self, options, outputDir, requiredSpecies,
                 referenceAlgorithm, minimumBlockDegree):
        Target.__init__(self, cpu=1, memory=8000000000)
        self.requiredSpecies = requiredSpecies
        self.outputDir = outputDir
        self.referenceAlgorithm = referenceAlgorithm
        self.minimumBlockDegree = minimumBlockDegree
        self.options = options
    
    def run(self):
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)
        cactusAlignmentName = "cactusAlignment"
        outputFile = os.path.join(self.outputDir, cactusAlignmentName)
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
                                                 sequences=self.options.haplotypeSequences.split(), 
                                                 newickTreeString=self.options.newickTree, 
                                                 requiredSpecies=self.requiredSpecies,
                                                 databaseName=cactusAlignmentName,
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
            #Copy across the final alignment
            localCactusDisk = os.path.join(self.getLocalTempDir(), cactusAlignmentName)
            #Compress and copy
            cwd = os.getcwd()
            os.chdir(self.getLocalTempDir())
            system("tar -cvf preReference.tar ./%s" % cactusAlignmentName)
            os.chdir(cwd)
            system("cp %s/preReference.tar %s/" % (self.getLocalTempDir(), self.outputDir))
            #Now add the reference sequence
            runCactusAddReferenceSequence(cactusDiskDatabaseString=getCactusDiskString(localCactusDisk), referenceEventString=self.options.referenceSpecies.split()[0]) #Moved this due to problems
            #Move the final db
            system("mv %s %s" % (localCactusDisk, outputFile))
            #Compute the stats
            system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir, self.outputDir))
            #We're done!
        self.addChildTarget(MakeStats(outputFile, self.outputDir, self.options))

class AddReferenceSequence(MakeAlignment):
    def run(self):
        #Now add the reference sequence
        localCactusDisk = os.path.join(self.getLocalTempDir(), cactusAlignmentName)
        runCactusAddReferenceSequence(cactusDiskDatabaseString=getCactusDiskString(localCactusDisk), referenceEventString=self.options.referenceSpecies.split()[0])
        #Compute the stats
        system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir, self.outputDir))
        #Copy across the final alignment
        system("mv %s %s" % (localCactusDisk, outputFile))
        #We're done!
        self.addChildTarget(MakeStats(outputFile, self.outputDir, self.options))
    
class MakeAlignments(Target):
    """Makes alignments using pipeline.
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options
        
    def run(self):
        statsFiles = []
        statsNames = []
        for requiredSpecies in (None, self.options.requiredSpecies):
            for referenceAlgorithm in self.options.referenceAlgorithms.split():
                for minimumBlockDegree in [ int(i) for i in self.options.rangeOfMinimumBlockDegrees.split() ]:
                    os.path.exists(self.options.outputDir)
                    def fn(i):
                        if i == None:
                            return "no-required-species"
                        return "required-species"
                    jobOutputDir = "%s-%s-%s" % (fn(requiredSpecies), referenceAlgorithm, minimumBlockDegree)
                    statsNames.append(jobOutputDir)
                    absJobOutputDir = os.path.join(self.options.outputDir, jobOutputDir)
                    statsFiles.append(os.path.join(absJobOutputDir, "treeStats.xml"))
                    self.addChildTarget(MakeAlignment(self.options, absJobOutputDir, 
                                                      requiredSpecies, referenceAlgorithm, minimumBlockDegree))

class MakeStats(Target):
    """Builds basic stats and the maf alignment.
    """
    def __init__(self, alignment, outputDir, options, cpu=4, memory=8000000000):
        Target.__init__(self, cpu=cpu, memory=memory)
        self.alignment = alignment
        self.outputDir = outputDir
        self.options = options
    
    def runScript(self, binaryName, outputFile, specialOptions):
        if not os.path.exists(outputFile):
            tempOutputFile = getTempFile(rootDir=self.getLocalTempDir())
            os.remove(tempOutputFile)
            system("%s --cactusDisk '%s' --outputFile %s --minimumNsForScaffoldGap %s %s" % 
            (os.path.join(getRootPathString(), "bin", binaryName),
             getCactusDiskString(self.alignment),
             tempOutputFile, 
             self.options.minimumNsForScaffoldGap, specialOptions))
            system("mv %s %s" % (tempOutputFile, outputFile))
        
    def run(self):
        for outputFile, program in (("treeStats.xml", runCactusTreeStats),
                                    ("alignment.maf", runCactusMAFGenerator)):
            outputFile = os.path.join(self.outputDir, outputFile)
            if not os.path.exists(outputFile):
                tempFile = os.path.join(self.getLocalTempDir(), "temp")
                program(tempFile, getCactusDiskString(self.alignment))
                system("mv %s %s" % (tempFile, outputFile))
        
        #Now build the different stats files..
        for outputFile, program in (("coverageStats.xml", "coverageStats"), 
                                    ("copyNumberStats.xml", "copyNumberStats")):
            outputFile = os.path.join(self.outputDir, outputFile)
            self.runScript(program, outputFile, "--referenceEventString %s" % self.options.referenceSpecies.split()[0])
        
        for outputFile, program, in (("pathStats_%s.xml", "pathStats"), 
                                     ("contiguityStats_%s.xml", "contiguityStats"), 
                                     ("snpStats_%s.xml", "snpStats")):
            for reference in self.options.referenceSpecies.split():
                self.runScript(program, os.path.join(self.outputDir, outputFile % reference), "--referenceEventString %s" % reference)
        
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
    parser.add_option("--referenceSpecies", dest="referenceSpecies")
    parser.add_option("--minimumNsForScaffoldGap", dest="minimumNsForScaffoldGap")
    
    Stack.addJobTreeOptions(parser)
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))
    
    Stack(MakeAlignments(options)).startJobTree(options)
    logger.info("Done with job tree")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from referenceScripts.bin.pipeline import *
    _test()
    main()

