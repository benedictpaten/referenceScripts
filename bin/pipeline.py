#!/usr/bin/env python

"""Script runs cactus to build a bunch of reference genome assemblies for a given locus
"""

import os
import xml.etree.ElementTree as ET
import xml
import sys
from threading import Thread
from optparse import OptionParser
import time
import random

from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack

from sonLib.bioio import logger
from sonLib.bioio import setLoggingFromOptions

from cactus.shared.config import CactusWorkflowExperiment
from cactus.shared.common import runCactusWorkflow

from cactusTools.shared.common import runCactusTreeStats
from cactusTools.shared.common import runCactusMAFGenerator

from sonLib.bioio import getTempFile, getTempDirectory
from sonLib.bioio import fastaRead, fastaWrite, cigarRead, cigarWrite
from sonLib.bioio import system

from jobTree.src.common import runJobTreeStatusAndFailIfNotComplete

from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.configWrapper import ConfigWrapper
from cactus.progressive.cactus_createMultiCactusProject import cleanEventTree
from cactus.pipeline.ktserverControl import runKtserver
from cactus.pipeline.ktserverControl import killKtServer
from cactus.pipeline.ktserverControl import blockUntilKtserverIsRunnning
from cactus.pipeline.cactus_workflow import CactusPhasesTarget

from cactus.pipeline.ktserverJobTree import addKtserverDependentChild

def getRootPathString():
    """
    function for finding external location
    """
    import os
    import referenceScripts.bin.pipeline
    i = os.path.abspath(referenceScripts.bin.pipeline.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def getCactusDiskString(options, outputDir):
    dbConfElem = ET.Element("st_kv_database_conf", type="kyoto_tycoon")
    ET.SubElement(dbConfElem, "kyoto_tycoon", host=options.databaseHost, port="1978", database_dir=outputDir, in_memory="true", snapshot="true")                
    return dbConfElem

class A:
    def __init__(self, cactusWorkflowExperiment, experimentFile):
        self.experimentNode = cactusWorkflowExperiment
        self.experimentFile = experimentFile

class MakeAlignment(Target):
    """Target runs the alignment.
    """
    def __init__(self, options,
                 sequences, 
                 constraints,
                 outputDir, 
                 referenceAlgorithm, minimumBlockDegree, 
                 blastAlignmentString, baseLevel, maxNumberOfChains, permutations,
                 theta, useSimulatedAnnealing, heldOutSequence, pruneOutStubAlignments, gapGamma,
                 singleCopyIngroup):
        Target.__init__(self, cpu=int(options.maxThreads), memory=8000000000)
        self.sequences = sequences
        self.constraints = constraints
        self.outputDir = outputDir
        self.referenceAlgorithm = referenceAlgorithm
        self.minimumBlockDegree = int(minimumBlockDegree)
        self.blastAlignmentString = blastAlignmentString
        self.baseLevel = baseLevel
        self.maxNumberOfChains = maxNumberOfChains
        self.permutations = permutations
        self.theta = theta
        self.useSimulatedAnnealing = useSimulatedAnnealing
        self.options = options
        self.heldOutSequence = heldOutSequence
        self.pruneOutStubAlignments = pruneOutStubAlignments
        self.gapGamma = gapGamma
        self.singleCopyIngroup = singleCopyIngroup
    
    def run(self):
        if not os.path.isdir(self.outputDir):
            os.mkdir(self.outputDir)
        cactusAlignmentName = "cactusAlignment"
        outputFile = os.path.join(self.outputDir, cactusAlignmentName)
        if not os.path.exists(outputFile):
            #Set up the database
            dbConfElem = getCactusDiskString(self.options, self.outputDir)
            config = ET.parse(os.path.join(getRootPathString(), "lib", "cactus_workflow_config.xml")).getroot()
            
            #Set the reference algorithm
            config.find("reference").attrib["matchingAlgorithm"] = self.referenceAlgorithm
            
            #Do the minimum block degree configuration
            blastIteration = config.find("caf") 
            baseIteration = config.find("bar")
            
            minimumBlastBlockDegree = self.minimumBlockDegree
            if minimumBlastBlockDegree <= 1:
                minimumBlastBlockDegree = 2
            blastIteration.attrib["minimumBlockDegree"] = str(minimumBlastBlockDegree)
            if not self.baseLevel:
                baseIteration.attrib["runBar"] = "0"
            baseIteration.attrib["minimumBlockDegree"] = str(self.minimumBlockDegree)
            baseIteration.attrib["pruneOutStubAlignments"] = str(int(self.pruneOutStubAlignments))
            baseIteration.attrib["gapGamma"] = str(float(self.gapGamma))
            
            #Set the blast string
            blastIteration.attrib["lastzArguments"] = blastIteration.attrib["lastzArguments"].replace("PARAMETERS", self.blastAlignmentString)
            blastIteration.attrib["singleCopyIngroup"] = str(int(self.singleCopyIngroup))
            
            #Set the number of chains to allow in a level, during promotion
            config.find("normal").attrib["maxNumberOfChains"] = str(self.maxNumberOfChains)
            
            #Set the number of chains to order per round of the matching algorithm
            config.find("reference").attrib["permutations"]  = str(self.permutations)
            
            #Set the chain weight function
            if bool(self.useSimulatedAnnealing):
                config.find("reference").attrib["useSimulatedAnnealing"]="1"
                
            config.find("reference").attrib["theta"] = str(self.theta)
            
            #Write the config file
            configFile = os.path.join(self.outputDir, "config.xml")
            fileHandle = open(configFile, 'w')
            tree = ET.ElementTree(config)
            tree.write(fileHandle)
            fileHandle.close()
            
            #Make the supporting temporary files
            experimentFile = os.path.join(self.outputDir, "experiment.xml")
            jobTreeDir = os.path.join(self.outputDir, "jobTree")
            c2hFile = os.path.join(self.outputDir, "out.c2h")
            fastaFile = os.path.join(self.outputDir, "out.fa")
            
            #Make the logfile for jobTree
            jobTreeLogFile = os.path.join(self.outputDir, "log.txt")
            
            #Make the experiment file
            cactusWorkflowExperiment = CactusWorkflowExperiment(
                                                 sequences=self.sequences.split(), 
                                                 newickTreeString=self.options.newickTree, 
                                                 outgroupEvents = self.options.outgroupEvent,
                                                 databaseName=cactusAlignmentName,
                                                 halFile=c2hFile,
                                                 fastaFile=fastaFile,
                                                 outputDir=self.outputDir,
                                                 configFile=configFile,
                                                 constraints=self.constraints,
                                                 databaseConf = dbConfElem)
            cactusWorkflowExperiment.writeExperimentFile(experimentFile)
        
            #Now run cactus workflow
            runCactusWorkflow(experimentFile=experimentFile, jobTreeDir=jobTreeDir, 
                              buildAvgs=False, buildReference=True,
                              maxThreads=int(self.options.maxThreads), jobTreeStats=True,
                              defaultMemory=int(self.options.defaultMemory),
                              batchSystem=self.options.batchSystemForAlignments,
                              retryCount=2,
                              logLevel="DEBUG",
                              logFile = jobTreeLogFile,
                              extraJobTreeArgumentsString="--parasolCommand '%s'" % self.options.parasolCommandForAlignment)
            logger.info("Ran the workflow")
            #Check if the jobtree completed sucessively.
            runJobTreeStatusAndFailIfNotComplete(jobTreeDir)
            logger.info("Checked the job tree dir")
            #Now copy the true assembly back to the output
            system("cp %s %s/constraints.cig" % (self.constraints, self.outputDir))
            halFile = os.path.join(self.outputDir, "out.hal")
            system("halAppendCactusSubtree %s %s '%s' %s" % (c2hFile, fastaFile, self.options.newickTree, halFile))
            
            #Compute the stats
            #system("jobTreeStats --jobTree %s --outputFile %s/jobTreeStats.xml" % (tempJobTreeDir,  v))

        #self.addChildTarget(MakeStats(outputFile, self.outputDir, self.options))
        self.options.constraints = self.constraints #re-assign the constraint file
        self.addChildTarget(MakeAlignment2(self.outputDir, self.options))
    
class MakeAlignment2(Target):
    def __init__(self, outputDir, options):
        Target.__init__(self)
        self.outputDir = outputDir
        self.options = options
    
    def run(self):
        time.sleep(random.random()*100)
        experimentFile = os.path.join(self.outputDir, "experiment.xml")
        addKtserverDependentChild(self, MakeStats(A(ET.parse(experimentFile).getroot(), experimentFile), self.outputDir, self.options), isSecondary = False)

def makeHeldOutAlignments(self, options, outputDir, 
                 referenceAlgorithm, minimumBlockDegree, 
                 blastAlignmentString, baseLevel, maxNumberOfChains, permutations,
                 theta, useSimulatedAnnealing, pruneOutStubAlignments, gapGamma, singleCopyIngroup):
    nullSequence = os.path.join(self.getGlobalTempDir(), "nullSequence.fa")
    open(nullSequence, 'w').close()
    for heldoutSequence in self.options.heldOutSequences.split():
        heldOutSequenceNames = set()
        heldOutOutputDir = outputDir + "_" + heldoutSequence
        def fn(i):
            if heldoutSequence == i.split("/")[-1]:
                #Add sequence to list of held out sequence names
                fileHandle = open(i, 'r')
                for header, sequence in fastaRead(fileHandle):
                    name = header.split()[0]
                    assert name not in heldOutSequenceNames
                    heldOutSequenceNames.add(name)
                fileHandle.close()
                return nullSequence
            return i
        heldOutSequences = " ".join([ fn(i) for i in options.haplotypeSequences.split() ])
        #Filter the constraints
        heldOutConstraints = os.path.join(self.getGlobalTempDir(), "constraints_%s.cig" % heldoutSequence)
        heldOutConstraintsFileHandle = open(heldOutConstraints, 'w')
        constraintFileHandle = open(options.constraints, 'r')
        for cigar in cigarRead(constraintFileHandle):
            if cigar.contig1 not in heldOutSequenceNames and cigar.contig2 not in heldOutSequenceNames:
                cigarWrite(heldOutConstraintsFileHandle, cigar, withProbs=False)
        heldOutConstraintsFileHandle.close()
        constraintFileHandle.close()
        
        self.addChildTarget(MakeAlignment(options, 
                      heldOutSequences,
                      heldOutConstraints,
                      heldOutOutputDir, 
                      referenceAlgorithm, minimumBlockDegree, 
                      blastAlignmentString, baseLevel, maxNumberOfChains, permutations,
                      theta, useSimulatedAnnealing, heldoutSequence, pruneOutStubAlignments, gapGamma, singleCopyIngroup))
    self.addChildTarget(MakeAlignment(options, 
                  options.haplotypeSequences,
                  options.constraints,
                  outputDir, 
                  referenceAlgorithm, minimumBlockDegree, 
                  blastAlignmentString, baseLevel, maxNumberOfChains, permutations,
                  theta, useSimulatedAnnealing, None, pruneOutStubAlignments, gapGamma, singleCopyIngroup))

class MakeAlignments(Target):
    """Makes alignments using pipeline.
    """
    def __init__(self, options):
        Target.__init__(self)
        self.options = options
    
    def run(self):
        statsFiles = []
        statsNames = []
        for gapGamma in self.options.gapGamma.split():
            for singleCopyIngroup in [ int(i) for i in self.options.singleCopyIngroup.split() ]:
                for pruneOutStubAlignments in (True,): # False):
                    for referenceAlgorithm in self.options.referenceAlgorithms.split():
                        for minimumBlockDegree in [ int(i) for i in self.options.rangeOfMinimumBlockDegrees.split() ]:
                            blastAlignmentStrings = self.options.blastAlignmentStrings.split("%")
                            for blastAlignmentStringIndex in xrange(len(blastAlignmentStrings)):
                                for baseLevel in [ bool(int(i)) for i in self.options.baseLevel.split() ]:
                                    for maxNumberOfChains in [ int(i) for i in self.options.maxNumberOfChains.split() ]:
                                        for permutations in [ int(i) for i in self.options.permutations.split() ]:
                                            for theta in [ float(i) for i in self.options.theta.split() ]:
                                                for useSimulatedAnnealing in [ bool(int(i)) for i in self.options.useSimulatedAnnealing.split() ]:
                                                    os.path.exists(self.options.outputDir)
                                                    jobOutputDir = "%s-%s-%s-%s-%s-%s-%s-%s-%s-%s-%s" % (referenceAlgorithm, minimumBlockDegree, blastAlignmentStringIndex, \
                                                                                                            baseLevel, maxNumberOfChains, permutations, theta, useSimulatedAnnealing, pruneOutStubAlignments, gapGamma, singleCopyIngroup)
                                                    statsNames.append(jobOutputDir)
                                                    absJobOutputDir = os.path.join(self.options.outputDir, jobOutputDir)
                                                    statsFiles.append(os.path.join(absJobOutputDir, "treeStats.xml"))
                                                    makeHeldOutAlignments(self, self.options, absJobOutputDir, 
                                                                          referenceAlgorithm, minimumBlockDegree, 
                                                                          blastAlignmentStrings[blastAlignmentStringIndex], 
                                                                          baseLevel, maxNumberOfChains, permutations, 
                                                                          theta, useSimulatedAnnealing, pruneOutStubAlignments, gapGamma, singleCopyIngroup)

class MakeStats(CactusPhasesTarget):
    """Builds basic stats and the maf alignment.
    """
    def __init__(self, cactusWorkflowArguments, outputDir, options, cpu=1, memory=4000000000):
        Target.__init__(self, cpu=cpu, memory=memory)
        self.cactusWorkflowArguments = cactusWorkflowArguments
        self.outputDir = outputDir
        self.options = options
    
    def runScript(self, binaryName, outputFile, specialOptions):
        self.addChildTarget(RunScript(self.cactusWorkflowArguments, self.outputDir, self.options, binaryName, outputFile, specialOptions))
        
    def run(self):
        self.addChildTarget(MakeStats1(self.cactusWorkflowArguments, self.outputDir, self.options))    
        self.addChildTarget(MakeStats2(self.cactusWorkflowArguments, self.outputDir, self.options))    
        self.addChildTarget(MakeStats3(self.cactusWorkflowArguments, self.outputDir, self.options))
        self.addChildTarget(MakeStats4(self.cactusWorkflowArguments, self.outputDir, self.options))
        self.addChildTarget(MakeStats5(self.cactusWorkflowArguments, self.outputDir, self.options))
        self.addChildTarget(MakeStats6(self.cactusWorkflowArguments, self.outputDir, self.options))    
        self.addChildTarget(MakeStats7(self.cactusWorkflowArguments, self.outputDir, self.options))    
        self.addChildTarget(MakeAssemblyHub(self.cactusWorkflowArguments, self.outputDir, self.options)) 
        self.addChildTarget(MakeCheckConstraints(self.cactusWorkflowArguments, self.outputDir, self.options))

class RunScript(MakeStats):
    """Builds basic stats and the maf alignment.
    """
    def __init__(self, cactusWorkflowArguments, outputDir, options, binaryName, outputFile, specialOptions):
        MakeStats.__init__(self, cactusWorkflowArguments, outputDir, options)
        self.binaryName = binaryName
        self.outputFile = outputFile
        self.specialOptions = specialOptions
  
    def run(self):
        if not os.path.exists(self.outputFile):
            tempOutputFile = getTempFile(rootDir=self.getGlobalTempDir())
            os.remove(tempOutputFile)
            system("%s --cactusDisk '%s' --outputFile %s --minimumNsForScaffoldGap %s --sampleNumber %s %s" % 
            (os.path.join(getRootPathString(), "bin", self.binaryName),
             self.cactusWorkflowArguments.cactusDiskDatabaseString, #self.alignment),
             tempOutputFile, 
             self.options.minimumNsForScaffoldGap, self.options.sampleNumber, self.specialOptions))
            system("mv %s %s" % (tempOutputFile, self.outputFile))
        
class MakeAssemblyHub(MakeStats):
    def run(self):          
        #Make the assembly hub if faToTwoBit is installed
        makeAssemblyHub = True
        try:
            system("which faToTwoBit")
        except RuntimeError:
            makeAssemblyHub = False
        if makeAssemblyHub:
            cwd = os.getcwd()
            os.chdir(self.outputDir)
            system("hal2assemblyHub.py out.hal outBrowser --lod --shortLabel='%s' --longLabel='%s'" % \
                   (self.outputDir[-10:], self.outputDir))
            os.chdir(cwd)
        else:
            self.logToMaster("Not building assembly hub")

class MakeCheckConstraints(MakeStats):
    def run(self):
        outputFile = os.path.join(self.outputDir, "constraintStats.txt")
        program = "checkConstraints"
        if not os.path.exists(outputFile):
            tempOutputFile = getTempFile(rootDir=self.getLocalTempDir())
            os.remove(tempOutputFile)
            system( "%s --cactusDisk '%s' --outputFile %s --constraintFile %s" % 
             ( os.path.join(getRootPathString(), "bin", program),
               self.cactusWorkflowArguments.cactusDiskDatabaseString, #self.alignment),
               tempOutputFile,
               self.options.constraints ) )
            system("mv %s %s" % (tempOutputFile, outputFile))

class MakeStats1(MakeStats):
    def run(self):  
        for outputFile, program in (("treeStats.xml", runCactusTreeStats),
                                    ("alignment.maf", runCactusMAFGenerator),
                                    ("alignment_substitutionsOnly.maf", runCactusMAFGenerator)):
            outputFile = os.path.join(self.outputDir, outputFile)
            if not os.path.exists(outputFile):
                tempFile = os.path.join(self.getLocalTempDir(), "temp")
                if "alignment_substitutionsOnly.maf" in outputFile:
                    program(tempFile, self.cactusWorkflowArguments.cactusDiskDatabaseString, showOnlySubstitutionsWithRespectToTheReference=True)
                else:
                    program(tempFile, self.cactusWorkflowArguments.cactusDiskDatabaseString)
                system("mv %s %s" % (tempFile, outputFile))

class MakeStats2(MakeStats):
    def run(self):  
        #Now build the different stats files..
        for outputFile, program in (("coverageStats.xml", "coverageStats"), 
                                    ("copyNumberStats.xml", "copyNumberStats"),
                                    ("filterNonComponentSequences.xml", "filterNonComponentSequences")):
            outputFile = os.path.join(self.outputDir, outputFile)
            ref1, ref2 = self.options.referenceSpecies.split()
            self.runScript(program, outputFile, "--referenceEventString %s --otherReferenceEventString %s --outgroupEventString %s" % (ref1, ref2, self.options.outgroupEvent))
        #self.addChildTarget(MakeStats3(self.alignment, self.outputDir, self.options))

class MakeStats3(MakeStats):
    def run(self):     
        outputFile = "contiguityStats_%s.xml"
        program =  "contiguityStats"
        specialOptions = ""  
        for reference in self.options.referenceSpecies.split():
            self.runScript(program, os.path.join(self.outputDir, outputFile % reference), "--referenceEventString %s %s" % (reference, specialOptions))
        ref1, ref2 = self.options.referenceSpecies.split()
        outputFile = "contiguityStats_intersection_%s.xml"
        self.runScript(program, os.path.join(self.outputDir, outputFile % ref1), "--referenceEventString %s --otherReferenceEventString %s %s" % (ref1, ref2, specialOptions))
        self.runScript(program, os.path.join(self.outputDir, outputFile % ref2), "--referenceEventString %s --otherReferenceEventString %s %s" % (ref2, ref1, specialOptions))
        #self.addChildTarget(MakeStats4(self.alignment, self.outputDir, self.options))

class MakeStats4(MakeStats):
    def run(self):          
        for outputFile, program, specialOptions in (("pathStats_%s.xml", "pathStats", ""), 
                                                    ("pathStats_ignoreAdjacencies_%s.xml", "pathStats", "--ignoreAdjacencies")
                                     ):
            ref1, ref2 = self.options.referenceSpecies.split()
            self.runScript(program, os.path.join(self.outputDir, outputFile % ref1), "--referenceEventString %s %s --otherReferenceEventString %s" % (ref1, specialOptions, ref2))
            self.runScript(program, os.path.join(self.outputDir, outputFile % ref2), "--referenceEventString %s %s --otherReferenceEventString %s" % (ref2, specialOptions, ref1))
            #for reference in self.options.referenceSpecies.split():
            #    self.runScript(program, os.path.join(self.outputDir, outputFile % reference), "--referenceEventString %s %s" % (reference, specialOptions))
        self.setFollowOnTarget(MakeStats4B(self.cactusWorkflowArguments, self.outputDir, self.options))

class MakeStats4B(MakeStats):
    def run(self):
        for reference in self.options.referenceSpecies.split():
            system("python %s %s" % (os.path.join(getRootPathString(), "src", "scripts", "indelIntersection.py"), os.path.join(self.outputDir, "pathStats_%s.xml") % reference))
            system("python %s %s" % (os.path.join(getRootPathString(), "src", "scripts", "indelIntersection.py"), os.path.join(self.outputDir, "pathStats_ignoreAdjacencies_%s.xml") % reference))

class MakeStats7(MakeStats):
    def run(self):          
        for outputFile, program, specialOptions in (("snpStats_%s.xml", "snpStats", "--reportDistanceMatrix"),
                                                    ("snpStats_filtered_%s.xml", "snpStats", "--ignoreFirstNBasesOfBlock 5"),
                                                    ("snpStats_%s_recurrent.xml", "snpStats", "--minimumRecurrence 2"),
                                                    ("snpStats_filtered_%s_recurrent.xml", "snpStats", "--ignoreFirstNBasesOfBlock 5 --minimumRecurrence 2 --reportDistanceMatrix")
                                     ):
            ref1, ref2 = self.options.referenceSpecies.split()
            self.runScript(program, os.path.join(self.outputDir, outputFile % ref1), "--referenceEventString %s %s" % (ref1, specialOptions))
            self.runScript(program, os.path.join(self.outputDir, outputFile % ref2), "--referenceEventString %s %s " % (ref2, specialOptions))
            #for reference in self.options.referenceSpecies.split():
            #    self.runScript(program, os.path.join(self.outputDir, outputFile % reference), "--referenceEventString %s %s" % (reference, specialOptions))
        self.setFollowOnTarget(MakeStats7B(self.cactusWorkflowArguments, self.outputDir, self.options))

class MakeStats7B(MakeStats):
    def run(self):  
        for reference in self.options.referenceSpecies.split():
            system("python %s %s" % (os.path.join(getRootPathString(), "src", "scripts", "snpIntersection.py"), os.path.join(self.outputDir, "snpStats_%s.xml") % reference))
            system("python %s %s" % (os.path.join(getRootPathString(), "src", "scripts", "snpIntersection.py"), os.path.join(self.outputDir, "snpStats_filtered_%s.xml") % reference))
            system("python %s %s" % (os.path.join(getRootPathString(), "src", "scripts", "snpIntersection.py"), os.path.join(self.outputDir, "snpStats_%s_recurrent.xml") % reference))
            system("python %s %s" % (os.path.join(getRootPathString(), "src", "scripts", "snpIntersection.py"), os.path.join(self.outputDir, "snpStats_filtered_%s_recurrent.xml") % reference))
        #self.addChildTarget(MakeStats5(self.alignment, self.outputDir, self.options))

class MakeStats5(MakeStats):
    def run(self):       
        ref1, ref2 = self.options.referenceSpecies.split()
        self.runScript("snpStats", os.path.join(self.outputDir, "snpStatsIntersection_%s.xml" % ref1), "--referenceEventString %s --otherReferenceEventString %s" % (ref1, ref2))
        self.runScript("snpStats", os.path.join(self.outputDir, "snpStatsIntersection_%s.xml" % ref2), "--referenceEventString %s --otherReferenceEventString %s" % (ref2, ref1))   
        self.runScript("snpStats", os.path.join(self.outputDir, "snpStats_ignoreSitesWithOtherReferencePresent_%s.xml" % ref1), "--referenceEventString %s --otherReferenceEventString %s --ignoreSitesWithOtherReferencePresent" % (ref1, ref2))
        self.runScript("snpStats", os.path.join(self.outputDir, "snpStats_ignoreSitesWithOtherReferencePresent_%s.xml" % ref2), "--referenceEventString %s --otherReferenceEventString %s --ignoreSitesWithOtherReferencePresent" % (ref2, ref1))                 

class MakeStats6(MakeStats):
    def run(self):
        self.runScript("danielAlignment", os.path.join(self.outputDir, "danielAlignment.txt"), "--referenceEventString hg19 --otherReferenceEventString NA12891")
        self.runScript("sequenceCoverages", os.path.join(self.outputDir, "sequenceCoverages.txt"), "--referenceEventString reference")

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
    parser.add_option("--rangeOfMinimumBlockDegrees", dest="rangeOfMinimumBlockDegrees")
    parser.add_option("--referenceSpecies", dest="referenceSpecies")
    parser.add_option("--minimumNsForScaffoldGap", dest="minimumNsForScaffoldGap")
    parser.add_option("--blastAlignmentStrings", dest="blastAlignmentStrings")
    parser.add_option("--baseLevel", dest="baseLevel")
    parser.add_option("--maxNumberOfChains", dest="maxNumberOfChains")
    parser.add_option("--permutations", dest="permutations")
    parser.add_option("--theta", dest="theta")
    parser.add_option("--useSimulatedAnnealing", dest="useSimulatedAnnealing")
    parser.add_option("--sampleNumber", dest="sampleNumber")
    parser.add_option("--heldOutSequences", dest="heldOutSequences")
    parser.add_option("--outgroupEvent", dest="outgroupEvent")
    parser.add_option("--gapGamma", dest="gapGamma")
    parser.add_option("--constraints", dest="constraints")
    parser.add_option("--singleCopyIngroup", dest="singleCopyIngroup")
    
    parser.add_option("--databaseHost", default="localhost")
    parser.add_option("--batchSystemForAlignments", default="singleMachine")
    parser.add_option("--parasolCommandForAlignment", default="parasol")
    
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
