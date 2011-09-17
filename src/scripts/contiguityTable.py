import sys
from tex import *
import xml.etree.ElementTree as ET

ref = ET.parse(sys.argv[1]).getroot()
hg19 = ET.parse(sys.argv[2]).getroot()
fileHandle = open(sys.argv[3], 'w')

writeDocumentPreliminaries(fileHandle)
writePreliminaries(5, fileHandle)

def fn(x):
    return "%.2f" % (100*x)

#writeRow(("samples", "sequence", "\% mapped", "\% mapped and contiguous", "\% contigious that mapped"), fileHandle)

writeLine(5, 1, (("Contiguity Statistics", 0, 4, 0, 0),), fileHandle)

writeLine(5, 2, (("samples", 0, 0, 0, 1), 
              ("sequence", 1, 1, 0, 1), 
              ("\% mapped", 2, 2, 0, 1),
              ("\% mapped and contiguous", 3, 3, 0, 1),
              ("\% contigious that mapped", 4, 4, 0, 1)), fileHandle)

samples = [ i.attrib["sampleName"] for i in ref.findall("statsForSample") if i.attrib["sampleName"] not in ("hg19", "reference", "", "ROOT")]
samples.sort()

refTotalSamples = 0.0
refTotalAligned = 0.0
refTotalContiguous = 0.0
hg19TotalSamples = 0.0
hg19TotalAligned = 0.0
hg19TotalContiguous = 0.0

for sample in samples:
    refSample = [ i for i in ref.findall("statsForSample") if i.attrib["sampleName"] == sample ][0]
    hg19Sample = [ i for i in hg19.findall("statsForSample") if i.attrib["sampleName"] == sample ][0]
    refTotalSamples += float(refSample.attrib["totalSamples"])
    refTotalAligned += float(refSample.attrib["totalAligned"])
    refTotalContiguous += float(refSample.attrib["totalCorrect"])
    hg19TotalSamples += float(hg19Sample.attrib["totalSamples"])
    hg19TotalAligned += float(hg19Sample.attrib["totalAligned"])
    hg19TotalContiguous += float(hg19Sample.attrib["totalCorrect"])
    writeLine(5, 2, ((sample, 0, 0, 0, 1), 
                     ("reference", 1, 1, 0, 0), 
                     ("hg19", 1, 1, 1, 1), 
                     (fn(float(refSample.attrib["totalAligned"])/float(refSample.attrib["totalSamples"])), 2, 2, 0, 0), 
                     (fn(float(hg19Sample.attrib["totalAligned"])/float(hg19Sample.attrib["totalSamples"])), 2, 2, 1, 1),
                     (fn(float(refSample.attrib["correctPerSample"])), 3, 3, 0, 0), 
                     (fn(float(hg19Sample.attrib["correctPerSample"])), 3, 3, 1, 1),
                     (fn(float(refSample.attrib["correctPerAligned"])), 4, 4, 0, 0), 
                     (fn(float(hg19Sample.attrib["correctPerAligned"])), 4, 4, 1, 1),
                     ), fileHandle, trailingLines=1)

writeLine(5, 2, (("aggregate", 0, 0, 0, 1), 
                     ("reference", 1, 1, 0, 0), 
                     ("hg19", 1, 1, 1, 1), 
                     (fn(refTotalAligned/refTotalSamples), 2, 2, 0, 0), 
                     (fn(hg19TotalAligned/hg19TotalSamples), 2, 2, 1, 1),
                     (fn(refTotalContiguous/refTotalSamples), 3, 3, 0, 0), 
                     (fn(hg19TotalContiguous/hg19TotalSamples), 3, 3, 1, 1),
                     (fn(refTotalContiguous/refTotalAligned), 4, 4, 0, 0), 
                     (fn(hg19TotalContiguous/hg19TotalAligned), 4, 4, 1, 1),
                     ), fileHandle, trailingLines=1)

writeEnd(fileHandle, "contiguityTable", "Statistics on correct contiguity.")
writeDocumentEnd(fileHandle)