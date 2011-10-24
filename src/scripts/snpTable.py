import sys
from tex import *
import xml.etree.ElementTree as ET

def fn(file):
    l = {}
    for line in [ line.split() for line in open(file, "r").readlines()[2:] ]:
        l[line[0]] = line[1:]
    referenceLine = l.pop("reference")
    aggregateLine = l.pop("aggregate")
    chimpLine = l.pop("panTro3")
    k = l.keys()
    k.sort()
    samples = 0
    truePositives = 3
    sampleTruePositives = 10
    sampleTrueNegatives = 12
    for key in k:
        if len(l[key]) > 8:
            print "gooo", key, l[key][samples], l[key][truePositives], l[key][sampleTruePositives], l[key][sampleTrueNegatives]
            yield (key, int(l[key][samples]), int(l[key][truePositives]), int(l[key][sampleTruePositives]), int(l[key][sampleTrueNegatives]))
        else:
            yield (key, int(l[key][samples]), int(l[key][truePositives]), None, None)
    yield "panTro3", int(chimpLine[samples]), int(chimpLine[truePositives]), None, None
    yield "aggregate", int(aggregateLine[samples]), int(aggregateLine[truePositives]), None, None
    yield "reference", int(referenceLine[samples]), int(referenceLine[truePositives]), None, None
    
fileHandle = open(sys.argv[3], "w")

writeDocumentPreliminaries(fileHandle)
writePreliminaries(9, fileHandle)

#writeRow(("samples", "sequence", "\% mapped", "\% mapped and contiguous", "\% contigious that mapped"), fileHandle)

writeLine(9, 1, (("Single Nucleotide Polymorphisms", 0, 8, 0, 0),), fileHandle)


writeLine(9, 2, (("Sample", 0, 0, 0, 1), 
              ("Unfiltered", 1, 4, 0, 0), 
              ("T\#", 1, 1, 1, 1), 
              ("TP", 2, 2, 1, 1), 
              ("STP", 3, 3, 1, 1), 
              ("SFN", 4, 4, 1, 1),
              ("Filtered", 5, 8, 0, 0), 
              ("T\#", 5, 5, 1, 1), 
              ("TP", 6, 6, 1, 1), 
              ("STP", 7, 7, 1, 1), 
              ("SFN", 8, 8, 1, 1)), fileHandle)

for sampleName, samples, truePositives, \
    sampleTruePositives, sampleTrueNegatives, \
    filteredSamples, filteredTruePositives,\
    filteredSampleTruePositives, filteredSampleTrueNegatives in [ tuple(list(i) + list(j[1:])) for i, j in zip(fn(sys.argv[1]), fn(sys.argv[2])) if i[0] == j[0] ]:
    def fn2(i, j=samples):
        if i != None:
            return "%.0f" % (100.0*float(i)/float(j))
        return "NA"
    def fn3(i, j=filteredSamples):
        return fn2(i, j)
    def fn4(i, j):
        if i == None:
            return "NA"
        return "%.0f" % (100.0*float(j)/(float(i) + float(j)))
    writeLine(9, 1, ((sampleName, 0, 0, 0, 0), 
                     (samples, 1, 1, 0, 0),
                     (fn2(truePositives), 2, 2, 0, 0),
                     (fn2(sampleTruePositives), 3, 3, 0, 0),
                     (fn4(sampleTruePositives, sampleTrueNegatives), 4, 4, 0, 0),
                     (filteredSamples, 5, 5, 0, 0),
                     (fn3(filteredTruePositives), 6, 6, 0, 0),
                     (fn3(filteredSampleTruePositives), 7, 7, 0, 0),
                     (fn4(filteredSampleTruePositives, filteredSampleTrueNegatives), 8, 8, 0, 0)), fileHandle, trailingLines=0)
       

writeEnd(fileHandle, "snpTable", "Unfiltered: all SNPs detected in each sample with respect to HG19. \
Filtered: as unfiltered, but excluding SNPs detected within 5 bps of an indel within the MSA. \
T\#: Total number of SNPs. \
TP: Percentage true positives, as validated by a matching SNP in dbSNP. \
STP: Percentage (sample) true positives, as validated by those reported for the sample in question. \
SFN: Percentage (sample) false negatives, as validated by those reported for the sample in question. \
An NA entry denotes that the data was not available. \
Aggregate row: gives the total SNPs in human samples (excluding chimp). \
Reference row: gives SNPs in our reference with respect to HG19")
writeDocumentEnd(fileHandle)
fileHandle.close()