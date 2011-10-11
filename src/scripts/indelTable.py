import sys
from tex import *
import xml.etree.ElementTree as ET

def fn(file):
    l = {}
    for line in [ line.split() for line in open(file, "r").readlines()[1:] ]:
        l[line[0]] = line[1:]
    referenceLine = l.pop("reference")
    #aggregateLine = l.pop("aggregate")
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
    #yield "aggregate", int(aggregateLine[samples]), int(aggregateLine[truePositives]), None, None
    yield "reference", int(referenceLine[samples]), int(referenceLine[truePositives]), None, None
    
fileHandle = open(sys.argv[3], "w")

writeDocumentPreliminaries(fileHandle)

writePreliminaries(8, fileHandle)

#writeRow(("samples", "sequence", "\% mapped", "\% mapped and contiguous", "\% contigious that mapped"), fileHandle)

writeLine(8, 1, (("Short Insertion Polymorphisms", 0, 7, 0, 0),), fileHandle)


writeLine(8, 2, (("Sample", 0, 0, 0, 1), 
              ("T\#", 1, 1, 0, 1), 
              ("Exact", 2, 4, 0, 0), 
              ("TP", 2, 2, 1, 1), 
              ("STP", 3, 3, 1, 1), 
              ("STN", 4, 4, 1, 1),
              ("Wobble", 5, 7, 0, 0), 
              ("TP", 5, 5, 1, 1), 
              ("STP", 6, 6, 1, 1), 
              ("STN", 7, 7, 1, 1)), fileHandle)

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
    writeLine(8, 1, ((sampleName, 0, 0, 0, 0), 
                     (samples, 1, 1, 0, 0),
                     (fn2(truePositives), 2, 2, 0, 0),
                     (fn2(sampleTruePositives), 3, 3, 0, 0),
                     (fn4(sampleTruePositives, sampleTrueNegatives), 4, 4, 0, 0),
                     (fn3(filteredTruePositives), 5, 5, 0, 0),
                     (fn3(filteredSampleTruePositives), 6, 6, 0, 0),
                     (fn4(filteredSampleTruePositives, filteredSampleTrueNegatives), 7, 7, 0, 0)), fileHandle, trailingLines=0)
       

writeEnd(fileHandle, "indelTable", "Exact: insertions detected in each sample with respect to HG19, matched precisely to insertions in dbsnp (location and length). \
Wobble: as exact, but allowing a match to an insertion within 5 bases of its location in dbSNP. \
T\#: Total number of insertions. \
TP: Percentage true positives, as validated by a match in dbSNP. \
STP: Percentage (sample) true positives, as validated by those reported for the sample in question. \
STN: Percentage (sample) false negatives, as validated by those reported for the sample in question. \
An NA entry denotes that the data was not available. \
Aggregate row: gives the total insertions in human samples (excluding chimp). \
Reference row: gives insertions in our reference with respect to HG19")
writeDocumentEnd(fileHandle)
fileHandle.close()