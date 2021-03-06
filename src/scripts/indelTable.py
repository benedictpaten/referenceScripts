import sys
from tex import *
import xml.etree.ElementTree as ET

def fn(file, insertionOrDeletion):
    l = {}
    for line in [ line.split() for line in open(file, "r").readlines()[2:] if line.split()[0] == insertionOrDeletion ]:
        l[line[1]] = line[1:]
    referenceLine = l.pop("reference")
    aggregateLine = l.pop("aggregate")
    chimpLine = l.pop("panTro3")
    k = l.keys()
    k.sort()
    samples = 1
    truePositives = 4
    sampleTruePositives = 11
    sampleTrueNegatives = 13
    for key in k:
        if len(l[key]) > 8:
            print "gooo", key, l[key][samples], l[key][truePositives], l[key][sampleTruePositives], l[key][sampleTrueNegatives]
            yield (key, int(l[key][samples]), int(l[key][truePositives]), int(l[key][sampleTruePositives]), int(l[key][sampleTrueNegatives]))
        else:
            yield (key, int(l[key][samples]), int(l[key][truePositives]), None, None)
    yield "panTro3", int(chimpLine[samples]), int(chimpLine[truePositives]), None, None
    yield "aggregate", int(aggregateLine[samples]), int(aggregateLine[truePositives]), None, None
    yield "C. Ref.", int(referenceLine[samples]), int(referenceLine[truePositives]), None, None
    
fileHandle = open(sys.argv[3], "w")

writeDocumentPreliminaries(fileHandle)

#writeRow(("samples", "sequence", "\% mapped", "\% mapped and contiguous", "\% contigious that mapped"), fileHandle)

for type in ("insertion", "deletion"): 
    writePreliminaries(8, fileHandle)
    writeLine(8, 1, (("Short %s Polymorphisms" % type, 0, 7, 0, 0),), fileHandle)
    
    writeLine(8, 2, (("Sample", 0, 0, 0, 1), 
                  ("T\#", 1, 1, 0, 1), 
                  ("All", 2, 4, 0, 0), 
                  ("TP", 2, 2, 1, 1), 
                  ("STP", 3, 3, 1, 1), 
                  ("SFN", 4, 4, 1, 1),
                  ("No wobble", 5, 7, 0, 0), 
                  ("TP", 5, 5, 1, 1), 
                  ("STP", 6, 6, 1, 1), 
                  ("SFN", 7, 7, 1, 1)), fileHandle)
    
    for sampleName, samples, truePositives, \
        sampleTruePositives, sampleTrueNegatives, \
        filteredSamples, filteredTruePositives,\
        filteredSampleTruePositives, filteredSampleTrueNegatives in [ tuple(list(i) + list(j[1:])) for i, j in zip(fn(sys.argv[1], type), fn(sys.argv[2], type)) if i[0] == j[0] ]:
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

    writeEnd(fileHandle, "%sTable" % type, "All: %ss detected in each sample with respect to GRCh37, allowing a match to an %s within 5 bases of its location in dbSNP. \
    No wobble: as All, matched precisely to insertions in dbsnp (location and length) \
    T\#: Total number of %ss. \
    TP: Percentage true positives, as validated by a match in dbSNP. \
    STP: Percentage (sample) true positives, as validated by those reported for the sample in question. \
    SFN: Percentage (sample) false negatives, as validated by those reported for the sample in question. \
    An NA entry denotes that the data was not available. \
    Aggregate row: gives the total %ss in human samples (excluding chimp). \
    C. Ref. row: gives %ss in C. Ref. with respect to GRCh37" % (type, type, type, type, type))
    
writeDocumentEnd(fileHandle)
fileHandle.close()