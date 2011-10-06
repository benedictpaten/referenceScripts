import sys
from tex import *
import xml.etree.ElementTree as ET

def fn(file):
    l = {}
    for line in [ line.split() for line in open(file, "r").readlines()[1:] ]:
        l[line[0]] = line[1:]
    l.pop("reference")
    aggregateLine = l.pop("aggregate")
    chimpLine = l.pop("panTro3")
    k = l.keys()
    k.sort()
    for key in k:
        yield (key, int(l[key][0]), int(l[key][1]))
    yield "panTro3", int(chimpLine[0]), int(chimpLine[1])
    yield "aggregate", int(aggregateLine[0]), int(aggregateLine[1])
    
fileHandle = open(sys.argv[3], "w")

writeDocumentPreliminaries(fileHandle)
writePreliminaries(5, fileHandle)

#writeRow(("samples", "sequence", "\% mapped", "\% mapped and contiguous", "\% contigious that mapped"), fileHandle)

writeLine(3, 1, (("Single Nucleotide Polymorphisms", 0, 2, 0, 0),), fileHandle)

writeLine(3, 1, (("Sample", 0, 0, 0, 0), 
              ("Total", 1, 1, 0, 0), 
              ("Total Filtered", 2, 2, 0, 0)), fileHandle)

for sampleName, samples, truePositives, filteredSamples, filteredTruePositives in [ tuple(list(i) + list(j[1:])) for i, j in zip(fn(sys.argv[1]), fn(sys.argv[2])) if i[0] == j[0] ]:
    writeLine(3, 1, ((sampleName, 0, 0, 0, 0), 
                     ("%s (%.1f)" % (samples, 100.0*truePositives/samples), 1, 1, 0, 0), 
                     ("%s (%.1f)" % (filteredSamples, 100.0*filteredTruePositives/filteredSamples), 2, 2, 0, 0))
                     , fileHandle, trailingLines=0)

writeEnd(fileHandle, "snpTable", "Total: the number of SNPs detected in each sample with respect to HG19. \
In brackets is the percentage confirmed true positives, already present in dbSNP.\
Total filtered: the number of SNPs detected in each sample not within 5 bps of an indel. \
Aggregate row: gives the total SNPs in human samples (excluding chimp).")
writeDocumentEnd(fileHandle)
fileHandle.close()