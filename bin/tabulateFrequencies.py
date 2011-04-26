import sys
import xml.etree.ElementTree as ET

"""Tabulates the frequencies of distributions.
"""

outputFileHandle = open(sys.argv[1], 'w')
def fn(node, path):
    if len(path) > 0:
        return fn(node.find(path[0]), path[1:])
    return node.text.split()
contigLengthsForAlgorithms = [ [ int(i) for i in fn(ET.parse(statsFile).getroot(), sys.argv[2].split(".")) ] for statsFile in sys.argv[3::2] ]
contigLengthsSet = set()
for contigLengths in contigLengthsForAlgorithms:
    contigLengthsSet = contigLengthsSet.union(set(contigLengths))
if 0 in contigLengthsSet:
    contigLengthsSet.remove(0)
contigLengthsSet = list(contigLengthsSet)
contigLengthsSet.sort()
outputFileHandle.write(("\t"*(2*len(contigLengthsForAlgorithms))).join(("Counts", "Cummulative_counts", "Summed_Cumulative_Counts\n")))
outputFileHandle.write("\t\t".join([ "\t".join([ "quantity\t" + algorithmName for algorithmName in sys.argv[4::2] ]) for i in xrange(3) ]) + "\n")
pVals = [ 0 ] * len(contigLengthsForAlgorithms)
pSumVals = [ 0 ] * len(contigLengthsForAlgorithms)
for contigLength in contigLengthsSet:
    outputFileHandle.write("\t".join([ str(contigLength) + "\t" + str(contigLengths.count(contigLength)) for contigLengths in contigLengthsForAlgorithms ]) + "\t\t")
    pVals = [ contigLengthsForAlgorithms[i].count(contigLength) + pVals[i] for i in xrange(len(contigLengthsForAlgorithms)) ]
    outputFileHandle.write("\t".join([ str(contigLength) + "\t" + str(pVal) for pVal in pVals ]) + "\t\t")
    pSumVals = [ contigLengthsForAlgorithms[i].count(contigLength)*contigLength + pSumVals[i] for i in xrange(len(contigLengthsForAlgorithms)) ]
    outputFileHandle.write("\t".join([ str(contigLength) + "\t" + str(pVal) for pVal in pSumVals ]) + "\n")
outputFileHandle.close()