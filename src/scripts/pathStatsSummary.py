import xml.etree.ElementTree as ET
import sys

def getSampleNames(fileName):
    sampleNames = set()
    pathStatsTag = ET.parse(fileName).getroot()
    for pathStatForSampleTag in pathStatsTag.findall("statsForSample"):
        sampleName = pathStatForSampleTag.attrib["sampleName"]
        if sampleName != "":
           sampleNames.add(sampleName)
    return sampleNames

sampleNames = getSampleNames(sys.argv[1])
for fileName in sys.argv[2:]:
    sampleNames = sampleNames.intersection(getSampleNames(fileName))

sampleNamesToComparisons = {}
for sampleName in sampleNames:
    sampleNamesToComparisons[sampleName] = []

for fileName in sys.argv[1:]:
    pathStatsTag = ET.parse(fileName).getroot()
    for pathStatForSampleTag in pathStatsTag.findall("statsForSample"):
        sampleName = pathStatForSampleTag.attrib["sampleName"]
        if sampleName in sampleNames:
            sampleNamesToComparisons[sampleName].append(pathStatForSampleTag)
        
print "args", sys.argv[1:]
        
statNames = ("blockN50", "contigPathN50", "scaffoldPathN50", "totalInsertion", "totalDeletion")
        
for sampleName in sampleNames:
    print sampleName + "\t" + "\t".join([ "\t".join([ i.attrib[statName] for i in sampleNamesToComparisons[sampleName] ]) for statName in statNames ])

print "aggregate" + "\t" + "\t".join([ "\t".join([ str(sum([ float(sampleNamesToComparisons[sampleName][i].attrib[statName]) for sampleName in sampleNames ])/len(sampleNames)) for i in xrange(len(sys.argv[1:])) ]) for statName in statNames ])