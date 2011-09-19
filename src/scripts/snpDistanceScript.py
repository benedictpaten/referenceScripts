import sys
import xml.etree.ElementTree as ET

samples = {}
#NODE_3393_length_512_cov_4_662109.NA.536.0.536.1 509 c reference.26 520 t
for i in ET.parse(sys.argv[1]).getroot().findall("statsForSample"):
    if i.attrib["sampleName"] not in { "ROOT", "" }:
        samples[i.attrib["sampleName"]] = set([ (j[3], int(j[4]), j[5].upper()) for j in \
        [ j.split() for j in i.text.split("\n")  if len(j.split()) == 6 ] ])        

samplesList = samples.keys()
samplesList.sort()

for i in xrange(len(samplesList)):
    sample1 = samplesList[i]
    for j in xrange(i+1,len(samplesList)):
        sample2 = samplesList[j]
        commonSnps = samples[sample1].intersection(samples[sample2])
        allSnps = samples[sample1].union(samples[sample2])
        print "\t".join([ sample1, sample2, str(len(commonSnps)), str(len(allSnps)), str(len(samples[sample1])), str(len(samples[sample2])), str(float(len(commonSnps)) / len(allSnps)) ])
