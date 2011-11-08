import sys
import xml.etree.ElementTree as ET

samples = 0
correct = 0
aligned = 0
sampleNo = 0
for i in ET.parse(sys.argv[1]).getroot().findall("statsForSample"):
    if i.attrib["sampleName"] not in ("reference", "hg19", "panTro3", ""):
        samples += int(i.attrib["totalSamples"])
        correct += int(i.attrib["totalCorrect"])
        aligned += int(i.attrib["totalAligned"])
        sampleNo += 1

print "samples", samples, "correct", correct, "aligned", aligned, "sampleNumber", sampleNo