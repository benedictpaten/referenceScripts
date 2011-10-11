import xml.etree.ElementTree as ET
import sys

print sys.argv
i = ET.parse(sys.argv[1]).getroot()
l = []
l2 = set()
totalInsertion = 0
totalDeletion = 0
totalInsertionAndDeletion = 0
for j in i.findall("statsForSample"):
    if "panTro" not in j.attrib["sampleName"] and "reference" not in j.attrib["sampleName"]:
        for k in j.text.split("\n"):
            key = " ".join(k.split()[:6])
            if key not in l2:
                l2.add(key)
                l.append(k)
                if int(k.split()[5]) > 0:
                    totalInsertion += 1
                    if int(k.split()[13]) > 0:
                        totalDeletion += 1
                        totalInsertionAndDeletion += 1
                elif int(k.split()[13]) > 0:
                    totalDeletion += 1
                    
                
j = { "totalInsertion":str(totalInsertion),
      "totalDeletion":str(totalDeletion),
      "totalInsertionAndDeletion":str(totalInsertionAndDeletion),
      "sampleName":"aggregate",
      "referenceName":"hg19" }

ET.SubElement(i, "statsForSample", attrib=j).text = "\n".join(l)
fH = open(sys.argv[1][:-4] + "_withAggregates.xml", 'w')
ET.ElementTree(i).write(fH)
fH.close()