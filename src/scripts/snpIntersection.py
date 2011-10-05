import xml.etree.ElementTree as ET
import sys

print sys.argv
i = ET.parse(sys.argv[1]).getroot()
l = []
l2 = set()
for j in i.findall("statsForSample"):
    for k in j.text.split("\n"):
        key = " ".join(k.split()[2:])
        if key not in l2:
            l2.add(key)
            l.append(k)
j = i.find("statsForSample").attrib.copy()
j["substitutionNumber"] = str(len(l))
j["sampleNumber"] = "NaN"
j["substitutionRate"] = "NaN"
j["sampleName"] = "aggregate"
ET.SubElement(i, "statsForSample", attrib=j).text = "\n".join(l)
fH = open(sys.argv[1][:-4] + "_withAggregates.xml", 'w')
ET.ElementTree(i).write(fH)
fH.close()