import sys
import xml.etree.ElementTree as ET

for i in [ i for i in ET.parse(sys.argv[1]).getroot().findall("statsForSample") if i.attrib["sampleName"] == sys.argv[2] ][0].text.split("\n"):
    i = i.split()
    if len(i) > 0:
        j = i[3].split(".")
        #NODE_1822_length_25843_cov_10_938745.NA.25867.0.606.1 585 C hg19.chr6.171115067.28468848.5047277.1 4908635 g
        print "\t".join(i[:3] + [ "%s.%s" % (j[0], j[1]), str(int(j[3]) + int(i[4])) ] + i[5:])
