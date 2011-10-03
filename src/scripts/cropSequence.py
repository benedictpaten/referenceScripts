import sys
import xml.etree.ElementTree as ET
from sonLib.bioio import fastaRead, fastaWrite
node = ET.parse(sys.argv[1]).getroot()
fH = open(sys.argv[3], 'w')
seqs = [ i for i in fastaRead(open(sys.argv[2], 'r')) ]
assert(len(seqs) == 1)
for name, sequence in seqs:
    #>hg19.chr6.171115067.28377796.5150977.1
    i = name.split(".")
    j = int(node.attrib["minOtherReferenceCoordinate"])
    k = int(node.attrib["maxOtherReferenceCoordinate"])
    fastaWrite(fH, ".".join(i[0:3] + [ str(int(i[3]) + j), str(k - j)] + i[-1:]), sequence[j:k])
fH.close()