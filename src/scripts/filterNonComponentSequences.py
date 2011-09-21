import sys
import xml.etree.ElementTree as ET
from sonLib.bioio import fastaRead, fastaWrite
i = set([ i for i in ET.parse(sys.argv[1]).getroot().text.split() ])
fH = open(sys.argv[3], 'w')
for name, sequence in fastaRead(open(sys.argv[2], 'r')):
        if name not in i:
            fastaWrite(fH, name, sequence)
fH.close()