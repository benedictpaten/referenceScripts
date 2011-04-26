import sys
import xml.etree.ElementTree as ET

"""Makes a scatter plot of X and Y
"""

outputFileHandle = open(sys.argv[1], 'w')
def fn(node, path):
    if len(path) > 0:
        return fn(node.find(path[0]), path[1:])
    return node.text.split()
statsNode = ET.parse(sys.argv[4]).getroot()
outputFileHandle.write("X\tY\t%s\t%s\n" % (sys.argv[2], sys.argv[3]))
X = [ int(i) for i in fn(statsNode, sys.argv[2].split(".")) ]
Y = [ int(i) for i in fn(statsNode, sys.argv[3].split(".")) ]
assert len(X) == len(Y)
for x, y in zip(X, Y):
    outputFileHandle.write("%s\t%s\t%f\n" % (x, y, float(x)/(y+0.0000000001)))
outputFileHandle.close()
