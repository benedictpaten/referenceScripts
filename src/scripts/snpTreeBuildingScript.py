import sys
import xml.etree.ElementTree as ET
from sonLib.tree import BinaryTree
from sonLib.tree import njI
from sonLib.tree import upgmaI
from sonLib.tree import DistancePair
from sonLib.bioio import printBinaryTree

l = {}
def fn(eventName):
    if not l.has_key(eventName):
        l[eventName] = BinaryTree(0.0, False, None, None, eventName)
    return l[eventName]
distancePairs = [ DistancePair(float(i.attrib["substitutionRate"]), fn(i.attrib["eventName1"]), 1, fn(i.attrib["eventName2"]), 1) for i in ET.parse(sys.argv[1]).getroot().findall("distancesForSamples") ] 
distancePairs += [ DistancePair(i.distance, i.leaf2, 1, i.leaf1, 1) for i in distancePairs ]

print len(distancePairs), l
print "NJ", printBinaryTree(njI(distancePairs, len(l.keys())), includeDistances=True)
print "UPGMA", printBinaryTree(upgmaI(distancePairs, len(l.keys())), includeDistances=True)
