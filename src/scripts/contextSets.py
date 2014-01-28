import sys
from sonLib.bioio import addNodeToGraph, addEdgeToGraph, setupGraphFile, finishGraphFile

"""Script to calculate context sets
"""

def getReverseComplement(string):
    if string == None:
        return None
    return "".join([ { 'A':'T', 'C':'G', 'G':'C', 'T':'A'}[i] for i in string[::-1] ])

def getSubstrings(string):
    return [ string[i:] for i in xrange(len(string))]

def getSubstringsOfGenome(genome):
    subStrings = []
    for string in genome:
        subStrings += getSubstrings(string) + getSubstrings(getReverseComplement(string))
    return subStrings

def getContextString(genome, string):
    startPoints = getSubstringsOfGenome(genome)
    index = 0
    while len(startPoints) > 1 and index < len(string):
        startPoints = [ matchedString for matchedString in startPoints if index < len(matchedString) and matchedString[index] == string[index]]
        index += 1
    if len(startPoints) == 1:
        return string[:index]
    return None

fileHandle = open(sys.argv[3], 'w')

genome = sys.argv[1].split()

rightContextSets = {}

setupGraphFile(fileHandle)

def refNodeName(i, j):
    return "n%ip%ipn" % (abs(i), abs(j))

def queryNodeName(i):
    return "n%in" % i

for j in xrange(len(genome)):
    string = genome[j]
    for i in xrange(len(string)):
        leftContextString = getReverseComplement(getContextString(genome, getReverseComplement(string[:(i+1)])))
        rightContextString = getContextString(genome, string[i:])
        
        if leftContextString != None:
            assert getReverseComplement(leftContextString) not in rightContextSets
            rightContextSets[getReverseComplement(leftContextString)] = (j, -i)
            leftContextString = leftContextString[:-1]
        
        if rightContextString != None:
            assert rightContextString not in rightContextSets
            rightContextSets[rightContextString] = (j, i)
            rightContextString = rightContextString[1:]
        
        #getRightContextSet(genome, getReverseComplement(string[:-(i+1)]))
        print "Genome-position", i, "left-context-set", leftContextString, "right-context-set", rightContextString
        
        #Graph vis
        addNodeToGraph(nodeName=refNodeName(i, j), graphFileHandle=fileHandle, shape="square", label=string[i])
        if i > 0:
            addEdgeToGraph(parentNodeName=refNodeName(i-1, j), 
                           childNodeName=refNodeName(i, j), graphFileHandle=fileHandle, colour="green", weight="100")
  
def getMatchingContextSet(contextSets, string):
    for contextString in contextSets.keys():
        if len(string) >= len(contextString) and contextString == string[:len(contextString)]:
            return contextSets[contextString]
    return None
        

string = sys.argv[2]
for i in xrange(len(string)):
    leftString = getReverseComplement(string[:i+1])
    leftContextSet = getMatchingContextSet(rightContextSets, leftString)
    if leftContextSet != None:
        j, k = leftContextSet
        leftContextSet = j, -k
    
    rightString = string[i:]
    rightContextSet = getMatchingContextSet(rightContextSets, rightString)
    
    print "String-position", i, "left-mapping", leftContextSet, "right-mapping", rightContextSet
    
    #Graph vis
    addNodeToGraph(nodeName=queryNodeName(i), graphFileHandle=fileHandle, label=string[i])
    if i > 0:
        addEdgeToGraph(parentNodeName=queryNodeName(i-1), 
                       childNodeName=queryNodeName(i), graphFileHandle=fileHandle, weight="100")
        
    if rightContextSet != None:
        j, k = rightContextSet
        addEdgeToGraph(parentNodeName=queryNodeName(i), 
                       childNodeName=refNodeName(k, j), graphFileHandle=fileHandle, weight="0.0", colour="blue", dir="forward")
    
    if leftContextSet != None:
        j, k = leftContextSet
        addEdgeToGraph(parentNodeName=queryNodeName(i), 
                       childNodeName=refNodeName(k, j), graphFileHandle=fileHandle, weight="0.0", colour="red", dir="forward")

finishGraphFile(fileHandle)   
fileHandle.close()    
    