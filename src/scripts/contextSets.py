import sys
from sonLib.bioio import addNodeToGraph, addEdgeToGraph, setupGraphFile, finishGraphFile
from optparse import OptionParser

"""Script to calculate and display context sets
"""

def getReverseComplement(string):
    if string == None:
        return None
    def fn(i):
        m = { 'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        if i in m:
            return m[i]
        return i
    return "".join([ fn(i) for i in string[::-1] ])

def getSubstrings(string):
    return [ string[i:] for i in xrange(len(string))]

def getSubstringsOfGenome(genome):
    subStrings = []
    for string in genome:
        subStrings += getSubstrings(string) + getSubstrings(getReverseComplement(string))
    return subStrings

def getRightContextString(genome, string, mismatches, minContextLength):
    startPoints = zip(getSubstringsOfGenome(genome), [ 0 ]*len(getSubstringsOfGenome(genome)))
    index = 0
    while len(startPoints) > 1 and index < len(string):
        startPoints = [ (matchedString, diff + (matchedString[index] != string[index])) for matchedString, diff in startPoints if index < len(matchedString) and (matchedString[index] == string[index] or diff < mismatches) ]
        index += 1
    if len(startPoints) == 1:
        if index < minContextLength:
            if len(string) < minContextLength:
                return None
            index = minContextLength
        return string[:index]
    return None

def hamdist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

def getMatchingRightContextSet(rightContextSets, string, mismatches):
    for contextString in rightContextSets.keys():
        if len(string) >= len(contextString) and hamdist(contextString, string[:len(contextString)]) <= mismatches:
            return rightContextSets[contextString]
    return None

def main():
    ##########################################
    #Construct the arguments.
    ##########################################    
    
    usage = "usage: %prog [options] <query> <target> <output>\n\n" + \
            "    <genome>:  strings representing genome\n" + \
            "    <string>: string to match\n" 
    description = "Script to demonstrate context sets and unique matching"
    parser = OptionParser(usage=usage, description=description)

    #output stuff
    parser.add_option("--graphVizFile", dest="graphVizFile", type="string",
                     help="File to do dump graphViz representation in (dot format)",
                     default="/dev/null")
    
    parser.add_option("--showContextSets", dest="showContextSets",  action="store_true",
                     help="Show the context sets",
                     default=False)
    
    parser.add_option("--mismatches", dest="mismatches", type="int", 
                     help="Number of mismatches to allow",
                     default=0)
    
    parser.add_option("--minContextLength", dest="minContextLength", type="int", 
                     help="Minimum length of a context set",
                     default=0)
    
    
    options, args = parser.parse_args()
    
    if len(args) == 0 or len(args) > 2:
        parser.print_help()
        return 1
    
    graphVizFileHandle = open(options.graphVizFile, 'w')
    
    genome = args[0].split()
    
    rightContextSets = {}
    
    setupGraphFile(graphVizFileHandle)
    
    def refNodeName(i, j):
        return "n%ip%ipn" % (abs(i), abs(j))
    
    def queryNodeName(i):
        return "n%in" % i
    
    id = 0
    for j in xrange(len(genome)):
        string = genome[j]
        for i in xrange(len(string)):
            leftContextString = getReverseComplement(getRightContextString(genome, getReverseComplement(string[:(i+1)]), mismatches=options.mismatches, minContextLength=options.minContextLength))
            rightContextString = getRightContextString(genome, string[i:], mismatches=options.mismatches, minContextLength=options.minContextLength)
            
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
            if options.showContextSets:
                addNodeToGraph(nodeName=refNodeName(i, j), graphFileHandle=graphVizFileHandle, shape="record", label="{ ID=%i | 5\'=%s | %s | 3\'=%s }" % (id, leftContextString, string[i], rightContextString))
            else:
                addNodeToGraph(nodeName=refNodeName(i, j), graphFileHandle=graphVizFileHandle, shape="record", label="{ ID=%i | %s }" % (id, string[i]))
            id += 1
            
            if i > 0:
                addEdgeToGraph(parentNodeName=refNodeName(i-1, j), 
                               childNodeName=refNodeName(i, j), graphFileHandle=graphVizFileHandle, colour="black", weight="100", dir="both, arrowtail=inv, arrowhead=normal", style="solid", length="1")
      
    if len(args) > 1:
        string = args[1]
        for i in xrange(len(string)):
            leftString = getReverseComplement(string[:i+1])
            leftContextSet = getMatchingRightContextSet(rightContextSets, leftString, mismatches=options.mismatches)
            if leftContextSet != None:
                j, k = leftContextSet
                leftContextSet = j, -k
            
            rightString = string[i:]
            rightContextSet = getMatchingRightContextSet(rightContextSets, rightString, mismatches=options.mismatches)
            
            print "String-position", i, "left-mapping", leftContextSet, "right-mapping", rightContextSet
            
            #Graph vis
            addNodeToGraph(nodeName=queryNodeName(i), graphFileHandle=graphVizFileHandle, label=string[i], shape="record")
            if i > 0:
                addEdgeToGraph(parentNodeName=queryNodeName(i-1), 
                               childNodeName=queryNodeName(i), graphFileHandle=graphVizFileHandle, weight="100", style="solid", dir="both, arrowtail=inv, arrowhead=normal", length="1")
            
            if rightContextSet != None and leftContextSet != None:
                if leftContextSet == rightContextSet:
                    j, k = rightContextSet
                    addEdgeToGraph(parentNodeName=queryNodeName(i), 
                               childNodeName=refNodeName(k, j), graphFileHandle=graphVizFileHandle, weight="0.0", colour="red", style="dotted", dir="forward", length="1", label="B")
                else:
                    j, k = rightContextSet
                    addEdgeToGraph(parentNodeName=queryNodeName(i), 
                               childNodeName=refNodeName(k, j), graphFileHandle=graphVizFileHandle, weight="0.0", colour="blue", style="dotted", dir="forward", length="1", label="3")
                    j, k = leftContextSet
                    addEdgeToGraph(parentNodeName=queryNodeName(i), 
                               childNodeName=refNodeName(k, j), graphFileHandle=graphVizFileHandle, weight="0.0", colour="blue", style="dotted", dir="forward", length="1", label="5") 
                    
            elif rightContextSet != None:
                j, k = rightContextSet
                addEdgeToGraph(parentNodeName=queryNodeName(i), 
                               childNodeName=refNodeName(k, j), graphFileHandle=graphVizFileHandle, weight="0.0", colour="black", style="dotted", dir="forward", length="1", label="3")
                
            elif leftContextSet != None:
                j, k = leftContextSet
                addEdgeToGraph(parentNodeName=queryNodeName(i), 
                               childNodeName=refNodeName(k, j) + ":bottom", graphFileHandle=graphVizFileHandle, weight="0.0", colour="black", style="dotted", dir="forward", length="1", label="5")
    
    finishGraphFile(graphVizFileHandle)   
    graphVizFileHandle.close()  

if __name__ == '__main__':
    exit(main())  
    