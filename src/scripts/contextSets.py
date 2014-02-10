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

class BasePosition:
    def __init__(self, id, base):
        self.base = base
        self.id = id
        
    def getDotNodeName(self):
        return "n%sn" % self.bP.id

class Side:
    def __init__(self, bP, orientation):
        self.adj = []
        self.otherSide = None
        self.bP = bP
        self.orientation = orientation
    
    def base(self):
        if self.orientation:
            return self.bP.base
        return getReverseComplement(self.bP.base)
    
    def enumerateThreads(self):
    
def getRightContextString(sG, string, mismatches, minContextLength):
    startPoints = sG.sides
    while len(startPoints) > 1 and index < len(string):
        l = []
        for side, diff in startPoints:
            diff += side.base() == string[index]
            if diff < mismatches:
                for adjSide in side.otherSide.adj:
                    l.append((adjSide, diff)) 
        startPoints = l
        index += 1
    if len(startPoints) == 1:
        if index < minContextLength:
            if len(string) < minContextLength:
                return None
            index = minContextLength
        return string[:index]

def getRightContextStrings(sG, side, mismatches, minContextLength):
    #Enumerate the threads
    return []
        
class SequenceGraph:
    def __init__(self):
        self.id = 0
        self.sides = []
        self.contextSets = {}
    
    def merge(self, side1, side2):
        assert side1 in self.sides and side2 in self.sides
        self.sides.remove(side2)
        #
        
    def _computeContextSets(self):
        pass
    
    def positiveSides(self):
        return [ side for side in self.sides if side.orientation ]
        
    def getMatch(self, side):
        pass
        
    def addString(self, string):
        return None #Return the end of the string.
    
    def renumber(self, startID):
        pass #Give the sides IDs in accordance with their
    
    def printDotFile(self, fileHandle, showContextSets):
        for side in self.positiveSides():
            #Graph vis
            if options.showContextSets:
                addNodeToGraph(nodeName=refNodeName(i, j), graphFileHandle=graphVizFileHandle, shape="record", label="{ ID=%i | L=%s | %s | R=%s }" % (side.basePosition.getDotNodeName(), 
                                                                                                                                                       " ".join(side.leftContextStrings), side.basePosition.base, " ".join(side.rightContextStrings)))
            else:
                addNodeToGraph(nodeName=refNodeName(i, j), graphFileHandle=graphVizFileHandle, shape="record", label="{ ID=%i | %s }" % (side.basePosition.getDotNodeName(), side.basePosition.base))
        
        seen = set()
        for side in self.sides:
            for adjSide in side.adjacencies:
                if not (adjSide, side) in seen:
                    #Add in the orientations to adjacencies
                    addEdgeToGraph(parentNodeName=adjSide.basePosition.getDotNodeName(), 
                                   childNodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, colour="black", weight="100", dir="both, arrowtail=inv, arrowhead=normal", style="solid", length="1")
                    seen.add((side, adjSide))

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
    
    parser.add_option("--showContextSets", dest="showContextSets", type="string",
                     help="Show the context sets for the selected graphs (enumerated starting from 0)",
                     default="")
    
    parser.add_option("--mergeContigs", dest="mergeContigs", type="string",
                     help="Merge the contig in the selected graphs (enumerated starting from 0)",
                     default="")
    
    parser.add_option("--mismatches", dest="mismatches", type="int", 
                     help="Number of mismatches to allow",
                     default=0)
    
    parser.add_option("--minContextLength", dest="minContextLength", type="int", 
                     help="Minimum length of a context set",
                     default=0)
    
    parser.add_option("--showDoubleMaps", dest="showDoubleMaps", action="store_true"
                     help="Show doubly mapped match edges",
                     default=False)
    
    parser.add_option("--showOnlyLowestMaps", dest="showOnlyLowestMaps", action="store_true"
                     help="Show maps to one level in the hierarchy",
                     default=False)
    
    options, args = parser.parse_args()
    
    if len(args) == 0 or len(args) > 2:
        parser.print_help()
        return 1
    
    mergeContigs = [ int(i) for i in options.mergeContigs.split() ]   
    
    #First create the sequence graphs for each input graph
    sequenceGraphs = []
    
    for index in xrange(len(args)):
        assembly = args[index]
        sG = SequenceGraph()
        sequenceGraphs.append(sG)
        for string in assembly.split():
            if index in mergeContigs: #Merge the new contig into the previous contigs
                sG2 = SequenceGraph()
                sG2.addString(string)
                matches = []
                for side in sG2.positiveSides():
                    leftMatch = sGTarget.getMatch(side)
                    rightMatch = sGTarget.getMatch(side.otherSide)
                    if leftMatch != None:
                        if rightMatch == None || leftMatch.oppositeSide == rightMatch:
                            matches.append((side, leftSide))
                    elif rightMatch != None:
                            matches.append((side, leftSide))
                sG.addSequenceGraph(sG2)
                for targetSide, inputSide in matches:
                    sG.merge(targetSide, inputSide)
            else:
                sG.addString(string)
     
    
    #Now reindex them and print them
    showContextSets = [ int(i) for i in options.showContextSets.split() ]  
    graphVizFileHandle = open(options.graphVizFile, 'w')           
    
    setupGraphFile(graphVizFileHandle)
    
    i = 0
    for index in len(sequenceGraphs):
        sG = sequenceGraphs[index]
        sG.renumber(i)
        i += len(sG.positiveSides())
        sG.printDotFile(graphVizFileHandle, showContextSets[index])
        
    #Now print the matching edges between the graphs
    for index in xrange(1, len(sequenceGraphs)):
        sGInput = sequenceGraphs[i]
        for side in sGInput.positiveSides():
            haveMatched = False
            for pIndex in xrange(index-1, 0, -1):
                sGTarget = sequenceGraphs[pIndex]
                
                leftMatch = sGTarget.getMatch(side)
                rightMatch = sGTarget.getMatch(side.otherSide)
                
                def addMatchEdge(colour, label, matchSide):
                    if not showOnlyLowestMaps or not haveMatched:
                        addEdgeToGraph(parentNodeName=matchingSide.getDotName(), 
                                       childNodeName=side.getDotName(), graphFileHandle=graphVizFileHandle, colour=colour, weight="100", dir="%s, arrowtail=inv, arrowhead=normal" % label, style="solid", length="1")
                if leftMatch != None:
                    if rightMatch != None:
                        if leftMatch.oppositeSide == rightMatch: 
                            addMatchEdge("red", "B", leftMatch)
                            haveMatched = True
                        else:
                            if options.showDoubleMaps:
                                addMatchEdge("grey", "L", leftMatch)
                                addMatchEdge("grey", "R", rightMatch)
                    else:
                        addMatchEdge("blue", "L", leftMatch)
                        haveMatched = True
                else:
                    if rightMatch != None:
                        addMatchEdge("blue", "R", rightMatch)
                        haveMatched = True
            
        
    finishGraphFile(graphVizFileHandle)   
    graphVizFileHandle.close()  
    
    
    
    id = 0
    for j in xrange(len(genome)):
        string = genome[j].upper()
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
                addNodeToGraph(nodeName=refNodeName(i, j), graphFileHandle=graphVizFileHandle, shape="record", label="{ ID=%i | L=%s | %s | R=%s }" % (id, leftContextString, string[i], rightContextString))
            else:
                addNodeToGraph(nodeName=refNodeName(i, j), graphFileHandle=graphVizFileHandle, shape="record", label="{ ID=%i | %s }" % (id, string[i]))
            id += 1
            
            if i > 0:
                addEdgeToGraph(parentNodeName=refNodeName(i-1, j), 
                               childNodeName=refNodeName(i, j), graphFileHandle=graphVizFileHandle, colour="black", weight="100", dir="both, arrowtail=inv, arrowhead=normal", style="solid", length="1")
    
    """  
    if len(args) > 1:
        strings = args[1].split()
        for stringIndex in xrange(len(strings)):
            string = strings[stringIndex].upper()
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
                addNodeToGraph(nodeName=queryNodeName(stringIndex,i), graphFileHandle=graphVizFileHandle, label=string[i], shape="record")
                if i > 0:
                    addEdgeToGraph(parentNodeName=queryNodeName(stringIndex,i-1), 
                                   childNodeName=queryNodeName(stringIndex,i), graphFileHandle=graphVizFileHandle, weight="100", style="solid", dir="both, arrowtail=inv, arrowhead=normal", length="1")
                
                if rightContextSet != None and leftContextSet != None:
                    if leftContextSet == rightContextSet:
                        j, k = rightContextSet
                        addEdgeToGraph(parentNodeName=queryNodeName(stringIndex,i), 
                                   childNodeName=refNodeName(k, j), graphFileHandle=graphVizFileHandle, weight="0.0", colour="red", style="dotted", dir="forward", length="1", label="B")
                    else:
                        j, k = rightContextSet
                        addEdgeToGraph(parentNodeName=queryNodeName(stringIndex,i), 
                                   childNodeName=refNodeName(k, j), graphFileHandle=graphVizFileHandle, weight="0.0", colour="blue", style="dotted", dir="forward", length="1", label="R")
                        j, k = leftContextSet
                        addEdgeToGraph(parentNodeName=queryNodeName(stringIndex,i), 
                                   childNodeName=refNodeName(k, j), graphFileHandle=graphVizFileHandle, weight="0.0", colour="blue", style="dotted", dir="forward", length="1", label="L") 
                        
                elif rightContextSet != None:
                    j, k = rightContextSet
                    addEdgeToGraph(parentNodeName=queryNodeName(stringIndex,i), 
                                   childNodeName=refNodeName(k, j), graphFileHandle=graphVizFileHandle, weight="0.0", colour="black", style="dotted", dir="forward", length="1", label="R")
                    
                elif leftContextSet != None:
                    j, k = leftContextSet
                    addEdgeToGraph(parentNodeName=queryNodeName(stringIndex,i), 
                                   childNodeName=refNodeName(k, j) + ":bottom", graphFileHandle=graphVizFileHandle, weight="0.0", colour="black", style="dotted", dir="forward", length="1", label="L")
    """
        
    finishGraphFile(graphVizFileHandle)   
    graphVizFileHandle.close()  

if __name__ == '__main__':
    exit(main())  
    