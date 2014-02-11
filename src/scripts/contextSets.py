import sys
from sonLib.bioio import addNodeToGraph, addEdgeToGraph, setupGraphFile, finishGraphFile
from optparse import OptionParser
import random

"""Script to calculate and display reference genome hierarchies
"""

class BasePosition:
    def __init__(self, id, base):
        self.base = base
        self.id = id
        
    def getDotNodeName(self):
        return "n%sn" % self.id

class Side:
    def __init__(self, basePosition, orientation):
        self.adjacencies = set()
        self.otherSide = None
        self.basePosition = basePosition
        self.orientation = orientation #A True orientation means on the left, otherwise on the right, by convention
       
    @staticmethod 
    def getReverseComplement(string):
        def fn(i):
            m = { 'A':'T', 'C':'G', 'G':'C', 'T':'A'}
            if i in m:
                return m[i]
            return i
        return "".join([ fn(i) for i in string[::-1] ])
    
    def base(self):
        if self.orientation:
            return self.basePosition.base
        return self.getReverseComplement(self.basePosition.base)
    
    def enumerateThreads(self, fn):
        stack = [ (self, self.otherSide.base()) ] #Add preceding base
        while len(stack) > 0:
            side, prefix = stack.pop()
            while len(side.adjacencies) == 1:
                adjSide = list(side.adjacencies)[0]
                prefix += adjSide.base()
                side = adjSide.otherSide
            if fn(prefix) or len(side.adjacencies) == 0:
                continue
            assert len(side.adjacencies) > 1
            for adjSide in side.adjacencies:
                stack.append((adjSide.otherSide, prefix + adjSide.base()))  

class ContextSet:
    def __init__(self, mismatches):
        self.mismatches = mismatches
        self.minimalUniqueStrings = set()
    
    def addString(self, string):
        """Add a string to the context set"""
        for i in xrange(1, len(string)):
            prefix = string[:i]
            if prefix in self.minimalUniqueStrings:
                return
        self.minimalUniqueStrings.add(string)
    
    @staticmethod
    def hamDist(str1, str2):
        """Count the # of differences between equal length strings str1 and str2"""
        diffs = 0
        for ch1, ch2 in zip(str1, str2):
            if ch1 != ch2:
                diffs += 1
        return diffs

    def prefixInContextSet(self, string):
        for i in xrange(1, len(string)):
            prefix = string[:i]
            for uniqueString in self.minimalUniqueStrings:
                if self.hamDist(prefix, uniqueString) <= self.mismatches:
                    return True
        return False  
        
class SequenceGraph:
    def __init__(self):
        self.id = 0
        self.sides = []
        self.contextSets = {}
        self.maxContextStringLength  = 0
        
    def computeContextSets(self, mismatches, minContextLength):
        """This function should be called on a graph before matching or printing is done"""
        self.contextSets = {}
        
        def getUniquePrefixWithRespectToGraph(string):
            startPoints = [ (side, 0) for side in self.sides ]
            index = 0
            while len(startPoints) > 1 and index < len(string):
                print "Recursion", len(startPoints), index, string[index]
                l = []
                for side, diff in startPoints:
                    diff += side.base() != string[index]
                    if diff <= mismatches:
                        for adjSide in side.otherSide.adjacencies:
                            l.append((adjSide, diff)) 
                startPoints = l
                index += 1
            print "End Recursion", len(startPoints), index
            if len(startPoints) == 1:
                if index < minContextLength:
                    if len(string) < minContextLength:
                        return None
                    index = minContextLength
                return string[:index]
            return None
        
        for side in self.sides:
            #Enumerate the threads
            contextSet = ContextSet(mismatches)
            def fn(string):
                uniquePrefix = getUniquePrefixWithRespectToGraph(string)
                print "For side ", side.basePosition.id, side.orientation, " we got the string ", string, " and the prefix ", uniquePrefix
                if uniquePrefix == None:
                    return False
                contextSet.addString(uniquePrefix)
                return True #Return true stops the traversal on that thread
            side.enumerateThreads(fn)
            self.contextSets[side] = contextSet
            #Set max
            if len(contextSet.minimalUniqueStrings) > 0:
                maxContextStringLength = max([len(i) for i in contextSet.minimalUniqueStrings ])
                if maxContextStringLength > self.maxContextStringLength:
                    self.maxContextSetLength = maxContextStringLength   
        
    def getMatch(self, side):
        matches = set()
        def fn(string):
            for otherSide in self.sides:
                if self.contextSets[otherSide].prefixInContextSet(string):
                    matches.add(otherSide)
                    return True
            return len(string) > self.maxContextStringLength
        side.enumerateThreads(fn)
        if len(matches) != 1:
            return None
        return matches.pop()
    
    def merge(self, side1, side2):
        self.sides.remove(side2)
        self.sides.remove(side2.otherSide)
        def fn(side1, side2):
            for adjSide in side2.adjacencies:
                assert side2 in adjSide.adjacencies
                adjSide.adjacencies.remove(side2)
                adjSide.adjacencies.add(side1)
                side1.adjacencies.add(adjSide)
        fn(side1, side2)
        fn(side1.otherSide, side2.otherSide)
                
    def positiveSides(self):
        return [ side for side in self.sides if side.orientation ]
    
    def mergeSequenceGraphs(self, sG2):
        self.sides += sG2.sides
        
    def addString(self, string):
        pSide = None
        for base in string:
            bP = BasePosition(0, base)
            leftSide = Side(bP, 1)
            rightSide = Side(bP, 0)
            leftSide.otherSide = rightSide
            rightSide.otherSide = leftSide
            if pSide != None:
                leftSide.adjacencies.add(pSide)
                pSide.adjacencies.add(leftSide)
            pSide = rightSide
            self.sides.append(leftSide)
            self.sides.append(rightSide)
        self.renumber(0) #Make sure everyone has a decent id.
    
    def renumber(self, startID):
        for side in self.positiveSides():
            side.basePosition.id = startID
            startID += 1
    
    def printDotFile(self, graphVizFileHandle, showContextSets):
        #Add nodes
        for side in self.positiveSides():
            #Graph vis
            if showContextSets:
                leftContextString = " ".join([ Side.getReverseComplement(i[1:]) for i in self.contextSets[side].minimalUniqueStrings ])
                if len(self.contextSets[side].minimalUniqueStrings) == 0:
                    leftContextString = "None"
                rightContextString = " ".join([ i[1:] for i in self.contextSets[side.otherSide].minimalUniqueStrings ])
                if len(self.contextSets[side.otherSide].minimalUniqueStrings) == 0:
                    rightContextString = "None"
                addNodeToGraph(nodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, 
                               shape="record", label="{ ID=%i | L=%s | %s | R=%s }" % (side.basePosition.id, 
                                                                                       leftContextString, side.basePosition.base, rightContextString))
            else:
                addNodeToGraph(nodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, shape="record", label="{ ID=%i | %s }" % (side.basePosition.id, side.basePosition.base))
        #Add edges
        seen = set()
        for side in self.sides:
            for adjSide in side.adjacencies:
                if not (adjSide, side) in seen:
                    assert (side, adjSide) not in seen
                    def arrowShape(side):
                        if side.orientation:
                            return "normal"
                        return "inv"
                    addEdgeToGraph(parentNodeName=side.basePosition.getDotNodeName(), 
                                   childNodeName=adjSide.basePosition.getDotNodeName(), 
                                   graphFileHandle=graphVizFileHandle, colour="black", #weight="1", 
                                   dir="both, arrowtail=%s, arrowhead=%s" % 
                                   (arrowShape(side), arrowShape(adjSide)), style="solid", length="10")
                    seen.add((side, adjSide))

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
    
    parser.add_option("--showDoubleMaps", dest="showDoubleMaps", action="store_true",
                     help="Show doubly mapped match edges",
                     default=False)
    
    parser.add_option("--showOnlyLowestMaps", dest="showOnlyLowestMaps", action="store_true",
                     help="Show maps to one level in the hierarchy",
                     default=False)
    
    options, args = parser.parse_args()
    
    if len(args) == 0:
        parser.print_help()
        return 1
    
    mergeContigs = [ int(i) for i in options.mergeContigs.split() ]   
    
    #First create the sequence graphs for each input graph
    sequenceGraphs = []
    
    for index in xrange(len(args)):
        print "Processing sequence graph", index
        assembly = args[index]
        sG = SequenceGraph()
        sequenceGraphs.append(sG)
        for string in assembly.split():
            print "Adding string:", string, " of assembly:", index
            if index in mergeContigs: #Merge the new contig into the previous contigs
                sG2 = SequenceGraph()
                sG2.addString(string)
                matches = []
                for side in sG2.positiveSides():
                    leftMatch = sG.getMatch(side)
                    rightMatch = sG.getMatch(side.otherSide)
                    if leftMatch != None:
                        if rightMatch == None or leftMatch.otherSide == rightMatch:
                            print "Got a left match", leftMatch, rightMatch
                            matches.append((leftMatch, side))
                    elif rightMatch != None:
                        print "Got a right match"
                        matches.append((rightMatch.otherSide, side))
                sG.mergeSequenceGraphs(sG2)
                print "We found %i merges" % len(matches)
                for targetSide, inputSide in matches:
                    sG.merge(targetSide, inputSide)
            else:
                sG.addString(string)
            print "Graph now has %i nodes" % len(sG.sides)
            sG.renumber(0)
            sG.computeContextSets(options.mismatches, options.minContextLength)
     
    #Now reindex them and print them
    showContextSets = [ int(i) for i in options.showContextSets.split() ]  
    graphVizFileHandle = open(options.graphVizFile, 'w')           
    setupGraphFile(graphVizFileHandle)
    
    i = 0
    for index in xrange(len(sequenceGraphs)):
        sG = sequenceGraphs[index]
        sG.renumber(i)
        print "Renumbering graph %i with %i sides" % (index, len(sG.positiveSides()))
        i += len(sG.positiveSides())
        sG.printDotFile(graphVizFileHandle, index in showContextSets)
        
    #Now print the matching edges between the graphs
    for index in xrange(1, len(sequenceGraphs)):
        sGInput = sequenceGraphs[index]
        for side in sGInput.positiveSides():
            haveMatched = False
            for pIndex in xrange(index-1, -1, -1):
                sGTarget = sequenceGraphs[pIndex]
                
                leftMatch = sGTarget.getMatch(side)
                rightMatch = sGTarget.getMatch(side.otherSide)
                
                def addMatchEdge(colour, label, matchingSide):
                    if not options.showOnlyLowestMaps or not haveMatched:
                        addEdgeToGraph(parentNodeName=matchingSide.basePosition.getDotNodeName(), 
                                       childNodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, colour=colour, weight="100", dir="%s, arrowtail=inv, arrowhead=normal" % label, style="solid", length="1")
                if leftMatch != None:
                    if rightMatch != None:
                        if leftMatch.otherSide == rightMatch: 
                            addMatchEdge("red", "B", leftMatch)
                            haveMatched = True
                        else:
                            if options.showDoubleMaps:
                                addMatchEdge("grey", "L", leftMatch)
                                addMatchEdge("grey", "R", rightMatch)
                    else:
                        addMatchEdge("blue", "L", leftMatch)
                        haveMatched = True
                elif rightMatch != None:
                    addMatchEdge("blue", "R", rightMatch)
                    haveMatched = True
        
    finishGraphFile(graphVizFileHandle)   
    graphVizFileHandle.close()  

if __name__ == '__main__':
    exit(main())  
    