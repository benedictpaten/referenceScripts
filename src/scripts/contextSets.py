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
        self.mappedSides = [ self ]
       
    @staticmethod 
    def getReverseComplement(string):
        def fn(i):
            m = { 'A':'T', 'C':'G', 'G':'C', 'T':'A'}
            if i in m:
                return m[i]
            return i
        return "".join([ fn(i) for i in string[::-1] ])
    
    def nonNAdjacencies(self):
        return [ adjSide for adjSide, ns in list(self.adjacencies) if ns == '' ]
    
    def base(self):
        if self.orientation:
            return self.basePosition.base
        return self.getReverseComplement(self.basePosition.base)
    
    def enumerateThreads(self, fn):
        stack = [ (side, side.otherSide.base()) for side in self.mappedSides ] #Add preceding base
        while len(stack) > 0:
            side, prefix = stack.pop()
            while len(side.nonNAdjacencies()) == 1:
                adjSide = side.nonNAdjacencies()[0]
                prefix += adjSide.base()
                side = adjSide.otherSide
            if fn(prefix) or len(side.nonNAdjacencies()) == 0:
                continue
            assert len(side.nonNAdjacencies()) > 1
            for adjSide in side.nonNAdjacencies():
                stack.append((adjSide.otherSide, prefix + adjSide.base()))  

def hamDist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

class ContextSet:
    def __init__(self):
        self.minimalUniqueStrings = set()
    
    def addString(self, string):
        """Add a string to the context set"""
        for i in xrange(1, len(string)):
            prefix = string[:i]
            if prefix in self.minimalUniqueStrings:
                return
        self.minimalUniqueStrings.add(string)

    def prefixInContextSet(self, string):
        for uniqueString in self.minimalUniqueStrings:
            if len(uniqueString) <= len(string) and hamDist(string, uniqueString) <= 0:
                return True
        return False  
        
class SequenceGraph:
    def __init__(self, usePhasedContexts):
        self.id = 0
        self.sides = []
        self.contextSets = {}
        self.maxContextStringLength  = 0
        self.usePhasedContexts = usePhasedContexts
        if self.usePhasedContexts:
            self.mappedSequenceGraph = SequenceGraph(usePhasedContexts=False)
        else:
            self.mappedSequenceGraph = self
            
    def getUniquePrefix(self, string, mismatches):
        startPoints = sum([ [ (side, mappedSide, 0)  for mappedSide in side.mappedSides ] for side in self.sides ], [])
        index = 0
        
        while len(startPoints) > 0 and index < len(string):
            l = []
            for side, mappedSide, diff in startPoints:
                diff += mappedSide.base() != string[index]
                if diff <= mismatches:
                    l.append((side, mappedSide, diff))
            index += 1
            
            if len(set([ side for side, mappedSide, diff in l ])) == 1: #We have a unique match to one node
                return string[:index], [ side for side, mappedSide, diff in l ][0].otherSide
               
            startPoints = []
            for side, mappedSide, diff in l:
                for adjMappedSide in mappedSide.otherSide.nonNAdjacencies():
                    startPoints.append((side, adjMappedSide, diff)) 
                    
        return (None, None)
        
    def computeContextSets(self, minContextLength, maxContextLength):
        """This function should be called on a graph before matching or printing is done"""
        self.contextSets = {}
        
        for side in self.sides:
            #Enumerate the threads
            contextSet = ContextSet()
            def fn(string):
                uniquePrefix, otherSide = self.getUniquePrefix(string, 0)
                if uniquePrefix == None:
                    return len(string) > maxContextLength #
                assert otherSide == side
                if len(uniquePrefix) < minContextLength:
                    if len(string) >= minContextLength:
                        contextSet.addString(string[:minContextLength])
                    else:
                        return False
                else:
                    contextSet.addString(uniquePrefix)
                return True #Return true stops the traversal on that thread
            side.enumerateThreads(fn)
            self.contextSets[side] = contextSet
            #Set max
            if len(contextSet.minimalUniqueStrings) > 0:
                maxContextStringLength = max([len(i) for i in contextSet.minimalUniqueStrings ])
                if maxContextStringLength > self.maxContextStringLength:
                    self.maxContextSetLength = maxContextStringLength   
        
    def getMatch(self, side, mismatches):
        assert side not in self.sides
        matches = set()
        def fn(string):
            uniquePrefix, otherSide = self.getUniquePrefix(string, mismatches)
            if uniquePrefix != None:
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
            for adjSide, ns in side2.adjacencies:
                assert (side2, ns) in adjSide.adjacencies
                adjSide.adjacencies.remove((side2, ns))
                assert (side2, ns) not in adjSide.adjacencies
                adjSide.adjacencies.add((side1, ns))
                side1.adjacencies.add((adjSide, ns))
        fn(side1, side2)
        fn(side1.otherSide, side2.otherSide)
        
        if self.usePhasedContexts:
            side1.mappedSides = side1.mappedSides + side2.mappedSides
            side1.otherSide.mappedSides = side1.otherSide.mappedSides + side2.otherSide.mappedSides
                
    def positiveSides(self):
        return [ side for side in self.sides if side.orientation ]
    
    def mergeSequenceGraphs(self, sG2):
        self.sides += sG2.sides
        if self.usePhasedContexts:
            assert sG2.usePhasedContexts
            self.mappedSequenceGraph.mergeSequenceGraphs(sG2.mappedSequenceGraph)
        
    def addString(self, string):
        pSide = None
        phasedPSide = None
        def tokenise(string):
            string = string.upper()
            ns = ""
            for i in string:
                if i != 'N':
                    yield i, ns
                    ns = ""
                else:
                    ns += 'N'
        for base, ns in tokenise(string):            
            def makeLinkedSides(base, ns, pSide, sides):
                bP = BasePosition(0, base)
                leftSide = Side(bP, 1)
                rightSide = Side(bP, 0)
                leftSide.otherSide = rightSide
                rightSide.otherSide = leftSide
                if pSide != None:
                    leftSide.adjacencies.add((pSide, ns))
                    pSide.adjacencies.add((leftSide, ns))
                sides.append(leftSide)
                sides.append(rightSide)
                return leftSide, rightSide
                    
            leftSide, rightSide = makeLinkedSides(base, ns, pSide, self.sides)
            pSide = rightSide
    
            if self.usePhasedContexts:
                mappedLeftSide, mappedRightSide = makeLinkedSides(base, ns, phasedPSide, self.mappedSequenceGraph.sides)
                leftSide.mappedSides = [ mappedLeftSide ]
                rightSide.mappedSides = [ mappedRightSide ]
                phasedPSide = mappedRightSide
                
        self.renumber(0) #Make sure everyone has a decent id.
    
    def renumber(self, startID=0, prefix=""):
        for side in self.positiveSides():
            if prefix != "":
                side.basePosition.id = str(prefix) + "_" + str(startID)
            else:
                side.basePosition.id = str(startID)
            startID += 1
        if self.usePhasedContexts:
            self.mappedSequenceGraph.renumber(startID=startID, prefix=prefix)
    
    def printDotFile(self, graphVizFileHandle, showContextSets, showIDs, number, label):
        graphVizFileHandle.write('subgraph cluster_%s {\nstyle=filled;\ncolor=lightgrey;label = "%s";\n' % (number,label))
        #Add nodes
        for side in self.positiveSides():
            #Graph vis
            if showContextSets:
                leftContextString = ",".join([ Side.getReverseComplement(i[1:]) for i in self.contextSets[side].minimalUniqueStrings ])
                if len(self.contextSets[side].minimalUniqueStrings) == 0:
                    leftContextString = "None"
                rightContextString = ",".join([ i[1:] for i in self.contextSets[side.otherSide].minimalUniqueStrings ])
                if len(self.contextSets[side.otherSide].minimalUniqueStrings) == 0:
                    rightContextString = "None"
                addNodeToGraph(nodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, 
                               shape="record", label="ID=%s | L=%s | %s | R=%s" % (side.basePosition.id, 
                                                                                       leftContextString, side.basePosition.base, rightContextString))
            elif showIDs:
                addNodeToGraph(nodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, shape="record", label="ID=%s | %s" % (side.basePosition.id, side.basePosition.base))
            else:
                addNodeToGraph(nodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, shape="record", label="%s" % (side.basePosition.base))
        #Add edges
        seen = set()
        for side in self.sides:
            for adjSide, ns in side.adjacencies:
                if not (adjSide, side, ns) in seen:
                    assert (side, adjSide, ns) not in seen
                    def arrowShape(side):
                        if side.orientation:
                            return "normal"
                        return "crow"
                    addEdgeToGraph(parentNodeName=side.basePosition.getDotNodeName(), 
                                   childNodeName=adjSide.basePosition.getDotNodeName(), 
                                   graphFileHandle=graphVizFileHandle, colour="black", #weight="1", 
                                   dir="both, arrowtail=%s, arrowhead=%s" % 
                                   (arrowShape(side), arrowShape(adjSide)), style="solid", length="10", label=ns)
                    seen.add((side, adjSide, ns))
        graphVizFileHandle.write("}\n")

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
    
    parser.add_option("--showIDs", dest="showIDs", type="string",
                     help="Show the IDs for the selected graphs (enumerated starting from 0)",
                     default="")
    
    parser.add_option("--mergeContigs", dest="mergeContigs", type="string",
                     help="Merge the contig in the selected graphs (enumerated starting from 0)",
                     default="")
    
    parser.add_option("--mismatches", dest="mismatches", type="int", 
                     help="Number of mismatches to allow",
                     default=0)
    
    parser.add_option("--minContextLength", dest="minContextLength", type="int", 
                     help="Minimum length of a string in a context set",
                     default=0)
    
    parser.add_option("--maxContextLength", dest="maxContextLength", type="int", 
                     help="Maximum length of a string in a context set",
                     default=10)
    
    parser.add_option("--showDoubleMaps", dest="showDoubleMaps", action="store_true",
                     help="Show doubly mapped match edges",
                     default=False)
    
    parser.add_option("--showOnlyLowestMaps", dest="showOnlyLowestMaps", action="store_true",
                     help="Show maps to one level in the hierarchy",
                     default=False)
    
    parser.add_option("--usePhasedContexts", dest="usePhasedContexts", action="store_true",
                     help="For each graph only allow context sets to be defined by the underlying sequences that serve as input to construct the sequence graph",
                     default=False)
    
    options, args = parser.parse_args()
    
    if len(args) == 0:
        parser.print_help()
        return 1
    
    mergeContigs = [ int(i) for i in options.mergeContigs.split() ]   
    
    #First create the sequence graphs for each input graph
    sequenceGraphs = []
    
    for index in xrange(len(args)):
        print "Processing sequence graph", index, options.mismatches
        assembly = args[index]
        sG = SequenceGraph(options.usePhasedContexts)
        sequenceGraphs.append(sG)
        for string in assembly.split():
            print "Adding string:", string, " of assembly:", index
            if index in mergeContigs: #Merge the new contig into the previous contigs
                sG2 = SequenceGraph(options.usePhasedContexts)
                sG2.addString(string)
                matches = []
                for side in sG2.positiveSides():
                    leftMatch = sG.getMatch(side, mismatches=options.mismatches)
                    rightMatch = sG.getMatch(side.otherSide, mismatches=options.mismatches)
                    def fn(side):
                        if side is None:
                            return "None"
                        return "%s_%s" % (side.basePosition.id, side.orientation)
                    if leftMatch != None:
                        if rightMatch == None or leftMatch.otherSide == rightMatch:
                            matches.append((leftMatch, side))
                    elif rightMatch != None:
                        matches.append((rightMatch.otherSide, side))
                sG.mergeSequenceGraphs(sG2)
                for targetSide, inputSide in matches:
                    sG.merge(targetSide, inputSide)
            else:
                sG.addString(string)
            print "Graph now has %i nodes" % len(sG.sides)
     
    #Now reindex them and print them
    showContextSets = [ int(i) for i in options.showContextSets.split() ]  
    showIDs = [ int(i) for i in options.showIDs.split() ]  
    graphVizFileHandle = open(options.graphVizFile, 'w')      
    setupGraphFile(graphVizFileHandle)
    graphVizFileHandle.write("splines=false;\n")   
    graphVizFileHandle.write("rankdir=LR;\n")   
    
    #i = 0
    for index in xrange(len(sequenceGraphs)):
        sG = sequenceGraphs[index]
        sG.renumber(prefix=index)
        print "Renumbering graph %i with %i sides" % (index, len(sG.positiveSides()))
        if showContextSets:
            if options.mismatches == 0:
                print "Calculating context sets for graph %i with %i sides" % (index, len(sG.positiveSides()))
                sG.computeContextSets(minContextLength=options.minContextLength, maxContextLength=options.maxContextLength)
            else:
                print "Can't display context sets with mismatches, they end up really big"
        #i += len(sG.positiveSides())
        sG.printDotFile(graphVizFileHandle, index in showContextSets and options.mismatches == 0, index in showIDs, index, "Sequence Graph %s" % index)
        
    #Now print the matching edges between the graphs
    for index in xrange(1, len(sequenceGraphs)):
        sGInput = sequenceGraphs[index]
        for side in sGInput.positiveSides():
            haveMatched = False
            for pIndex in xrange(index-1, -1, -1):
                sGTarget = sequenceGraphs[pIndex]
                
                leftMatch = sGTarget.getMatch(side, mismatches=options.mismatches)
                rightMatch = sGTarget.getMatch(side.otherSide, mismatches=options.mismatches)
                
                def addMatchEdge(colour, label, matchingSide):
                    if not options.showOnlyLowestMaps or not haveMatched:
                        addEdgeToGraph(parentNodeName=matchingSide.basePosition.getDotNodeName(), 
                                       childNodeName=side.basePosition.getDotNodeName(), graphFileHandle=graphVizFileHandle, colour=colour, 
                                       weight="100", label=label, dir="both, arrowtail=normal, arrowhead=none", style="solid", length="1")
                if leftMatch != None:
                    if rightMatch != None:
                        if leftMatch.otherSide == rightMatch: 
                            addMatchEdge("red", "B", leftMatch)
                            haveMatched = True
                        else:
                            if options.showDoubleMaps:
                                addMatchEdge("orange", "L", leftMatch)
                                addMatchEdge("orange", "R", rightMatch)
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
    