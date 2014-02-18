import sys
from sonLib.bioio import addNodeToGraph, addEdgeToGraph, setupGraphFile, finishGraphFile
from optparse import OptionParser
import random

"""Script to calculate and display reference genome hierarchies and mappings
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
            i = 0
            while len(side.nonNAdjacencies()) == 1 and i < 10:
                adjSide = side.nonNAdjacencies()[0]
                prefix += adjSide.base()
                side = adjSide.otherSide
                i += 1
            if fn(prefix) or len(side.nonNAdjacencies()) == 0:
                continue
            #assert len(side.nonNAdjacencies()) > 1
            for adjSide in side.nonNAdjacencies():
                stack.append((adjSide.otherSide, prefix + adjSide.base()))  

def hamDist(str1, str2):
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

def sharedPrefix(str1, str2):
    """Return shared prefix of the two strings"""
    for i in xrange(min(len(str1), len(str2))):
        if str1[i] != str2[i]:
            return str1[:i]
    return str1[:min(len(str1), len(str2))]

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
    
    def maxLength(self):
        return max([ 0 ] + [ len(i) for i in self.minimalUniqueStrings ])
    
    def maxSharedPrefixLength(self, otherContextSet):
        i = ""
        for uniqueString in self.minimalUniqueStrings:
            for otherUniqueString in otherContextSet.minimalUniqueStrings:
                j = sharedPrefix(uniqueString, otherUniqueString)
                if len(j) > len(i):
                    i = j
        return len(i) 
        
class SequenceGraph:
    def __init__(self, usePhasedContexts, label=""):
        self.id = 0
        self.sides = []
        self.contextSets = {}
        self.maxContextStringLength  = 0
        self.usePhasedContexts = usePhasedContexts
        if self.usePhasedContexts:
            self.mappedSequenceGraph = SequenceGraph(usePhasedContexts=False)
        else:
            self.mappedSequenceGraph = self
        self.label = label
            
    def getUniquePrefix(self, string, mismatches, dissimilarity=1):
        assert dissimilarity >= 1
        startPoints = sum([ [ (side, mappedSide, 0)  for mappedSide in side.mappedSides ] for side in self.sides ], [])
        index = 0
        
        while len(startPoints) > 0 and index < len(string):
            l = []
            for side, mappedSide, diff in startPoints:
                diff += mappedSide.base() != string[index]
                if diff <= mismatches + dissimilarity - 1:
                    l.append((side, mappedSide, diff))
            index += 1
            
            s2 = set([ side for side, mappedSide, diff in l if diff <= mismatches + dissimilarity - 1 ]) 
            s = set([ side for side, mappedSide, diff in l if diff <= mismatches ]) 
            if len(s2) == 1 and len(s) == 1: #We have a unique match to one node
                return string[:index], s.pop().otherSide
               
            startPoints = []
            for side, mappedSide, diff in l:
                for adjMappedSide in mappedSide.otherSide.nonNAdjacencies():
                    startPoints.append((side, adjMappedSide, diff)) 
                    
        return (None, None)
        
    def computeContextSets(self, minContextLength, maxContextLength, dissimilarity=1):
        """This function should be called on a graph before matching or printing is done"""
        self.contextSets = {}
        
        for side in self.sides:
            #Enumerate the threads
            contextSet = ContextSet()
            def fn(string):
                uniquePrefix, otherSide = self.getUniquePrefix(string, 0, dissimilarity=dissimilarity)
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
            if contextSet.maxLength() > self.maxContextStringLength:
                self.maxContextSetLength = contextSet.maxLength()  
        
    def getMatch(self, side, mismatches, maxMatchLength=100, dissimilarity=1):
        assert side not in self.sides
        matches = set()
        def fn(string):
            uniquePrefix, otherSide = self.getUniquePrefix(string, mismatches, dissimilarity=dissimilarity)
            if uniquePrefix != None:
                matches.add(otherSide)
                return True
            return len(string) > maxMatchLength
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
        if side1.basePosition.id > side2.basePosition.id:
            side1.basePosition.id = side2.basePosition.id
        
        if self.usePhasedContexts:
            side1.mappedSides = side1.mappedSides + side2.mappedSides
            side1.otherSide.mappedSides = side1.otherSide.mappedSides + side2.otherSide.mappedSides
                
    def positiveSides(self):
        return [ side for side in self.sides if side.orientation ]
    
    def mergeSequenceGraphs(self, sG2):
        sG2.renumber(startID = len(self.positiveSides()))
        self.sides += sG2.sides
        if self.usePhasedContexts:
            assert sG2.usePhasedContexts
            self.mappedSequenceGraph.mergeSequenceGraphs(sG2.mappedSequenceGraph)
        
    def addString(self, string):
        pSide = None
        phasedPSide = None
        firstSide = None
        circular = string[-1:] == '!'
        if circular:
            string = string[:-1]
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
                bP = BasePosition(len(sides)/2, base)
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
            if firstSide == None:
                firstSide = leftSide
            pSide = rightSide
    
            if self.usePhasedContexts:
                mappedLeftSide, mappedRightSide = makeLinkedSides(base, ns, phasedPSide, self.mappedSequenceGraph.sides)
                leftSide.mappedSides = [ mappedLeftSide ]
                rightSide.mappedSides = [ mappedRightSide ]
                phasedPSide = mappedRightSide
        
        if circular: #Need to add an extra adjacency
            firstSide.adjacencies.add((pSide, ''))
            pSide.adjacencies.add((firstSide, ''))
            if self.usePhasedContexts:
                firstSide.mappedSides[0].adjacencies.add((phasedPSide, ''))
                phasedPSide.adjacencies.add((firstSide.mappedSides[0], ''))
    
    def renumber(self, startID=0, prefix=""):
        ids = [ int(side.basePosition.id) for side in self.positiveSides() ]
        print "oh dear", ids
        assert len(ids) == len(set(ids))
        ids.sort()
        for side in self.positiveSides():
            i = int(side.basePosition.id)
            assert i in ids
            i = ids.index(i) + startID
            if prefix != "":
                side.basePosition.id = str(prefix) + "_" + str(i)
            else:
                side.basePosition.id = i
        if self.usePhasedContexts:
            self.mappedSequenceGraph.renumber(startID=startID, prefix=prefix)
    
    def printDotFile(self, graphVizFileHandle, showContextSets, showIDs, number, displayAsSubgraph):
        if displayAsSubgraph:
            graphVizFileHandle.write('subgraph cluster_%s {\nstyle=filled;\ncolor=lightgrey;label = "%s";\n' % (number,self.label))
        else:
            graphVizFileHandle.write('subgraph cluster_%s {\nlabel = "%s";\n' % (number,self.label))
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
        #if displayAsSubgraph:
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
                     help="Number of mismatches to allow between context strings",
                     default=0)
    
    parser.add_option("--dissimilarity", dest="dissimilarity", type="int", 
                     help="Number of differences to require to be defined as unique (e.g. dissimilarity=2 means a context string must be a hamming distance of 2 from all other context strings). Default=1",
                     default=1)
    
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
    
    parser.add_option("--mergeSymmetric", dest="mergeSymmetric", type="string",
                     help="For each sequence graph indexed collapse to given m",
                     default="")
    
    options, args = parser.parse_args()
    
    if len(args) == 0:
        parser.print_help()
        return 1
    
    mergeContigs = [ int(i) for i in options.mergeContigs.split() ]  
    mergeSymmetric = {} 
    for i in options.mergeSymmetric.split():
        mergeSymmetric[int(i.split("=")[0])] = int(i.split("=")[1])
    
    #First create the sequence graphs for each input graph
    sequenceGraphs = []
    
    for index in xrange(0, len(args), 2):
        assembly = args[index]
        label = args[index + 1]
        print "Processing sequence graph", index, options.mismatches, label
        sG = SequenceGraph(options.usePhasedContexts, label=label)
        sequenceGraphs.append(sG)
        for string in assembly.split():
            print "Adding string:", string, " of assembly:", index
            if index in mergeContigs: #Merge the new contig into the previous contigs
                sG2 = SequenceGraph(options.usePhasedContexts)
                sG2.addString(string)
                matches = []
                for side in sG2.positiveSides():
                    leftMatch = sG.getMatch(side, mismatches=options.mismatches, dissimilarity=options.dissimilarity)
                    rightMatch = sG.getMatch(side.otherSide, mismatches=options.mismatches, dissimilarity=options.dissimilarity)
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
            
        if index in mergeSymmetric.keys():
            m = mergeSymmetric[index]
            sG.computeContextSets(minContextLength=options.minContextLength, maxContextLength=options.maxContextLength)
            while True:
                #Try and get node with context string longer than m
                validSides = [ side for side in sG.sides if len([ side2 for side2 in sG.sides if (side2.base() == side.base() and side2 != side and side2.otherSide != side) ]) > 0  ] #Get sides which have a potential merge partner
                if len(validSides) == 0:
                    break
                random.shuffle(validSides)
                side1 = max(validSides, key=lambda side : sG.contextSets[side].maxLength())
                contextSet1 = sG.contextSets[side1]
                if contextSet1.maxLength()-1 < m:
                    break #We are done
                #Find other node with maximum suffix of long context string.
                l = [ side for side in sG.sides if side != side1 and side != side1.otherSide and side.base() == side1.base() ]
                assert len(l) != 0
                #print "frarrr 2", [ (side.basePosition.base, side.basePosition.id, sG.contextSets[side].minimalUniqueStrings) for side in sG.sides ]
                side2 = max(l, key=lambda side : contextSet1.maxSharedPrefixLength(sG.contextSets[side]))
                #Merge sides
                #print "Going to merge", side1.basePosition.id, side2.basePosition.id, contextSet1.maxSharedPrefixLength(sG.contextSets[side2])
                sG.merge(side1, side2)
                #Recompute context sets
                sG.computeContextSets(minContextLength=options.minContextLength, maxContextLength=options.maxContextLength)
     
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
        sG.printDotFile(graphVizFileHandle, index in showContextSets and options.mismatches == 0, index in showIDs, index, len(sequenceGraphs) > 1)
        
    #Now print the matching edges between the graphs
    for index in xrange(1, len(sequenceGraphs)):
        sGInput = sequenceGraphs[index]
        for side in sGInput.positiveSides():
            haveMatched = False
            for pIndex in xrange(index-1, -1, -1):
                sGTarget = sequenceGraphs[pIndex]
                
                leftMatch = sGTarget.getMatch(side, mismatches=options.mismatches, dissimilarity=options.dissimilarity)
                rightMatch = sGTarget.getMatch(side.otherSide, mismatches=options.mismatches, dissimilarity=options.dissimilarity)
                
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
    