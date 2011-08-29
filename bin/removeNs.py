
"""Replaces all runs of Ns greater than M in length with M Ns.
"""

import sys, re
from sonLib.bioio import fastaRead, fastaWrite

if len(sys.argv) == 0:
    print "fasta-file-in fasta-file-out minimum-length-of-ns-to-mask"
    sys.exit()
    
class Header:
    def __init__(self, header, lenSeq):
        items = header.split('.')
        if len(items) < 6:
            self.name = header.replace(".", "_")
            self.chr = 'NA'
            self.chrSize = str( lenSeq )
            self.start = 0
            self.fragSize = lenSeq
            self.strand = '1'
        else:
            self.name = items[0]
            self.chr = items[1]
            self.chrSize = items[2]
            self.start = int( items[3] )
            self.fragSize = int( items[4] )
            self.strand = items[5]

    def getStr( self ):
        return '.'.join( [self.name, self.chr, self.chrSize, str(self.start), str(self.fragSize), self.strand] )

def fn( header, sequence, lengthOfNs ):
    pattern = "(?P<Nstr>[Nn]+)"
    searchedSeq = ""
    m = re.search( pattern, sequence )
    while m:
        lenNs = len( m.group('Nstr') )
        if lenNs < lengthOfNs:
            searchedSeq += sequence[: m.start() + lenNs]
        else:
            subSequence = searchedSeq + sequence[: m.start()]
            header.fragSize = len( subSequence )
            newheader = header.getStr()
            yield newheader, subSequence
            searchedSeq = ""
            #Update the start coordinate:
            header.start += len( subSequence ) + lenNs
        
        sequence = sequence[m.start() + lenNs: ]
        m = re.search( pattern, sequence )
    if len(sequence) > 0:
        header.fragSize = len( sequence )
        yield header.getStr(), sequence

#=========== MAIN ====================
fH = open(sys.argv[1], 'r')
fH2 = open(sys.argv[2], 'w')
lengthOfNs = int(sys.argv[3])

for name, sequence in fastaRead(fH):
    header = Header( name, len(sequence) )
    for newheader, subsequence in fn( header, sequence, lengthOfNs ):
        if len( subsequence ) > 0:
            fastaWrite(fH2, newheader, subsequence)
        
fH.close()
fH2.close()
