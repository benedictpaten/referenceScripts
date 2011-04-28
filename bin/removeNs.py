
"""Replaces all runs of Ns greater than M in length with M Ns.
"""

import sys
from sonLib.bioio import fastaRead, fastaWrite

if len(sys.argv) == 0:
    print "fasta-file-in fasta-file-out minimum-length-of-ns-to-mask"
    sys.exit()
    
def fn(sequence, start, pattern):
    while pattern in sequence:
        i = sequence.index(pattern)
        if i > 0:
            yield sequence[:i], start
        start += i+len(pattern)
        sequence = sequence[i+len(pattern):]
    if len(sequence) > 0:
        yield sequence, start

fH = open(sys.argv[1], 'r')
fH2 = open(sys.argv[2], 'w')
lengthOfNs = int(sys.argv[3])
for name, sequence in fastaRead(fH):
    sequence = ("N"*lengthOfNs).join([ subsequence for subsequence, start in fn(sequence, 0, "N"*lengthOfNs) ])
    if len(sequence) > 0:
        fastaWrite(fH2, name, sequence)
        
fH.close()
fH2.close()