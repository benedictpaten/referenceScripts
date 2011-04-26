
"""Removes runs of N's from longer files, splitting them into seperate contigs.

#Ngan's coordinates system
>apd.chr6_apd_hap1.4622290.0.4622290.1
contig-name / contig-length / offset (zero based, from 5 end) / length / strand (0 = negative / 1 = positive)

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
    strand = int(name.split(".")[-1])
    length = int(name.split(".")[-2])
    start = int(name.split(".")[-3])
    rest = ".".join((name.split(".")[:-3]))
    for subSequence, start in fn(sequence, start, "N"*lengthOfNs):
        while len(subSequence) > 0 and subSequence[0] == 'N':
            subSequence = subSequence[1:]
            start += 1 #Increment the start coordinate by 1.
        if len(subSequence) > 0:
            fastaWrite(fH2, ".".join((rest, str(start), str(len(subSequence)), str(strand))), subSequence)
        
fH.close()
fH2.close()