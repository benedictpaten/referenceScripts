from sonLib.bioio import fastaRead, fastaWrite
import sys
import random
fH = open(sys.argv[2], "w")
def fn(k, i, j):
    if k != i:
        l = random.choice(j)
        if k == k.upper():
            return l.upper()
        return l.lower()
    else:
        return k
for name, seq in fastaRead(open(sys.argv[1]), "r"):
    for i, j in [ ("W", ("A", "T")),
                 ("S", ("C", "G")), 
                 ("M", ("A", "C")),
                 ("K", ("G", "T")),
                 ("R", ("A", "G")),
                 ("Y", ("C", "T")),
                 ("B", ("C", "G", "T")),
                 ("D", ("A", "G", "T")),
                 ("H", ("A", "C", "T")),
                 ("V", ("A", "C", "G")) ]:
        seq = "".join([ fn(k, i, j) for k in seq ])
    fastaWrite(fH, name, seq)
fH.close()