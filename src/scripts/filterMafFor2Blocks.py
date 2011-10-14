import sys
l = []
fH = open(sys.argv[2], 'w')
for line in open(sys.argv[1], 'r').readlines():
    if line[:1] == 'a':
        if len(l) == 4 and "hg19" in " ".join(l):
            print "making block"
            for i in l:
                fH.write(i)
            fH.write("\n")
        else:
            print "invalid block", len(l)
        l = [ line ]
    elif line[:1] == 's':
        l.append(line)
fH.close()