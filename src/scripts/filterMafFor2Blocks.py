import sys
l = []
fH = open(sys.argv[2], 'w')
for line in open(sys.argv[1], 'r').readlines():
    if line[:2] == 'a ':
        if len(l) == 4:
            for line in l:
                fH.write(l)
            fH.write("\n")
        l = [ line ]
    elif line[:2] == 's ':
        l.append(l)
fH.close()