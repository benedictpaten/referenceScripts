import sys

fH = open(sys.argv[2], 'w')
pL = {}
l = []
for line in open(sys.argv[1], 'r'):
    if line[:2] == 'a ':
        if l != []:
            pL2 = {}
            l.sort()
            fH.write("%i\n" % len(l))
            for i in l:
                i = i.split()
                b = 0
                b2 = 0
                b3 = 0
                pL2[i[1]] = i
                if pL.has_key(i[1]):
                    b = 1
                    if i[4] == pL[i[1]][4]:
                        b2 = 1
                        if int(i[2]) == int(pL[i[1]][2]) + int(pL[i[1]][3]):
                            b3 = 1
                fH.write("%s\t%s\t%i\t%i\t%i\t%s\n" % (i[1][:10] + " " * (10 - len(i[1][:10])), i[2], b, b2, b3, i[-1][:30]))
            fH.write("\n")
            l = []
            pL = pL2
    elif line[0] == 's':
        l.append(line)
fH.close()
        
                