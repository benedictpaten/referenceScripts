import sys
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

#Run command: python2.7-32 ~/sync/eclipse/git/referenceScripts/src/scripts/contigPathPlot.py results1/pathStats_reference.xml 1 results2/pathStats_reference.xml 2 results3/pathStats_reference.xml 3 results4/pathStats_reference.xml 4 results5/pathStats_reference.xml 5 results6/pathStats_reference.xml 6 results7/pathStats_reference.xml 7 results8/pathStats_reference.xml 8 results1/pathStats_hg19.xml hg19

results = [ (j, [ (k.attrib["sampleName"], k.attrib["contigPathN50"]) \
                 for k in ET.parse(i).getroot().findall("statsForSample") \
                 if k.attrib["sampleName"] not in ("", "reference", "hg19", "ROOT") ]) \
           for i, j in zip(sys.argv[1::2], sys.argv[2::2]) ]

for i, j in results:
    j.sort()
    
print "samples\t", '\t'.join([ i for i, j in results[0][1] ] + [ "aggregate" ])
for i, j in results:
    j.append(("aggregate", sum([ float(l) for k, l in j ])/len(j)))
    print i,"\t", "\t".join([ str(l) for k, l in j ])

r = [ [ float(k[1]) for k in j ] for i, j in results ] 

i = range(len(r[0]))
plt.plot(i, r[0], "b_", 
         i, r[1], 'g_', 
         i, r[2], 'r_', 
         i, r[3], 'c_', 
         i, r[4], 'mo', 
         i, r[5], 'y_', 
         i, r[6], 'k_', 
         i, r[7], 'w_',
         i, r[8], 'wo')
plt.show()