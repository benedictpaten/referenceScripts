
import sys
import os.path
from optparse import OptionParser

parser = OptionParser()
 
parser.add_option("--branchLength", dest="branchLength", default=0.002)
parser.add_option("--outGroupTree", dest="outGroupTree", default="x;")
parser.add_option("--filePath", dest="filePath", default="")

options, args = parser.parse_args()

def fn(file):
    i = os.path.split(sequence)[-1]
    if ".fa" == i[-3:]:
        return i[:-3]
    return i

print "NEWICK_TREE:" + options.outGroupTree.replace("x", "(" + ",".join([ "%s:%f" % (fn(sequence), float(options.branchLength)) for sequence in args ]) + ")")
print "REQUIRED_SPECIES_ARGS: " + " ".join([ "%s" % fn(sequence) for sequence  in args ]) 
print "FILES: " + " ".join([ os.path.join(options.filePath, i) for i in args ])