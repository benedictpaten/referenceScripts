samples = {}
samplesList = samples.keys()
samplesList.sort()

for i in xrange(len(samplesList)):
    sample1 = samplesList[i]
    for j in xrange(i+1,len(samplesList)):
        sample2 = samplesList[j]
        commonSnps = samples[sample1].intersection(samples[sample2])
        allSnps = samples[sample1].union(samples[sample2])
        print "\t".join([ sample1, sample2, str(len(commonSnps)), str(len(allSnps)) ])
