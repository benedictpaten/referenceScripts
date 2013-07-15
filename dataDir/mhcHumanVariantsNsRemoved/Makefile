minimumNsForScaffoldGap=10
minimumLengthOfFragment=0
grcSequences = apd dbb mann qbl cox hg19 mcf panTro3 ssto
sequences= ${grcSequences} NA12892 NA19239 NA19238 NA19240 NA12878 nigerian venter yanhuang cgIndelContigs.fa 1kIndelContigs.fa snp135IndelContigs.fa
baseDir=../mhcHumanVariants
constraints= ${baseDir}/all.cigar ${baseDir}/1kIndelConstraints.cig ${baseDir}/cgIndelConstraints.cig ${baseDir}/snp135IndelConstraints.cig

all : ${sequences:%=../mhcHumanVariantsNsRemoved/%} 
#constraints.cig

../mhcHumanVariantsNsRemoved/% : 
	python ../../src/scripts/makeHaploid.py ${baseDir}/$* $*_haploid
	python ../../src/scripts/removeNs.py $*_haploid $* ${minimumNsForScaffoldGap} ${minimumLengthOfFragment} DEBUG
	rm $*_haploid

constraints.cig : ${constraints} ${sequences}
	cat ${constraints} > constraints.cig.pre
	python ../../src/scripts/cigarConvertRef.py constraints.cig.pre constraints.cig ${sequences}
	rm constraints.cig.pre

clean : 
	rm -rf ${sequences} constraints.cig