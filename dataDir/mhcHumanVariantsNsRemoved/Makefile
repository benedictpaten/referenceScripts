minimumNsForScaffoldGap=10
minimumLengthOfFragment=0
grcSequences = apd dbb mann qbl cox hg19 mcf panTro3 ssto
sequences= ${grcSequences} NA12892 NA19239 nigerian venter yanhuang NA19238 NA19240 

all : ${sequences:%=../mhcHumanVariantsNsRemoved/%} NA12878 constraints.cig

clean : 
	rm -rf ${sequences} constraints.cig NA12878

constraints.cig : ../mhcHumanVariants/allUniq.cigar ${grcSequences}
	python ../../src/scripts/cigarConvertRef.py ../mhcHumanVariants/allUniq.cigar constraints.cig ${grcSequences}

../mhcHumanVariantsNsRemoved/% : 
	python ../../src/scripts/removeNs.py ../mhcHumanVariants/$* $* ${minimumNsForScaffoldGap} ${minimumLengthOfFragment} DEBUG

NA12878 :
	python ../../src/scripts/makeHaploid.py ../mhcHumanVariants/NA12878 NA12878_haploid
	python ../../src/scripts/removeNs.py NA12878_haploid NA12878 ${minimumNsForScaffoldGap} ${minimumLengthOfFragment} DEBUG