minimumNsForScaffoldGap=10
minimumLengthOfFragment=0
sequences=NA12892 NA19239 apd dbb mann nigerian qbl venter yanhuang NA19238 NA19240 cox hg19 mcf panTro3 ssto

all : ${sequences:%=../mhcHumanVariantsNsRemoved/%} NA12878

../mhcHumanVariantsNsRemoved/% : 
	python ../../src/scripts/removeNs.py ../mhcHumanVariants/$* $* ${minimumNsForScaffoldGap} ${minimumLengthOfFragment} DEBUG

NA12878 :
	python ../../src/scripts/makeHaploid.py ../mhcHumanVariants/NA12878 NA12878_haploid
	python ../../src/scripts/removeNs.py NA12878_haploid NA12878 ${minimumNsForScaffoldGap} ${minimumLengthOfFragment} DEBUG