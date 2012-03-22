grcSequences = apd dbb mann qbl cox mcf panTro3 ssto
sequences=${grcSequences} NA12878 NA12892 NA19239 nigerian venter yanhuang  NA19238 NA19240

all : filterSeqs hg19 constraints.cig

clean : 
	rm -rf ${sequences} hg19 hg19_long constraints.cig

filterSeqs : ${sequences:%=../mhcHumanVariantsNsRemovedAndFiltered/%}

constraints.cig : ../mhcHumanVariants/allUniq.cigar ${grcSequences} hg19
	python ../../src/scripts/cigarConvertRef.py ../mhcHumanVariants/allUniq.cigar constraints.cig ${grcSequences} hg19

../mhcHumanVariantsNsRemovedAndFiltered/% : 
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/$* $*  

hg19 :
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/hg19 hg19_long  
	python ../../src/scripts/cropSequence.py filterNonComponentSequences.xml hg19_long hg19