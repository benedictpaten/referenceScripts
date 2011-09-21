minimumNsForScaffoldGap=15
minimumLengthOfFragment=100 
sequences= panTro2 cox qbl mann mcf dbb apd ssto hg19 NA12878.fa NA12891.fa NA12892.fa NA19238.fa NA19239.fa NA19240.fa

all : ${sequences:%=../mhcHumanVariantsNsRemovedAndFiltered/%}

../mhcHumanVariantsNsRemovedAndFiltered/% : 
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/$* $*  
