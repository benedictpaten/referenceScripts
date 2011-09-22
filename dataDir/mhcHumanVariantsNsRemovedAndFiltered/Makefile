sequences= panTro2 cox qbl mann mcf dbb apd ssto hg19 NA12878 NA12891 NA12892 NA19238 NA19239 NA19240

all : ${sequences:%=../mhcHumanVariantsNsRemovedAndFiltered/%}

../mhcHumanVariantsNsRemovedAndFiltered/% : 
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/$* $*  
