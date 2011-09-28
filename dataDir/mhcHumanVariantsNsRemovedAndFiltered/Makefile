sequences=NA12878 NA12892 NA19239 apd dbb mann nigerian qbl venter yanhuang NA12891 NA19238 NA19240 cox hg19 mcf panTro3 ssto watson
#panTro2 cox qbl mann mcf dbb apd ssto hg19 NA12878 NA12891 NA12892 NA19238 NA19239 NA19240

all : ${sequences:%=../mhcHumanVariantsNsRemovedAndFiltered/%}

../mhcHumanVariantsNsRemovedAndFiltered/% : 
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/$* $*  
