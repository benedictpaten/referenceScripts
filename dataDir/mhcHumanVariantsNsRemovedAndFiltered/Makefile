sequences=NA12878 NA12892 NA19239 apd dbb mann nigerian qbl venter yanhuang  NA19238 NA19240 cox mcf panTro3 ssto

all : filterSeqs hg19

filterSeqs : ${sequences:%=../mhcHumanVariantsNsRemovedAndFiltered/%}

../mhcHumanVariantsNsRemovedAndFiltered/% : 
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/$* $*  

hg19 :
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/hg19 hg19_long  
	python ../../src/scripts/cropSequence.py filterNonComponentSequences.xml hg19_long hg19