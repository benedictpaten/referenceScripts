minimumNsForScaffoldGap=15
minimumLengthOfFragment=100 
sequences= panTro2 cox qbl mann mcf dbb apd ssto hg19 NA12878.fa NA12891.fa NA12892.fa NA19238.fa NA19239.fa NA19240.fa

all : ${sequences:%=../mhcHumanVariantsNsRemovedAndFiltered/%}

../mhcHumanVariantsNsRemovedAndFiltered/% : 
	python ../../src/scripts/filterNonComponentSequences.py ../../output/main/mhcHumanVariantsNsRemoved/required-species-single-copy-species_2-maxWeight-2-0-True-10000-100-1e-05-False/filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/$* $*  
