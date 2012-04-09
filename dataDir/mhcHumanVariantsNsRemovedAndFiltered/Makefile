grcSequences = apd dbb mann qbl cox mcf panTro3 ssto
sequences= ${grcSequences} NA12892 NA19239 NA19238 NA19240 NA12878 nigerian venter yanhuang cgIndelContigs.fa 1kIndelContigs.fa snp135IndelContigs.fa
baseDir=../mhcHumanVariants
constraints= ${baseDir}/all.cigar ${baseDir}/1kIndelConstraints.cig ${baseDir}/cgIndelConstraints.cig ${baseDir}/snp135IndelConstraints.cig

all : filterSeqs hg19 constraints.cig

clean : 
	rm -rf ${sequences} hg19 hg19_long constraints.cig

filterSeqs : ${sequences:%=../mhcHumanVariantsNsRemovedAndFiltered/%}

../mhcHumanVariantsNsRemovedAndFiltered/% : 
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/$* $*  

constraints.cig : ${constraints} ${sequences} hg19
	cat ${constraints} > constraints.cig.pre
	python ../../src/scripts/cigarConvertRef.py constraints.cig.pre constraints.cig ${sequences} hg19
	rm constraints.cig.pre

hg19 :
	python ../../src/scripts/filterNonComponentSequences.py filterNonComponentSequences.xml ../mhcHumanVariantsNsRemoved/hg19 hg19_long  
	python ../../src/scripts/cropSequence.py filterNonComponentSequences.xml hg19_long hg19
	rm hg19_long