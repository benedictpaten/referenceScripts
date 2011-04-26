clean : 
	cd tests/little && make clean
	cd main/mhcHumanVariants && make clean
	cd main/mhcHumanVariantsNsRemoved && make clean

all :
	cd data/mhcHumanVariantsNsRemoved && make all

test : all
	cd tests/little && make all

run : mhcHumanVariants mhcHumanVariantsNsRemoved

mhcHumanVariants : all
	cd main//mhcHumanVariants && make all
	
mhcHumanVariantsNsRemoved : all
	cd main//mhcHumanVariantsNsRemoved && make all