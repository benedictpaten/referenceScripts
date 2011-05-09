

all : dataM srcM

srcM :
	cd src && make all

dataM :
	cd data/mhcHumanVariantsNsRemoved && make all

test : all
	cd tests/little && make all

run : mhcHumanVariantsNsRemoved
	
mhcHumanVariantsNsRemoved : all
	cd main/mhcHumanVariantsNsRemoved && make all

clean : 
	cd src && make clean
	cd tests/little && make clean
	cd tests/big && make clean
	cd main/mhcHumanVariantsNsRemoved && make clean