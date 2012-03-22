

all : dataM srcM

srcM :
	cd src && make all

dataM :
	cd dataDir/mhcHumanVariantsNsRemoved && make all
	cd dataDir/mhcHumanVariantsNsRemovedAndFiltered && make all

test : all
	cd tests && make all

run : mhcHumanVariantsNsRemoved

runFiltered : mhcHumanVariantsNsRemovedAndFiltered
	
mhcHumanVariantsNsRemoved : all
	cd main/mhcHumanVariantsNsRemoved && make all

mhcHumanVariantsNsRemovedAndFiltered : all
	cd main/mhcHumanVariantsNsRemovedAndFiltered && make all

clean : 
	cd src && make clean
	cd tests && make clean
	#cd main/mhcHumanVariantsNsRemoved && make clean