rootPath=../
cactusRootPath=${rootPath}/../cactus
cactusLibPath=${cactusRootPath}/lib
cactusToolsLibPath=${rootPath}/../cactusTools/lib
assemblaLibPath=${rootPath}/../assemblaLib/lib

include ${cactusRootPath}/include.mk

binPath=${rootPath}/bin

libSources = impl/*.c
libHeaders = inc/*.h
commonSources = impl/referenceCommon.c

extraLibs=${cactusToolsLibPath}/cactusMafs.a ${cactusToolsLibPath}/cactusTreeStats.a ${assemblaLibPath}/assemblaLib.a ${cactusLibPath}/cactusReference.a ${cactusLibPath}/cactusLib.a

programs = pathStats contiguityStats copyNumberStats snpStats coverageStats

all : ${programs:%=${binPath}/%}

${binPath}/%: ${libSources} ${libHeaders} ${basicLibsDependencies} ${extraLibs}
	${cxx} ${cflags} -I ${cactusLibPath} -I ${cactusToolsLibPath} -I ${assemblaLibPath} -I ${libPath} -I inc -o ${binPath}/$* impl/$*.c ${commonSources} ${extraLibs} ${basicLibs} 

clean : ${programs:%=%.clean} 
	rm -rf *.o ${binPath}/*.dSYM

%.clean : 
	rm -rf ${binPath}/$*
