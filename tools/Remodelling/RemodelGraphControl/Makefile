SHELL = /bin/sh
.SUFFIXES: .o .cxx .h .hxx
.PHONY: clean cleanall all remake

PLUGININC = ../../../gui/VoxFETools
VPATH=${PLUGININC}/itkReader:${PLUGININC}/strainReader

CXX=g++

CFLAGS   = -g
OPTS     = 
SRC      = main.cxx voxfeRemodelControl.cxx ReadVoxFEScript.cxx voxfeStrainData.cxx
HEADERS  = voxfeRemodelControl.h ReadVoxFEScript.hxx ReadVoxFEScript.h voxfeStrainData.h
OBJ      = main.o voxfeRemodelControl.o ReadVoxFEScript.o voxfeStrainData.o voxfeRecomputeStrainFilter.o
INCLUDES = -I.  -I${PLUGININC} -I${PLUGININC}/strainReader
LIBS     = -lm


all: voxfeRemodelControl

voxfeRemodelControl: ${OBJ} 
	${CXX} $(CFLAGS) ${INCLUDES} -o $@ ${OBJ} ${LIBS}

.cxx.o:
	${CXX} $(CFLAGS) ${INCLUDES} $(OPTS)  -c $<	

.h:

.hxx:

voxfeRemodelControl.cxx: voxfeRemodelControl.h
ReadVoxFEScript.cxx: ReadVoxFEScript.hxx ReadVoxFEScript.h
voxfeStrainData.cxx: voxfeStrainData.h
voxfeRecomputeStrainFilter.cxx: voxfeRecomputeStrainFilter.h

#just clean the stuff from the current source dir
clean:
	rm -f main.o voxfeRemodelControl.o voxfeRemodelControl

cleanall:
	rm -f core* *.o ${TARGETS}
	


