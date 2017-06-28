SHELL = /bin/sh
.SUFFIXES: .cpp .c .o .cxx
.PHONY: clean cleanall all

CC=/usr/bin/g++

#OPTS =  -O2 -funroll-loops
#OPTS =  -O0 -g -DDEBUG 
#PROFILING  = -pg 

TARGETS = convertVoxFEScript
READER_SRC = ReadVoxFEScript.cxx ReadVoxFEScript.hxx  ReadVoxFEScript.h
READER_OBJ = ReadVoxFEScript.o

C_SOURCES = convertVoxFEScriptMain.cpp convertVoxFEScript.hxx  convertVoxFEScript.h
C_OBJ = convertVoxFEScriptMain.o 

LIBPATH = 
CFLAGS  = $(PROFILING)
HEADERS = -I.
LIBS    = -lm

all: convertVoxFEScript

convertVoxFEScript:  ${C_OBJ} ${READER_OBJ}
	${CC}  $(CFLAGS) $(OPTS)  $(HEADERS) -o $@ ${C_OBJ} ${READER_OBJ} $(LIBS)


.cpp.o:
	${CC} $(CFLAGS) $(HEADERS) $(OPTS)  -c $<

.cxx.o:
	${CC} $(CFLAGS) $(HEADERS) $(OPTS)  -c $<	

.h:


clean:
	rm -f core* *.o


cleanall:
	rm -f core* *.o ${TARGETS}

convertVoxFEScript.o: ${C_SOURCES}
ReadVoxFEScript.o: ${READER_SRC}
GlobalAssembly.o: ${G_SOURCES}

convertVoxFEScript.hxx:
convertVoxFEScript.h:
ReadVoxFEScript.hxx:
ReadVoxFEScript.h:
GlobalAssembly.cpp:
GlobalAssembly.h:

