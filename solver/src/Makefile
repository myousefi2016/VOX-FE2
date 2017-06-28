# subdir Makefile

TOP=..
BIN=$(TOP)/bin

include $(TOP)/buildconfig.mk

HDRS = Commands.h Common.h VoxFE.h pFESolver.h  Node.h Elem.h Cons.h Material.h GradientMtx.h IntegrationMtx.h PropertyMtx.h LocalStiffnessMtx.h

CXXSRCS = main.cpp VoxFE.cpp pFESolver.cpp  Node.cpp Elem.cpp Cons.cpp Material.cpp GradientMtx.cpp IntegrationMtx.cpp PropertyMtx.cpp LocalStiffnessMtx.cpp

OBJS  = $(CXXSRCS:.cpp=.$(OBJECTEXT))

CXXFLAGS = $(GLOBAL_CXXFLAGS) $(CXXLIBS)

$(EXE): $(OBJS)
	$(LD) $(LDFLAGS) $^ -o $@

	mv $(EXE) $(BIN)

.PHONY : clean
clean :
	-$(RM) $(OBJS) $(EXE)
	-$(RM) $(BIN)/$(EXE)
