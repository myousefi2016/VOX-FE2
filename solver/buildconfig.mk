SHELL           = /bin/sh
CXX             = mpiCC
LD              = $(CXX)
EXE             = voxfe_solver

MAKE_OPTIMIZED  = 1

ifeq ($(MAKE_OPTIMIZED),1)
   DEBUG_FLAGS =
   OPT_FLAGS   = -O3 -g
else
    DEBUG_FLAGS = -g -O0
    OPT_FLAGS   =
endif

GLOBAL_CXXFLAGS	= $(OPT_FLAGS) $(DEBUG_FLAGS)
LDFLAGS		= $(OPT_FLAGS)
CXXLIBS		=

OBJECTEXT       = o

%.$(OBJECTEXT): %.cpp
	$(CXX) $(CXXFLAGS) -c $<
