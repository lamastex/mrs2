#======== modify the following values ==================

# The main files to be compiled
PROGRAMS = test testLevy testRosen testBiG

OBJ=dim2taylor.o itaylor.o Int.o

# The calc directory
#CALC = /home/tomas
#CALC = /Users/raazesh
CALC=/usr/local
#CALC=/home/gat41/sources

# The (un-)installation prefix
CXSCDIR=$(CALC)/cxsc

#========= the following commands should work on most Unix systems ========

# The Name of the C++ compiler
CXX=g++

# Optional flags to give to the C++ compiler
CXXOPTS=-Wall -Winline#-frounding-math

# Extra flags to give to the C++ compiler
CXSCFLAGS=$(CXSCINC) $(CXXOPTS)

# Additional include path
CXSCINC=-I$(CXSCDIR)/include -L$(CXSCDIR)/lib 

# Names of libraries
CXSCLIBRARIES=-lcxsc -lm

#========== you shouldn't modify anything below ===========================

.SUFFIXES:
.SUFFIXES: .cpp .cc .hpp .o

default:
	@echo
	@echo "C-XSC - C++ library for eXtended Scientific Computation"
	@echo
	@echo "Usage: make <target>, where targets are"
	@echo
	@echo "    XYKM"
	@echo ""
	@echo "    tidy / clean                 (tidy / clean up)"
	@echo "";

# -------------------------------------------------------------------------

all: $(OBJ) $(PROGRAMS)
obj: $(OBJ)

test: test.cc Int.o itaylor.o dim2taylor.o
	@echo "Building 'test'   --- Compute Integrals---"
	$(CXX) -o test $(CXSCFLAGS) test.cc Int.o dim2taylor.o $(CXSCLIBRARIES)

testLevy: testLevy.cc Int.o itaylor.o dim2taylor.o
	@echo "Building 'test'   --- Compute Integrals---"
	$(CXX) -o testLevy $(CXSCFLAGS) testLevy.cc Int.o dim2taylor.o $(CXSCLIBRARIES)

testRosen: testRosen.cc Int.o itaylor.o dim2taylor.o
	@echo "Building 'test'   --- Compute Integrals---"
	$(CXX) -o testRosen $(CXSCFLAGS) testRosen.cc Int.o dim2taylor.o $(CXSCLIBRARIES)

testBiG: testBiG.cc Int.o itaylor.o dim2taylor.o
	@echo "Building 'test'   --- Compute Integrals---"
	$(CXX) -o testBiG $(CXSCFLAGS) testBiG.cc Int.o dim2taylor.o $(CXSCLIBRARIES)

Int.o: Int.cc Int.h dim2taylor.o
	$(CXX) -c Int.cc $(CXSCFLAGS)

itaylor.o: itaylor.hpp itaylor.cpp
	$(CXX) -c itaylor.cpp $(CXSCFLAGS)

dim2taylor.o: dim2taylor.hpp dim2taylor.cpp
	$(CXX) -c dim2taylor.cpp $(CXSCFLAGS)



.PHONY: default all

# -----------------------------------------------------------------------
# TIDY UP
# -----------------------------------------------------------------------

tidy:
	@ - rm -f *.o *~

#-----------------------------------------------------------------------
# CLEAN UP
# -----------------------------------------------------------------------

clean:
	@ - rm -f *.o *~ $(PROGRAMS)

# -----------------------------------------------------------------------

