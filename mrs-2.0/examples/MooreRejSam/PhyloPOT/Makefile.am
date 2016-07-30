
noinst_PROGRAMS = PhyloPOT
PhyloPOT_SOURCES = PhyloPOT.cpp ../../targets/FPhyloPOT.cpp \
../../targets/FPhyloPOT.hpp ../../targets/PhyloTree.cpp \
../../targets/PhyloTree.hpp
EXTRA_DIST = commands 
PhyloPOT_CPPFLAGS = -I${top_builddir}/src -I../../targets
PhyloPOT_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
