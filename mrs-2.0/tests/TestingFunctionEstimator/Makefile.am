noinst_PROGRAMS = TestFunEst TestFunEstPaperExample

TARGETS = ../../examples/MappedTargets

TestFunEst_SOURCES = 	TestFunEst.cpp \
						simpleFobj1.cpp simpleFobj1.hpp \
						simpleFobj2.cpp simpleFobj2.hpp \
						$(TARGETS)/oscFobj1.cpp $(TARGETS)/oscFobj1.hpp

TestFunEstPaperExample_SOURCES = 	TestFunEst_PaperExample.cpp \
						$(TARGETS)/oscFobj1.cpp $(TARGETS)/oscFobj1.hpp

						
#EXTRA_DIST =

AM_CPPFLAGS = -I${top_builddir}/src -I$(TARGETS)
LDADD = -lcxsc ${top_builddir}/src/libmrs.a

