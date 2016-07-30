noinst_PROGRAMS = TestSOR

TARGETS = ../../MappedTargets

TestSOR_SOURCES = 	TestSOR.cpp TestSOR.hpp \
					TestSORTools.cpp TestSORTools.hpp \
					TestSORBVBasic.cpp \
					TestSORBVArithmetic.cpp \
					TestSORBVSlice.cpp \
					TestSORBVFunImageEst.cpp \
					TestSORArithmeticExtended.cpp \
					TestSORSliceExtended.cpp \
					$(TARGETS)/SphereFobj.cpp $(TARGETS)/SphereFobj.hpp

						
#EXTRA_DIST =

AM_CPPFLAGS = -I${top_builddir}/src -I$(TARGETS)
LDADD = -lcxsc ${top_builddir}/src/libmrs.a

