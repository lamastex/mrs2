noinst_PROGRAMS = TestPCF

TestPCF_SOURCES = 	TestPCF.cpp TestPCF.hpp \
					TestPCFTools.cpp TestPCFTools.hpp \
					TestPCFBasic.cpp \
					TestPCFFunEst.cpp \
					TestPCFPtWiseCoverage.cpp \
					TestPCFNormalise.cpp \
					TestPCFMarginalise.cpp \
					TestPCFSlice.cpp \
					TestPCFL1.cpp \
					TestPCFArithmetic.cpp \
					TestPCFMakeFromHist.cpp \
					TestPCFIAE.cpp \
					TestPCFLogLik.cpp \
					TestPCFRandom.cpp \
					TestPCFSmear.cpp \
					TestPCFArithmeticSpecial.cpp \
					simpleFobj1.cpp simpleFobj1.hpp \
					simpleFobj2.cpp simpleFobj2.hpp 

						
#EXTRA_DIST =

AM_CPPFLAGS = -I${top_builddir}/src #-I$(TARGETS)
LDADD = -lcxsc ${top_builddir}/src/libmrs.a

