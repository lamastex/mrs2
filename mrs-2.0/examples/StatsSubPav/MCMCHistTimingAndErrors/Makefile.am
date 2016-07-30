noinst_PROGRAMS = testHistDen 

# MixtureTargets is two directories down here
TARGETS = ../../MixtureTargets

				
testHistDen_SOURCES = testHistDen.cpp \
				testDenCommon.cpp testDenCommon.hpp \
				testDenTools.cpp testDenTools.hpp \
				../../MixtureTargets/mixture_mvn.cpp ../../MixtureTargets/mixture_mvn.hpp


AM_CPPFLAGS = -I${top_builddir}/src
#AM_CPPFLAGS = -I${top_builddir}/src -I$(TARGETS) -I$(CARVER)
LDADD = -lcxsc ${top_builddir}/src/libmrs.a

testHistDen_CPPFLAGS = $(AM_CPPFLAGS) -I../../MixtureTargets

