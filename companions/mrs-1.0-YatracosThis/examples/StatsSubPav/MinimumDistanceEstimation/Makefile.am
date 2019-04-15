noinst_PROGRAMS = MappedGaussianMDE 

TARGETS = ../../MappedTargetsTrunk

MappedGaussianMDE_SOURCES = MappedGaussianMDE.cpp \
						../../MappedTargetsTrunk/GaussianFobj.cpp \
						../../MappedTargetsTrunk/GaussianFobj.hpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp					

EXTRA_DIST = 
MappedGaussianMDE_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk
MappedGaussianMDE_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
