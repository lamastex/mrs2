noinst_PROGRAMS = MappedGaussianMDE MappedRosenbrockMDE MDETxt

TARGETS = ../../MappedTargetsTrunk

MappedGaussianMDE_SOURCES = MappedGaussianMDE.cpp \
						../../MappedTargetsTrunk/GaussianFobj.cpp \
						../../MappedTargetsTrunk/GaussianFobj.hpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp					

MappedRosenbrockMDE_SOURCES = MappedRosenbrockMDE.cpp \
						../../MappedTargetsTrunk/RosenDensityFobj.cpp \
						../../MappedTargetsTrunk/RosenDensityFobj.hpp \
						../../targets/FRosenbrock.cpp \
						../../targets/FRosenbrock.hpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp	

MDETxt_SOURCES = MDETxt.cpp \
						../../MappedTargetsTrunk/GaussianFobj.cpp \
						../../MappedTargetsTrunk/GaussianFobj.hpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp	
EXTRA_DIST = 
MappedGaussianMDE_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk
MappedGaussianMDE_LDADD = -lcxsc ${top_builddir}/src/libmrs.a

MappedRosenbrockMDE_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk -I../../targets
MappedRosenbrockMDE_LDADD = -lcxsc ${top_builddir}/src/libmrs.a

MDETxt_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk
MDETxt_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
