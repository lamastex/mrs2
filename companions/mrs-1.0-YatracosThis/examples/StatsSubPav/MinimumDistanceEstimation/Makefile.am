noinst_PROGRAMS = MappedScheffe MappedScheffeIAE MappedScheffeKL

TARGETS = ../../MappedTargetsTrunk

MappedScheffe_SOURCES = MappedScheffe.cpp \
						../../MappedTargetsTrunk/GaussianFobj.cpp \
						../../MappedTargetsTrunk/GaussianFobj.hpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp					
MappedScheffeIAE_SOURCES = MappedScheffeIAE.cpp \
						../../MappedTargetsTrunk/GaussianFobj.cpp \
						../../MappedTargetsTrunk/GaussianFobj.hpp
MappedScheffeKL_SOURCES = MappedScheffeKL.cpp \
						../../MappedTargetsTrunk/GaussianFobj.cpp \
						../../MappedTargetsTrunk/GaussianFobj.hpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp													

EXTRA_DIST = 
MappedScheffe_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk
MappedScheffe_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
MappedScheffeIAE_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk
MappedScheffeIAE_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
MappedScheffeKL_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk -I../../targets
MappedScheffeKL_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
