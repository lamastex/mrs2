noinst_PROGRAMS = MappedGaussianMDE MappedRosenbrockMDE UniformMDE MDETest CheckBounds

TARGETS = ../../MappedTargetsTrunk

MappedGaussianMDE_SOURCES = MappedGaussianMDE.cpp \
						../../MappedTargetsTrunk/GaussianFobj.cpp \
						../../MappedTargetsTrunk/GaussianFobj.hpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp	\
						mdeTools.cpp mdeTools.hpp				

MappedRosenbrockMDE_SOURCES = MappedRosenbrockMDE.cpp \
						../../MappedTargetsTrunk/RosenDensityFobj.cpp \
						../../MappedTargetsTrunk/RosenDensityFobj.hpp \
						../../targets/FRosenbrock.cpp \
						../../targets/FRosenbrock.hpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp	 \
						mdeTools.cpp mdeTools.hpp				

UniformMDE_SOURCES = UniformMDE.cpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp	\
						mdeTools.cpp mdeTools.hpp	\
						../../MappedTargetsTrunk/UniformFobj.cpp \
						../../MappedTargetsTrunk/UniformFobj.hpp

MDETest_SOURCES = MDETest.cpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp \
                                                mdeTools.cpp mdeTools.hpp			
						
							
CheckBounds_SOURCES = CheckBounds.cpp \
						testDenCommon.cpp testDenCommon.hpp \
						testDenTools.cpp testDenTools.hpp	\
						mdeTools.cpp mdeTools.hpp	\
						../../MappedTargetsTrunk/UniformFobj.cpp \
						../../MappedTargetsTrunk/UniformFobj.hpp \
						../../MappedTargetsTrunk/GaussianFobj.cpp \
						../../MappedTargetsTrunk/GaussianFobj.hpp \
						../../MappedTargetsTrunk/RosenDensityFobj.cpp \
						../../MappedTargetsTrunk/RosenDensityFobj.hpp \
						../../targets/FRosenbrock.cpp \
						../../targets/FRosenbrock.hpp

EXTRA_DIST = 
MappedGaussianMDE_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk
MappedGaussianMDE_LDADD = -lcxsc ${top_builddir}/src/libmrs.a

MappedRosenbrockMDE_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk -I../../targets
MappedRosenbrockMDE_LDADD = -lcxsc ${top_builddir}/src/libmrs.a

UniformMDE_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk -I../../targets
UniformMDE_LDADD = -lcxsc ${top_builddir}/src/libmrs.a

MDETest_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk
MDETest_LDADD = -lcxsc ${top_builddir}/src/libmrs.a

CheckBounds_CPPFLAGS = -I${top_builddir}/src \
-I../../MappedTargetsTrunk
CheckBounds_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
