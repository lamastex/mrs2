noinst_PROGRAMS = 	TestGaussianFobj \
					TestUniformFobj \
					TestRosenDensityFobj \
					TestLevyDensity2DFobj
					
TestGaussianFobj_SOURCES = TestGaussianFobj.cpp \
							GaussianFobj.cpp GaussianFobj.hpp
TestUniformFobj_SOURCES = TestUniformFobj.cpp \
							UniformFobj.cpp UniformFobj.hpp

TestRosenDensityFobj_SOURCES = TestRosenDensityFobj.cpp \
							RosenDensityFobj.cpp RosenDensityFobj.hpp
							
TestLevyDensity2DFobj_SOURCES = TestLevyDensity2DFobj.cpp \
							LevyDensityFobj2D.cpp LevyDensityFobj2D.hpp

#EXTRA_DIST =

AM_CPPFLAGS = -I${top_builddir}/src
LDADD = -lcxsc ${top_builddir}/src/libmrs.a
