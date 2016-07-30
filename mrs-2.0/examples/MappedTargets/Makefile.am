noinst_PROGRAMS = 	TestGaussianFobj \
					TestUniformFobj \
					TestOscFobj \
					TestRosenDensityFobj \
					TestRosenFobj \
					TestLevyDensity2DFobj \
					TestSphereFobj
					
TestGaussianFobj_SOURCES = TestGaussianFobj.cpp \
							GaussianFobj.cpp GaussianFobj.hpp

TestUniformFobj_SOURCES = TestUniformFobj.cpp \
							UniformFobj.cpp UniformFobj.hpp

TestOscFobj_SOURCES = TestOscFobj.cpp \
							oscFobj1.cpp oscFobj1.hpp

TestRosenDensityFobj_SOURCES = TestRosenDensityFobj.cpp \
							RosenDensityFobj.cpp RosenDensityFobj.hpp

TestRosenFobj_SOURCES = TestRosenFobj.cpp \
							RosenFobj.cpp RosenFobj.hpp
							
TestLevyDensity2DFobj_SOURCES = TestLevyDensity2DFobj.cpp \
							LevyDensityFobj2D.cpp LevyDensityFobj2D.hpp

TestSphereFobj_SOURCES = TestSphereFobj.cpp \
							SphereFobj.cpp SphereFobj.hpp


#EXTRA_DIST =

AM_CPPFLAGS = -I${top_builddir}/src
LDADD = -lcxsc ${top_builddir}/src/libmrs.a
