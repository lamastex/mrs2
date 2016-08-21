noinst_PROGRAMS = CVOptMAP
TARGETS = ../../MixtureTargets ../PointClouds
CVOptMAP_SOURCES = CVOptMAP.cpp ../PointClouds/dataprep.hpp ../PointClouds/testDenCommon.cpp ../PointClouds/testDenCommon.hpp ../PointClouds/testDenTools.cpp ../PointClouds/testDenTools.hpp ../../MixtureTargets/mixture_mvn.cpp ../../MixtureTargets/mixture_mvn.hpp
CVOptMAP_CPPFLAGS = -I${top_builddir}/src -I./  -I../../MixtureTargets -I../PointClouds
CVOptMAP_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
