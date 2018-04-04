noinst_PROGRAMS = MappedSPSamplingGaussian MappedSPSamplingTest
MappedSPSamplingGaussian_SOURCES = MappedSPSamplingGaussian.cpp \
GaussianFobj1D.cpp GaussianFobj1D.hpp GaussianFobj2D.cpp GaussianFobj2D.hpp \
GaussianFobj9D.cpp GaussianFobj9D.hpp \
GaussianFobj10D.cpp GaussianFobj10D.hpp \
GaussianFobj100D.cpp GaussianFobj100D.hpp \
RosenFobj2D.cpp RosenFobj2D.hpp \
RosenFobj10D.cpp RosenFobj10D.hpp \
../MCMCGRSims/auto_tools.hpp ../MCMCGRSims/auto_tools.cpp \
../BivGaussian/dataprep.hpp ../MCMCGR/MCMCGRtools.cpp ../MCMCGR/MCMCGRtools.hpp
MappedSPSamplingTest_SOURCES = MappedSPSamplingTest.cpp \
GaussianFobj1D.cpp GaussianFobj1D.hpp GaussianFobj2D.cpp GaussianFobj2D.hpp \
GaussianFobj10D.cpp GaussianFobj10D.hpp \
GaussianFobj9D.cpp GaussianFobj9D.hpp \
GaussianFobj100D.cpp GaussianFobj100D.hpp \
RosenFobj2D.cpp RosenFobj2D.hpp \
RosenFobj10D.cpp RosenFobj10D.hpp \
../MCMCGRSims/auto_tools.hpp ../MCMCGRSims/auto_tools.cpp \
../BivGaussian/dataprep.hpp ../MCMCGR/MCMCGRtools.cpp ../MCMCGR/MCMCGRtools.hpp



#EXTRA_DIST =
MappedSPSamplingGaussian_CPPFLAGS = -I${top_builddir}/srcGloria -I../BivGaussian -I../MCMCGR \
-I${top_builddir}/srcGloria -I../../targets -I../BivGaussian -I../MCMCGRSims
MappedSPSamplingGaussian_LDADD = -lcxsc ${top_builddir}/srcGloria/libmrs.a
MappedSPSamplingTest_CPPFLAGS = -I${top_builddir}/srcGloria -I../BivGaussian -I../MCMCGR \
-I${top_builddir}/srcGloria -I../../targets -I../BivGaussian -I../MCMCGRSims
MappedSPSamplingTest_LDADD = -lcxsc ${top_builddir}/srcGloria/libmrs.a
