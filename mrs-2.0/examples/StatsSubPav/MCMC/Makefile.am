noinst_PROGRAMS = MCMCuniform MCMCGaussian MCMCnewtest \
					MCMCFunctionSimGaussian

TARGETS = ../../MappedTargets

MCMCuniform_SOURCES = MCMCuniform.cpp dataprep.hpp

MCMCGaussian_SOURCES = MCMCGaussian.cpp dataprep.hpp

MCMCnewtest_SOURCES = MCMCnewtest.cpp dataprep.hpp

MCMCFunctionSimGaussian_SOURCES = MCMCFunctionSimGaussian.cpp \
						../../MappedTargets/GaussianFobj.cpp ../../MappedTargets/GaussianFobj.hpp

EXTRA_DIST =

AM_CPPFLAGS = -I${top_builddir}/src
LDADD = -lcxsc ${top_builddir}/src/libmrs.a

MCMCFunctionSimGaussian_CPPFLAGS = -I${top_builddir}/src -I../../MappedTargets

