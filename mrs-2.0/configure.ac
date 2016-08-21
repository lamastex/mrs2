
AC_INIT([mrs], [2.0], [raazesh.sainudiin@gmail.com])

#AC_PREREQ(2.61)

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/cxsc.hpp])

# added to make it shut-up about `subdir-objects`
AM_INIT_AUTOMAKE([subdir-objects]) 

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB

#AC_LANG_CPLUSPLUS

# Checks for libraries.
#AC_CHECK_LIB([cxsc], [main], [],
#    AC_MSG_WARN([require libcxsc for building programs]))

AC_CHECK_LIB([m], [sqrt], [],
    AC_MSG_ERROR([cannot find -m]))

AC_CHECK_LIB([gslcblas], [main], [],
    AC_MSG_ERROR([cannot find -lgslcblas]))

AC_CHECK_LIB([gsl], [gsl_version], [],
    AC_MSG_ERROR([cannot find -lgsl]))

AC_CHECK_LIB([cxsc], [main], [],
    AC_MSG_ERROR([cannot find -lcxsc]))

# FIXME: Not sure what functions this uses
#AC_CHECK_LIB([gslcblas], [main])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([modf])

# define warning flags 
# AM_CXXFLAGS = @WARNING_FLAGS@ 
# or mumble_CXXFLAGS = @WARNING_FLAGS@ can then be used in any Makefile.am 
WARNING_FLAGS="-Wall -pedantic"
AC_SUBST(WARNING_FLAGS) 

# Doxygen Settings
DX_HTML_FEATURE(ON)
DX_PDF_FEATURE(OFF)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN([mrs], [Doxyfile])

AM_INIT_AUTOMAKE

AC_CONFIG_FILES([Makefile
src/Makefile
examples/Makefile
examples/AIA/Makefile
examples/AIA/Exm_3_3/Makefile
examples/AIA/Exm_3_4/Makefile
examples/AIA/Exr_11_33/Makefile
examples/AIA/Exr_11_35/Makefile
examples/CXSC/Makefile
examples/GOpt/Makefile
examples/GOpt/Levy/Makefile
examples/GOpt/PhyloPOT/Makefile
examples/GOpt/Rosenbrock/Makefile
examples/MooreRejSam/Makefile
examples/MooreRejSam/CFN3/Makefile
examples/MooreRejSam/IsIt1or2Coins/Makefile
examples/MooreRejSam/BinomialPartition/Makefile
examples/MooreRejSam/Levy/Makefile
examples/MooreRejSam/PhyloPOT/Makefile
examples/MooreRejSam/Rosenbrock/Makefile
examples/MooreRejSam/Shiryaev/Makefile
examples/MooreRejSam/WitchesHat/Makefile
examples/newSubPav/Makefile
examples/newSubPav/Exm_3_3/Makefile
examples/newSubPav/Exm_3_4/Makefile
examples/newSubPav/Exr_11_33/Makefile
examples/newSubPav/Exr_11_35/Makefile
examples/StatsSubPav/Makefile
examples/StatsSubPav/Averaging/Makefile
examples/StatsSubPav/BivGaussian/Makefile
examples/StatsSubPav/PointClouds/Makefile
examples/StatsSubPav/CVOptMap/Makefile
examples/StatsSubPav/Levy/Makefile
examples/StatsSubPav/Uniform/Makefile
examples/StatsSubPav/COPERRandAIC/Makefile
examples/StatsSubPav/HistReport/Makefile
examples/StatsSubPav/CatalanCoeff/Makefile
examples/StatsSubPav/MCMC/Makefile
examples/StatsSubPav/MultiManager/Makefile
examples/StatsSubPav/MCMCHistTimingAndErrors/Makefile
examples/MappedTargets/Makefile
examples/MappedSP/Makefile
examples/MappedSP/Example1D/Makefile
examples/MappedSP/Example2D/Makefile
examples/MappedSP/ExampleSine/Makefile
examples/MappedSP/ExampleColour2D/Makefile
examples/MappedSP/ExampleColour3D/Makefile
examples/MappedSP/ExampleVector3D/Makefile
examples/MappedSP/ExampleOsc1D/Makefile
examples/MappedSP/Levy/Makefile
examples/MappedSP/Rosen/Makefile
examples/MappedSP/SpatialObjectRepresentation/Makefile
examples/MixtureTargets/Makefile
tests/Makefile
tests/TestingMappedSPs/Makefile
tests/TestingRealMapped/Makefile
tests/TestingFunctionEstimator/Makefile
tests/TestingPCF/Makefile
tests/TestingIrregularShapes/Makefile
tests/TestingNonMapped/Makefile
tests/TestingNonMapped/TestHist/Makefile
tests/TestingNonMapped/TestHistColl/Makefile
tests/TestingNonMapped/TestSPS/Makefile
tests/TestingNonMapped/TestSPColl/Makefile
tests/TestingNonMapped/TestSPnodes/Makefile
tests/TestingNonMapped/TestHistInput/Makefile
])

# these modules are currently in development
#examples/Voxels/Makefile
#examples/Voxels/ImageSetDiff/Makefile
#examples/MooreImpSam/Makefile
#examples/MooreImpSam/CesaroBellNo/Makefile
#examples/MooreImpSam/CFN3/Makefile
#examples/MooreImpSam/Rosenbrock/Makefile
#)

AC_OUTPUT
