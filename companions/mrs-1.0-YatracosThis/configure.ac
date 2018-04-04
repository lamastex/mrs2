
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
examples/StatsSubPav/Makefile
examples/StatsSubPav/MinimumDistanceEstimation/Makefile
examples/MappedTargetsTrunk/Makefile
examples/StatsSubPav/MappedSPSampling/Makefile
examples/targets/Makefile
])



# these modules are currently in development
#examples/Voxels/Makefile
#examples/Voxels/ImageSetDiff/Makefile
#examples/MooreImpSam/Makefile
#examples/MooreImpSam/CesaroBellNo/Makefile
#examples/MooreImpSam/CFN3/Makefile
#examples/MooreImpSam/Rosenbrock/Makefile
#examples/StatsSubPav/HullPropagation/Makefile
#)

AC_OUTPUT
