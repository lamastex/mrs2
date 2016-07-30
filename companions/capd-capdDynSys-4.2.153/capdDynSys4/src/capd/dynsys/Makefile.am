include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapddynsys.la
libcapddynsys_la_SOURCES= DynSys.cpp  FadTaylor.cpp Solver.cpp  SolverHOE.cpp C2Solver.cpp  CnSolver.cpp   
libcapddynsys_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapddynsys_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capddynsys_CFLAGS} ${capdAlg_CFLAGS}



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =
