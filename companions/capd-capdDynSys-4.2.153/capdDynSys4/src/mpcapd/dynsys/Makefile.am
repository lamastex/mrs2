include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libmpcapddynsys.la
libmpcapddynsys_la_SOURCES=C2Solver.cpp  CnSolver.cpp  DynSys.cpp  FadTaylor.cpp Solver.cpp
libmpcapddynsys_la_CXXFLAGS= ${AM_CXXFLAGS}  ${mpcapddynsys_CFLAGS} ${mpcapdAlg_CFLAGS} 



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =
