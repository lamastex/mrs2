noinst_PROGRAMS = Averaging
Averaging_SOURCES = Averaging.cpp ../BivGaussian/dataprep.hpp
EXTRA_DIST = BGAvBoxes.m BGHistsBoxes.m EvalGaussPDF.m
Averaging_CPPFLAGS = -I${top_builddir}/src -I../BivGaussian
Averaging_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
