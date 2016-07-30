noinst_PROGRAMS = TestHists

TestHists_SOURCES = TestHists.cpp HistArithmetic.cpp ../testing_tools.cpp ../testing_tools.hpp ../dataprep.hpp
EXTRA_DIST = testing_hists.sh
TestHists_CPPFLAGS = -I${top_builddir}/src -I../
TestHists_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
