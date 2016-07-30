noinst_PROGRAMS = TestCollNodes
TestCollNodes_SOURCES = TestCollatorNodes.cpp ../testing_tools.cpp ../testing_tools.hpp ../dataprep.hpp 
EXTRA_DIST = testing_colls.sh
TestCollNodes_CPPFLAGS = -I${top_builddir}/src -I../
TestCollNodes_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
