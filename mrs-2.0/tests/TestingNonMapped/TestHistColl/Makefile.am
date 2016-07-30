noinst_PROGRAMS = TestCollHists
TestCollHists_SOURCES = TestCollatorHists.cpp ../testing_tools.cpp ../testing_tools.hpp ../dataprep.hpp 
EXTRA_DIST = import_test_badbox1.txt \
			import_test_badbox2.txt \
			import_test_badbox3.txt \
			import_test_badnodelevels1.txt \
			import_test_badnodelevels2.txt \
			import_test_extrarange.txt \
			import_test_missingrange.txt \
			import_test_nogap.txt \
			testing_collhists.sh

TestCollHists_CPPFLAGS = -I${top_builddir}/src -I../
TestCollHists_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
