noinst_PROGRAMS = TestSPS
TestSPS_SOURCES = TestSPSnodes.cpp ../testing_tools.cpp ../testing_tools.hpp ../dataprep.hpp 

EXTRA_DIST = testing_sps.sh \
				Hist1SPS_Expected.txt \
				testSPSnodeByNodeOutputExpected.txt \
				unionOfHist1AndHist2PavingsExpected.txt

TestSPS_CPPFLAGS = -I${top_builddir}/src -I../
TestSPS_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
