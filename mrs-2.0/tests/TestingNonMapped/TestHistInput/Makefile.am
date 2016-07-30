noinst_PROGRAMS = TestInput
TestInput_SOURCES = input_test.cpp

EXTRA_DIST = inputFile1d_illegalformatting.txt \
				inputFile2d_extratabs.txt \
				inputFile2d_ints.txt \
				inputFile2d.txt \
				inputFile1d_ints.txt \
				inputFile2d_headers1.txt 
				inputFile2d_legalformatting.txt \
				inputFile3d1.txt \
				inputFile1d.txt \
				inputFile2d_headers2.txt \
				inputFile2d_missingelements1.txt \
				inputFile3d2.txt \
				inputFile2d_all_illegalformatting.txt \
				inputFile2d_headers3.txt \
				inputFile2d_missingelements2.txt \
				inputFile3d_missingelements1.txt \
				inputFile2d_blanklines1.txt \
				inputFile2d_headers4.txt \
				inputFile2d_missingelements3.txt \
				inputFileDosCR1.txt \
				inputFile2d_blanklines2.txt \
				inputFile2d_illegalformatting.txt \
				inputFile2d_missingelements_and_illegalformatting.txt

TestInput_CPPFLAGS = -I${top_builddir}/src

TestInput_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
