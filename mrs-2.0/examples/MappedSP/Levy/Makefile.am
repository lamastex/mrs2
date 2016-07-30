noinst_PROGRAMS = LevyDensity2D LevyDensity2D_PaperExample

TARGETS = ../../MappedTargets

LevyDensity2D_SOURCES = LevyDensity2D.cpp \
						$(TARGETS)/LevyDensityFobj2D.cpp \
						$(TARGETS)/LevyDensityFobj2D.hpp

LevyDensity2D_PaperExample_SOURCES = LevyDensity2D_PaperExample.cpp \
						$(TARGETS)/LevyDensityFobj2D.cpp \
						$(TARGETS)/LevyDensityFobj2D.hpp

#EXTRA_DIST =
AM_CPPFLAGS = -I${top_builddir}/src -I$(TARGETS)
LDADD = -lcxsc ${top_builddir}/src/libmrs.a
