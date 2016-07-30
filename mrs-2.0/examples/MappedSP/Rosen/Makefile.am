noinst_PROGRAMS = Rosen \
					RosenMultipleD \
					RosenMultipleD_Width

TARGETS = ../../MappedTargets

Rosen_SOURCES = Rosen.cpp \
				$(TARGETS)/RosenFobj.cpp $(TARGETS)/RosenFobj.hpp
RosenMultipleD_SOURCES = RosenMultipleD.cpp \
				$(TARGETS)/RosenFobj.cpp $(TARGETS)/RosenFobj.hpp
				
RosenMultipleD_Width_SOURCES = RosenMultipleD_Width.cpp \
				$(TARGETS)/RosenFobj.cpp $(TARGETS)/RosenFobj.hpp

#EXTRA_DIST =
AM_CPPFLAGS = -I${top_builddir}/src -I$(TARGETS)
LDADD = -lcxsc ${top_builddir}/src/libmrs.a

