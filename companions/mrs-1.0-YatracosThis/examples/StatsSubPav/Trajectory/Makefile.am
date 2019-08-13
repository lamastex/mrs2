noinst_PROGRAMS = Trajectory DynamicTrajectory
Trajectory_SOURCES = Trajectory.cpp
DynamicTrajectory_SOURCES = DynamicTrajectory.cpp

#EXTRA_DIST =
Trajectory_CPPFLAGS = -I${top_builddir}/src -I../../targets
Trajectory_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
DynamicTrajectory_CPPFLAGS = -I${top_builddir}/src -I../../targets
DynamicTrajectory_LDADD = -lcxsc ${top_builddir}/src/libmrs.a
