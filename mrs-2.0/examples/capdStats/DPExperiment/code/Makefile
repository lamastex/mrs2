#    A simple Makefile for the pendulum program.
#
#    Author: Warwick Tucker <warwick@math.uu.se>
#    Latest edit: Thu Mar 27 02:00:41 CET 2014

default :
	@-echo ""; exit 0
	@echo "****************** DOUBLE PENDULUM PROJECT ******************"
	@echo "";
	@echo "You must specify 'make target', where the 'target' can be:"
	@echo ""
	@echo "    pendulum     (The main file...)"
	@echo ""
	@echo "    tidy         (tidy up the directory)"
	@echo ""
	@echo "    clean        (clean up the directoy)"
	@echo ""
	@echo "";

# -----------------------------------------------------------------------

# The executable files:
PROGRAMS = pendulum

# Absolute CAPD4 binary path:
CAPD = /home/rsa64/all/svn/mrs/branches/raaz/git/mrs2/companions/capd-capdDynSys-4.2.153/bin

# CAPD4 include path.
INCLUDE = `$(CAPD)/capd-config --cflags`

# CAPD4 library path.
LIBS = `$(CAPD)/capd-config --libs`

# -----------------------------------------------------------------------
# CLEAN UP
# -----------------------------------------------------------------------

clean: 
	@ rm -f *.o *~ $(PROGRAMS) gmon* gprof*

# -----------------------------------------------------------------------

tidy:
	@ rm -f *~

# -----------------------------------------------------------------------
# ALL
# -----------------------------------------------------------------------

all: $(PROGRAMS)

# -----------------------------------------------------------------------
# PENDULUM
# -----------------------------------------------------------------------

pendulum: pendulum.cpp
	g++ -O2 $(INCLUDE) $@.cpp $(LIBS) -o $@

# -----------------------------------------------------------------------
# END OF FILE
# -----------------------------------------------------------------------

