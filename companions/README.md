To prepare for mrs-2.0 from companion librarys' source tarballs:

```%sh
$ cd companions/tarballs/

$ tar zxf cxsc-2-5-4.tar.gz 
$ mv cxsc-2-5-4 ../

$ tar zxf capd-capdDynSys-4.2.153.tar.gz 
$ mv capd-capdDynSys-4.2.153 ../

$ tar zxf gsl-2.1.tar.gz 
$ mv gsl-2.1 ../

$ cd ..
```

To configure and make capd-capdDynSys-4.2.153
=============================================

```%sh
cd capd-capdDynSys-4.2.153/
$ ./configure
$ make
```
We will use the path to this directory in Makefiles needing capd later: `...git/mrs2/companions/capd-capdDynSys-4.2.153/'


To configure, make and make install gsl-2.1
===========================================

```%sh
$ cd gsl-2.1/
$ vi README 
$ vi INSTALL 
$ pwd
# the output of pwd should be the directory path after "=" below:
$ ./configure --prefix=.../git/mrs2/companions/gsl-2.1
$ make
$ make install
----------------------------------------------------------------------
Libraries have been installed in:
   .../git/mrs2/companions/gsl-2.1/lib

If you ever happen to want to link against installed libraries
in a given directory, LIBDIR, you must either use libtool, and
specify the full pathname of the library, or use the `-LLIBDIR'
flag during linking and do at least one of the following:
   - add LIBDIR to the `LD_LIBRARY_PATH' environment variable
     during execution
   - add LIBDIR to the `LD_RUN_PATH' environment variable
     during linking
   - use the `-Wl,-rpath -Wl,LIBDIR' linker flag
   - have your system administrator add LIBDIR to `/etc/ld.so.conf'

See any operating system documentation about shared libraries for
more information, such as the ld(1) and ld.so(8) manual pages.
----------------------------------------------------------------------
```

To configure, make and make install c-xsc - see Installation (classical version) in README
===========================================

```%sh
$ cd cxsc-2-5-4/
$ vi README
```

Copy the `pwd` output of this directory for pasting later during `./install_cxsc`.

```
$ pwd 
.../git/mrs2/companions/cxsc-2-5-4

$ ./install_cxsc

##
##  CXSC is a C++ library for eXtended Scientific Computing
##
##  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
##                          Universitaet Karlsruhe, Germany
##            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
##                          Universitaet Wuppertal, Germany
##
##  This library is free software; you can redistribute it and/or
##  modify it under the terms of the GNU Library General Public
##  License as published by the Free Software Foundation; either
##  version 2 of the License, or (at your option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##  Library General Public License for more details.
##
##  You should have received a copy of the GNU Library General Public
##  License along with this library; if not, write to the Free
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##

Type 'yes' to accept this license offer.
Type 'no' to decline this license offer.

Do you accept the terms of either license? [no] yes

---------------------------------------------------------------------------------

checking host system type....i686
checking build system type...Linux

Configuring for a i686-pc-linux-gnu host.

Section: Compiler
=================
Possible values:    gnu    for the GNU C/C++ Compiler
                    intel  for the INTEL C/C++ Compiler
Which C++ compiler? [gnu] 
Would you like to select a special GNU C/C++ Compiler-Version? [no] 
Would you like to use compiler optimization options? [yes] 
Which optimization level (e.g. -O1, -O2, ...)? [-O3] 

Using GNU V.4.8 Compiler


Section: Rounding Operations
============================
Possible values: asm    for hardware support for IEEE 754 arithmetic implemented in Assembler
                 hard   for hardware support for IEEE 754 arithmetic
                 soft   for software emulations for directed rounded floating-point operations
Which method? [asm] 

cxscconf.h succesfully created


Section: Prefix
===============
(un-)installation prefix
e.g. /usr/local/cxsc or local home directory
Prefix [~/cxsc] .../git/mrs2/companions/cxsc-2-5-4 ## paste pwd output here!

...
...

Section: Library
================
Do you want to create a dynamic or static library? [dynamic]/static/both both

Create a dynamic and static library.


Section: Automatic Code Generation
==================================
Proceed with automatic code generation? [yes]/no yes

...
...

***********************************************************
*                      FINAL RESULT                       *
***********************************************************

Number of equalities :128
Number of subsets    :0
Number of supersets  :16
Other relations      :0

*********************************************************
*    All results are equal to the expected solutions.   *
*  It seems the libraries (C-XSC and CToolbox) are OK!  *
*********************************************************
make[1]: Leaving directory `/home/rsa64/all/svn/mrs/branches/raaz/git/mrs2/companions/cxsc-2-5-4/CToolbox'

Section: CXSC Results
=====================
Library successfully created
Examples successfully created
Library successfully installed


Section: CToolbox Results
=========================
Examples successfully created
Library and programs successfully installed and tested

```
