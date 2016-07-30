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
