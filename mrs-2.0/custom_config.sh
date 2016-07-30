#!/bin/sh

#Edit this to match your own location
CXSCDIR=/home/rsa64/all/software/cxsc/cxsc-2-5-4
GSLDIR=/home/rsa64/all/software/gsl/gsl-2.1
#CXSCDIR=/usr/local/cxsc
#CXSCDIR=/home/gat41/cxsc
#MYMRSDIR=/users/math/rsa64/mrs
#MYMRSDIR=/home/rsa64/software/mrs
#MYMRSDIR=/home/gat41/mrs
#MYMRSDIR=/Users/jah217/newsvn/mrs
#CXSCDIR=/Users/raazesh/cxsc
MYMRSDIR=/home/wynrs1/software/mrs
#MYMRSDIR=/Users/raazesh/newsvn/mrs

set -x;
#./configure CPPFLAGS="-I${CXSCDIR}/include" LDFLAGS="-L${CXSCDIR}/lib" --prefix="${MYMRSDIR}" $@
#./configure CPPFLAGS="-I${CXSCDIR}/include -I/usr/local/include" LDFLAGS="-L${CXSCDIR}/lib" --prefix="${MYMRSDIR}" $@

#using -system for cxsc include headers avoids compiler warning we cannot do anything about
#./configure CPPFLAGS="-isystem${CXSCDIR}/include -I/usr/local/include" LDFLAGS="-L${CXSCDIR}/lib" --prefix="${MYMRSDIR}" $@
./configure CPPFLAGS="-isystem${CXSCDIR}/include -I${GSLDIR}/include" LDFLAGS="-L${CXSCDIR}/lib -L${GSLDIR}/lib" --prefix="${MYMRSDIR}" $@


##sometimes you need to:
#export DYLD_LIBRARY_PATH=/Users/raazesh/cxsc/lib:$DYLD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/usr/local/cxsc/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/rsa64/all/software/cxsc/cxsc-2-5-4/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/rsa64/all/software/gsl/gsl-2.1/lib:$LD_LIBRARY_PATH
