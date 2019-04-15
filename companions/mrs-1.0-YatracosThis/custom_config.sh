#!/bin/sh

#Edit this to match your own location

CXSCDIR=/home/gat41/cxsc
MYMRSDIR=/home/gat41/mrs


set -x;
./configure CPPFLAGS="-I${CXSCDIR}/include" LDFLAGS="-L${CXSCDIR}/lib" --prefix="${MYMRSDIR}" $@

#./configure CPPFLAGS="-isystem${CXSCDIR}/include -I${GSLDIR}/include" LDFLAGS="-L${CXSCDIR}/lib -L${GSLDIR}/lib" --prefix="${MYMRSDIR}" $@


##gat41's configuration settings
LD_LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH
