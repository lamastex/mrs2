##Configuration
```
./bootstrap
./custom_config.sh
make
```

In `custom_config.sh:`
```
CXSCDIR=/home/gat41/cxsc
MYMRSDIR=/home/gat41/mrs
set -x;
./configure CPPFLAGS="-I${CXSCDIR}/include" LDFLAGS="-L${CXSCDIR}/lib" --prefix="${MYMRSDIR}" $@
LIBRARY_PATH=/usr/local/lib
export LD_LIBRARY_PATH
```

#####################################
In `configure.ac`:

##Uniform:
`Uniform/UnifMDE.cpp`: Application of Devroye, 2001's MDE for density estimation of mixtures 
for the uniform distribution. Returns the minimum distance estimate and the L1-error.

There are currently 4 types of mixtures:

A uniform mixture can be created by ...

*

*


Calls the following methods:

*

*

`RunUnifMDE.sh`: a shell script to run `UnifMDE`:

Syntax: `UnifMDE n d dataSeed mixShape simNum critLeaves maxCheck`

* `n`: sample size (note that the training set will have `n` data points, and the validation set will have `n/2` data points, i.e. the the number of points generated will be `3n/2`.
* `d`: dimensions
* `dataSeed`: seed for generating data  
* `shape`: which mixture to use
* `simNum`: the `simNum`-th simulation
* `critLeaves`: stop splitting when there are `critLeaves` leaves in the SRP
* `maxCheck`:



##MappedMDEGaussian


##MappedMDERosen




##RCode
