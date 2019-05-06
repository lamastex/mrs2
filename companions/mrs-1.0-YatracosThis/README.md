
## Minimum estimation density branch

## gat41's thesis branch from mrs-1

## To make files in path mrs2/companions/mrs-1.0-YatracosThis

```%sh
$ ./bootstrap

$ ./custom_config.sh

$ make
```

Upon successful compilation, the folder ``examples/StatsSubPav/MinimumDistanceEstimation`` will contain the following executables:

```%sh
CheckBounds
MappedGaussianMDE
MappedRosenbrockMDE
MDETest
UniformMDE
```  
Shell scripts have been set up for each executable. Details for each executable are given in the following.  

1. ...

2. [MDE for Input with Text Files]([https://github.com/lamastex/mrs2/blob/master/companions/mrs-1.0-YatracosThis/README.md#2-mde-for-input-with-text-files])

3. ...
---

## 1. MDE for Densities Built using Piecewise Constant Functions
	
### Mapped Gaussian densities: 
```%sh
./RunMappedGaussianMDE.sh
./RunMappedGaussianMDESimulations.sh
```

### Mapped Rosenbrock densities
```%sh
./RunMappedRosenbrockMDE.sh
./RunMappedRosenbrockMDESimulations.sh
```

### Uniform densities
```%sh
./RunUniformMDE.sh
```

## 2. MDE for Input with Text Files
The program ``MDETest`` takes in data from a text file and outputs a histogram estimate using minimum distance estimation with a hold-out scheme.

The columns of the data file should separated by white space (space or tabs), with no non-numeric characters and integers. Two example data files ``simulated_gaussian_data1.txt`` (*n* = 150 and *d* = 1) and ``simulated_rosenbrock_data1.txt``  (*n* = 150 and *d* = 2)  are available in the folder for code testing.

A root box based on the data set will be built. If you want to define your own root box, insert a comment to line **88**. Then,  remove the comments at lines **78 - 84** in ``MDETest.cpp``  and input your desired values for the root box.  The code at lines 78 - 84 will build a root box with equal widths at all dimensions. For example, if the given values are. a data set with 2 dimensions will be given the root box [-5, 5] x [-5, 5]. Compile the code before execution.

(Note: in the future, hopefully, a boolean can handle whether or not we want to build a root box based on the data, or build a self-defined root box.)

The shell script ``RunMDETest.sh`` allows us to input the parameters needed to run the code:

 - ``INPUTFILENAME``: the name of the txt file that contains the data set
 
 - ``HOLDOUTPERCENT``: % of data to be held out. 
	 - The number of data points to be held out is obtained by multiplying the hold out percent with the total number of points, rounded to the nearest integer.
 
 - ``CRITLEAVES``: this should an integer, and is the stopping criteria for the PQ. The splitting will stop when this number of leaves is reached.

 - ``NUM_CHECKS``: also an integer. 
	 - This is the number of histograms that will be collated  in the first iteration, along with the histogram without any splits (i.e. the root box). A Yatracos set based on these histograms will be obtained. 
	 - For example, if ``NUM_CHECKS = 10`` and ``CRITLEAVES = `100`, then ``100/10 = 10`` histograms will be collated in the first iteration, along with the histogram without any splits. The sequence of histograms obtained will be: 
```%sh		 
		 Hist 1 (histogram with 1 leaf node)
		 Hist 10 (histogram with 10 leaf nodes)
		 Hist 20 (histogram with 20 leaf nodes)
		 .
		 .
		 .
		 Hist 100 (histogram with 100 leaf nodes) 
```

 - ``NUM_ITERS``: the number of iterations for "zooming-in" 
	 - After the first iteration, ``NUMCHECKS + 1`` histograms are collated and the corresponding Yatracos set obtained. Delta_theta values for each histogram are stored in the vector ``vecMaxDelta``. The number of leaf nodes (also theta) corresponding to the *k* = 3 lowest Delta_theta values are obtained.  The value of *k* can be changed in line **131**.  
	 - In the example below, the 3 lowest Delta_theta values corresponding to the histograms with 10, 20 and 30 leaf nodes. Thus, in the second iteration, we will "zoom in" to obtain ``NUM_CHECKS = 10`` histograms between 10 and 30 leaf nodes will be collated with increments of  ``(30-10)/10 = 2``, along with the other histograms from the first iteration.
	 - The "zooming-in" action will stop after ``NUM_ITERS`` is reached or when the increment value becomes 1.  In the example given, we only needed to perform 1 iteration to obtain a increment value of 1, giving a sequence of
```%sh
Hist 1
Hist 10
Hist 11
Hist 12
Hist 13
Hist 14
Hist 16
Hist 18
Hist 20
Hist 22
Hist 24
Hist 26
Hist 28
Hist 30
Hist 40
.
.
.
---- Hist 100-----
```

**Output:**
 - ``deltas.txt``: the Delta_theta values for each histogram in the final sequence
 - ``sequence.txt``: the number of leaf nodes / thetas that were used
 - ``mdehist.txt``: the histogram estimate built using the number of leaf nodes / theta* that minimizes Delta_theta.

Note: to suppress any of these output files, comment the corresponding lines **211 - 234**.

```%sh
$ vim RunMDETest.sh
------Adjust parameters accordingly--------
#!/bin/bash
#File: RunMDETest.sh

#INPUTFILENAME="simulated_rosenbrock_data1.txt" 
INPUTFILENAME="simulated_gaussian_data1.txt"
HOLDOUTPERCENT=0.33  # % of data to be held out for validation
#the number of data points to be held out is obtained by multiplying
#the hold out percent with the total number of points, rounded to the nearest integer
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate this number of histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

./MDETest $INPUTFILENAME $HOLDOUTPERCENT $CRITLEAVES $NUM_CHECKS $NUM_ITERS
------------------------------------------

$ ./RunMDETest.sh
Processing file simulated_gaussian_data1.txt
End of reading data input file: 150 valid data points read in
Holding out 50 data points
A box is being made for the data.  The box is
[ -2.223954,  2.011870]

Running minimum distance estimation with hold-out...
50 points held out.
Increment by : 10
Perform 5 iterations

Iteration 0......
Calling prioritySplitAndEstimatePlain...
---- Hist 1-----
---- Hist 10-----
---- Hist 20-----
---- Hist 30-----
---- Hist 40-----
---- Hist 50-----
---- Hist 60-----
---- Hist 70-----
---- Hist 80-----
---- Hist 90-----
---- Hist 100-----

Iteration 1......
Calling prioritySplitAndEstimatePlain...
---- Hist 1-----
---- Hist 10-----
---- Hist 12-----
---- Hist 14-----
---- Hist 16-----
---- Hist 18-----
---- Hist 20-----
---- Hist 22-----
---- Hist 24-----
---- Hist 26-----
---- Hist 28-----
---- Hist 30-----
---- Hist 40-----
---- Hist 50-----
---- Hist 60-----
---- Hist 70-----
---- Hist 80-----
---- Hist 90-----
---- Hist 100-----

Run MDE with the final sequence...
Calling prioritySplitAndEstimatePlain...
---- Hist 1-----
---- Hist 10-----
---- Hist 11-----
---- Hist 12-----
---- Hist 13-----
---- Hist 14-----
---- Hist 16-----
---- Hist 18-----
---- Hist 20-----
---- Hist 22-----
---- Hist 24-----
---- Hist 26-----
---- Hist 28-----
---- Hist 30-----
---- Hist 40-----
---- Hist 50-----
---- Hist 60-----
---- Hist 70-----
---- Hist 80-----
---- Hist 90-----
---- Hist 100-----
Computing time for MDE: 0.099367 s.
The minimum max delta is 0.1625 at 14 leaf nodes.
MDE histogram output to mdehist.txt
Sequence of histograms output to sequence.txt
Delta theta values output to deltas.txt
```


## 3. Check bounds of Theorems 2 and 3
```%sh
./RunCheckBounds.sh
```
