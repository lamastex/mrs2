
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

1. [MDE for Densities Built using Piecewise Constant Functions](https://github.com/lamastex/mrs2/blob/master/companions/mrs-1.0-YatracosThis/README.md#1-mde-for-densities-built-using-piecewise-constant-functions)

2. [MDE for Input with Text Files](https://github.com/lamastex/mrs2/blob/master/companions/mrs-1.0-YatracosThis/README.md#2-mde-for-input-with-text-files)

3. [Check bounds of Theorems 2 and 3](https://github.com/lamastex/mrs2/blob/master/companions/mrs-1.0-YatracosThis/README.md#3-check-bounds-of-theorems-2-and-3)
---
---
## 1. MDE for Densities Built using Piecewise Constant Functions
	
The programs ``MappedGaussianMDE``, ``MappedRosenbrockMDE``, and ``UniformMDE`` will build the corresponding density object (``fobj``) which is then stored as a `PiecewiseConstantFunction` object.  To output the object, remove the comments at lines **162** and **169**.

Data is generated from this density object and stored in a `RVecData` container. 
The data is input into an `AdaptiveHistogramValidation` object. Toggle comments at lines **192 - 204** for output of the data set. 

The Yatracos set and Delta_theta values for the histogram estimates built using this data set are then computed by calling `AdaptiveHistogramValidation::prioritySplitAndEstimate`. 

The shell scripts `RunMappedGaussianMDE.sh` or `RunMappedRosenbrockMDE.sh` or `RunUniformMDE.sh` allows us to input the following parameters needed to run the code:

 - `DATASEED`: value of the seed for the random number generator
 - `D`: dimension of the data set to be generated
 - `N`: number of data points to be generated
 - `HOLDOUTPERCENT`: % of data to be held out.
 - `MAXLEAVESEST`: integer, maximum number of leaves in the function estimator for the `fobj` object
 - `CRITLEAVES`: this should an integer, and is the stopping criteria for the PQ. The splitting will stop when this number of leaves is reached.
 - `NUM_CHECKS`: see [(2)](https://github.com/lamastex/mrs2/blob/master/companions/mrs-1.0-YatracosThis/README.md#2-mde-for-input-with-text-files) for details.
 - `NUM_ITERS`: see [(2)](https://github.com/lamastex/mrs2/blob/master/companions/mrs-1.0-YatracosThis/README.md#2-mde-for-input-with-text-files) for details.

**Output:** 

 - `iaes_and_diffs{DATASEED}.txt`: there are 4 values in this txt file - 
	 - the IAE value corresponding to the minimum Delta_theta value
	 - the minimum IAE value
	 - the absolute difference between the values in first and second columns
 - `sequence{DATASEED}.txt`: the number of leaf nodes / thetas used to obtain the Yatracos set and thus Delta_theta values
 - `deltas{DATASEED}.txt`: the Delta_theta values corresponding to sequence.
 - `iaes{DATASEED}.txt`: the IAE values of each histogram obtained at each split starting from the root node to `CRITLEAVES`. There will `CRITLEAVES` IAE values.

Note: to （un)suppress any of these output files, toggle the corresponding comments (see lines **325 - 370**).

### An example using ``MappedGaussianMDE``
```%sh
$ vim RunMappedGaussianMDE.sh
------Adjust parameters accordingly-----------------
#!/bin/bash
#File: RunMappedGaussianMDE.sh

DATASEED=1
DIM=1
N=150
HOLDOUTPERCENT=0.33
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

./MappedGaussianMDE $DATASEED $D $N $HOLDOUTPERCENT $MAXLEAVESEST $CRITLEAVES $NUM_CHECKS $NUM_ITERS
------------------------------------------------------

$ ./RunMappedGaussianMDE.sh 
./MappedGaussianMDE : process id is 8864
Data seed is 1

Make the function estimator to 100 leaves
pq down to max leaves 119
Number of leaves in estimate: 119 s.
After split, getTotalAreaOfIntervalBand() =   0.041837
Computing time for pq split in estimate: 0.59 s.
Hull propagation
Priority merge to 100 leaves
Computing time for hull propagate and merge up in estimate: 0 s.
After propagation and priority merge, getTotalAreaOfIntervalBand() =   0.050950
number of leaves is = 100
Making estimate and normalising
estimate has integral   0.999871 before normalizing
estimate has integral   1.000000

Generating data for simulation
Computing time for simulating data: 0 s.
150 points generated

Running minimum distance estimation with hold-out...
50 points held out.
Increment by : 10
Perform 5 iterations

Iteration 0......
Calling prioritySplitAndEstimate for mapped functions...
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
Calling prioritySplitAndEstimate for mapped functions...
---- Hist 1-----
---- Hist 10-----
---- Hist 20-----
---- Hist 22-----
---- Hist 24-----
---- Hist 26-----
---- Hist 28-----
---- Hist 30-----
---- Hist 32-----
---- Hist 34-----
---- Hist 36-----
---- Hist 38-----
---- Hist 40-----
---- Hist 50-----
---- Hist 60-----
---- Hist 70-----
---- Hist 80-----
---- Hist 90-----
---- Hist 100-----

Run MDE with the final sequence...
Calling prioritySplitAndEstimate for mapped functions...
---- Hist 1-----
---- Hist 10-----
---- Hist 20-----
---- Hist 22-----
---- Hist 23-----
---- Hist 24-----
---- Hist 25-----
---- Hist 26-----
---- Hist 28-----
---- Hist 30-----
---- Hist 32-----
---- Hist 34-----
---- Hist 36-----
---- Hist 38-----
---- Hist 40-----
---- Hist 50-----
---- Hist 60-----
---- Hist 70-----
---- Hist 80-----
---- Hist 90-----
---- Hist 100-----
Computing time for MDE: 0.2 s.
The minimum max delta is 0.1675 at 22 leaf nodes with IAE  0.465737
The minimum IAE is  0.300462 at 10 leaf nodes.
Main results output to gaussian_1d_150n_iaes_and_diff1.txt
```

### Performing Simulations
The scripts `RunMappedGaussianMDESimulations.sh`/ `RunMappedRosenbrockMDESimulations.sh` / `RunUniformMDESimulations.sh`can be used to run repeated trials of the MDE. 

The output files are same as the above, indexed by the data seed. The `iaes_and_diff` files will be appended to one text file `n${N}critleaves${CRITLEAVES}all.txt` for the ease of comparison.

```%sh
$ vim RunMappedGaussianMDE.sh
------Adjust parameters accordingly-----------------
#!/bin/bash
#File: RunMappedGaussianMDESimulations.sh

rm *.txt #be careful not to remove txt files that you want to keep

NUM_SIMS=10; #How many simulations
HOLDOUTPERCENT=0.33
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

for DIM in 1 #may need more CRITLEAVES for higher dims 
do
	for N in 150 1500 
	do
		for DATASEED in `seq 1 ${NUM_SIMS}`
			do 
			echo Simulation ${NUM_SIMS} for n = $N and CRITLEAVES = ${CRITLEAVES}
			./MappedGaussianMDE ${DATASEED} ${DIM} $N ${HOLDOUTPERCENT} ${MAXLEAVESEST} ${CRITLEAVES} ${NUM_CHECKS} ${NUM_ITERS}  
			cat gaussian_${DIM}d_${N}n_iaes_and_diff${DATASEED}.txt >> gaussian_${DIM}d_${N}n_iaes_and_diff_${NUM_SIMS}sims.txt
			rm gaussian_${DIM}d_${N}n_iaes_and_diff${DATASEED}.txt
		done
		echo Simulation results appended to gaussian_${DIM}d_${N}n_iaes_and_diff_${NUM_SIMS}sims.txt
	done
done
----------------------------------

$ ./RunMappedGaussianMDESimulations.sh
#the console output will be similar to ./RunMappedGaussianMDE.sh
```

An example of `gaussian_1d_150n_iaes_and_diff_5sims.txt` for 5 simulations:
```%sh 
$ cat gaussian_1d_150n_iaes_and_diff_5sims.txt 
  0.465737	  0.300462	  0.165275
  0.247111	  0.230395	  0.016716
  0.426323	  0.363088	  0.063235
  0.422476	  0.232552	  0.189924
  0.363403	  0.363403	  0.000000
```
--- 
## 2. MDE for Input with Text Files
The program ``MDETest`` takes in data from a text file and outputs a histogram estimate using minimum distance estimation with a hold-out scheme.

The columns of the data file should separated by white space (space or tabs), with no non-numeric characters and integers. Two example data files ``simulated_gaussian_data1.txt`` (*n* = 150 and *d* = 1) and ``simulated_rosenbrock_data1.txt``  (*n* = 150 and *d* = 2)  are available in the folder for code testing.

A root box based on the data set will be built. If you want to define your own root box, insert a comment to line **88**. Then,  remove the comments at lines **78 - 84** in ``MDETest.cpp``  and input your desired values for the root box.  The code at lines 78 - 84 will build a root box with equal widths at all dimensions. For example, if the given values are. a data set with 2 dimensions will be given the root box [-5, 5] x [-5, 5]. Compile the code before execution.

(Note: in the future, hopefully, a boolean can handle whether or not we want to build a root box based on the data, or build a self-defined root box.)

The shell script ``RunMDETest.sh`` allows us to input the parameters needed to run the code:

 - ``INPUTFILENAME``: the name of the txt file that contains the data set
 
 - ``OUTPUTFILENAME``: a filename to be associated with the txt files that will be output from the program
 
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
 - ``theta_and_delta_theta.txt``: there are two columns: the number of leaf nodes / thetas that were used and the corresponding Delta_theta values for each histogram in the final sequence
 - ``mdehist.txt``: the histogram estimate built using the number of leaf nodes / theta* that minimizes Delta_theta.
 - ``coverage.txt``: the tail probabilities for each data point

Note: to suppress any of these output files, comment the corresponding lines **216 - 243**.

```%sh
$ vim RunMDETest.sh
------Adjust parameters accordingly--------
#!/bin/bash
#File: RunMDETest.sh

INPUTFILENAME="simulated_data_for_mdetest/simulated_gaussian_data1.txt"
OUTPUTFILENAME="gaussian"
HOLDOUTPERCENT=0.33  
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

./MDETest $INPUTFILENAME $OUTPUTFILENAME $HOLDOUTPERCENT $CRITLEAVES $NUM_CHECKS $NUM_ITERS
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
Computing time for MDE: 0.099367 s.
The minimum max delta is 0.1625 at 14 leaf nodes.
The minimum max delta is 0.1625 at 14 leaf nodes.
MDE histogram output to gaussian_mdehist.txt
Coverage values output to gaussian_coverage.txt
The delta_theta values and corresponding number of leaf nodes (theta) output to gaussian_theta_and_delta_theta.txt
```

---
## 3. Check bounds of Theorems 2 and 3
The program `CheckBounds` calls the function `AdaptiveHistogramValidation::getMDETheoremValues` to obtain the required values to check the bounds of Theorems 2 and 3 of the JJSD paper (in revision). Data generated from  `fobj` density objects will be used to check the bounds. The current set up is such that data points are generated from a standard uniform distribution. To use another density object, edit line *100* and the root box (line *96*) if necessary.

The required values for checking bounds:

 - `IAEforMinDelta`:  the IAE corresponding to the minimum distance estimate (for Theorems 2 and 3)
 - `minIAE`: the minimum IAE for the histogram built using `(n - m)`  data points (for Thoerem 2)
 - `maxDelta`: the max Delta value (for Theorem 2)
 - `minIAEAllPoints`: the minimum IAE for the histogram built using `n` data points (for Theorem 3) 
 
The shell script `RunCheckBounds.sh` allows us to input the parameters needed to run the code, These parameters are similar to those in [(1)](https://github.com/lamastex/mrs2/blob/master/companions/mrs-1.0-YatracosThis/README.md#1-mde-for-densities-built-using-piecewise-constant-functions).


**Output:** 

 - `theorem2_check.txt`
	 - For Theorem 2, we need to check that `IAEforMinDelta <= 3*minIAE + 4maxDelta`. 
	 - The text file has 3 columns: 
	 `(IAEforMinDelta) (3*minIAE + 4maxDelta) (Boolean for IAEforMinDelta <= 3*minIAE + 4maxDelta)`
	 
- `theorem3_iaes.txt`:  this txt file has 2 columns and will be appended to an overall file for the ease of computing the bounds of Theorem 3. The 2 columns are 
	- The IAE of the minimum distance estimate
	- The minimum IAE value based on `n` data points 
	 
 - `delta_theta_and_iaes{DATASEED}.txt`: this text file has 3 columns:
	 - the delta_theta values for each theta
	 - the IAE values based on `n-m` data points for each theta
	 - the IAE values based on `n` data points for each theta
 
 - `delta{DATASEED}.txt`: the Delta values (Theorem 2).

Note: toggle the corresponding comments to suppress/output the files needed. See lines *318 - 374*.

**Example**
```%sh
#!/bin/bash
#File: RunCheckBounds.sh

#Executes MDETest which will produce uniform random variates to obtain 
#the values needed to check the bounds of Theorems 2 and 3.

rm *.txt

NUM_SIMS=10; #How many simulations
HOLDOUTPERCENT=0.33
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ (this will also be Theta)
NUM_CHECKS=10 #collate num_checks histogram

for DIM in 3  
do
  for N in 1500
  do
    for DATASEED in `seq 1 ${NUM_SIMS}`
    do 
    	./CheckBounds $DATASEED $DIM $N $HOLDOUTPERCENT $MAXLEAVESEST $CRITLEAVES $NUM_CHECKS
	cat uniform_${DIM}d_${N}n_theorem2_check${DATASEED}.txt >> uniform_${DIM}d_${N}n_theorem2_check_${NUM_SIMS}sims.txt
	cat uniform_${DIM}d_${N}n_theorem3_iaes${DATASEED}.txt >> uniform_${DIM}d_${N}n_theorem3_iaes_${NUM_SIMS}sims.txt
	
	rm uniform_${DIM}d_${N}n_theorem2_check${DATASEED}.txt
	rm uniform_${DIM}d_${N}n_theorem3_iaes${DATASEED}.txt
    done
	echo Results for Theorem 2 appended to uniform_${DIM}d_${N}n_theorem2_check_${NUM_SIMS}sims.txt
	echo Results for Theorem 3 appended to uniform_${DIM}d_${N}n_theorem3_iaes_${NUM_SIMS}sims.txt
  done
done
```

All text files will be indexed by the data seed. In particular, the files `theorem2_check` and `theorem3_iaes` will each be appended to an "overall" file.

- The checks for Theorem 2 for `NUM_SIMS` simulations can be found in `uniform_1d_1500n_theorem2_check_5sims.txt `.
```%sh
$ cat uniform_3d_1500n_theorem2_check_5sims.txt 
  0.160199	  5.417866	1
  0.000000	  1.606439	1
  0.104104	  2.050189	1
  0.100124	  3.654924	1
  0.067164	  0.538889	1
```
- The values needed for Theorem 3 over `NUM_SIMS` simulations can be found in `uniform_3d_1500n_theorem3_iaes_5sims.txt `. We can either use a program or spreadsheet program to get the average values of the first and second columns, and perform the necessary computations to check the bounds. 
- Note that the average of the first column corresponds to the LHS of Theorem 3, while the average of the second column corresponds to the RHS.

```%sh
$ cat uniform_3d_1500n_theorem3_iaes_5sims.txt 
  0.160199	  0.000000
  0.000000	  0.000000
  0.104104	  0.000000
  0.100124	  0.000000
  0.067164	  0.000000
```
