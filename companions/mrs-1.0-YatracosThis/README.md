
## Minimum estimation density branch

## gat41's thesis branch from mrs-1

## To make files in path mrs2/companions/mrs-1.0-YatracosThis

```%sh
$ ./bootstrap

$ ./custom_config.sh

$ make
```
## Examples
### MDE for Mapped Gaussian Densities
```%sh
$ cd examples/StatsSubPav/MinimumDistanceEstimation

$ vim RunMappedGaussianMDE.sh
------Adjust parameters accordingly-----------------
#!/bin/bash
#File: RunMapedGaussianMDE.sh
DATASEED=1 #seed use to generate data
D=1 #dimensions
N=100 #number of points
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate the num_check-th histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

./MappedGaussianMDE $DATASEED $D $N $MAXLEAVESEST $CRITLEAVES $NUM_CHECKS $NUM_ITERS
------------------------------------------------------

$ ./RunMappedGaussianMDE.sh
```
Example output:
```%sh
./MappedGaussianMDE : process id is 18617
Data seed is 1

Make the function estimator to 100 leaves
X	10.00000	1.48672E-006	0.39894	-5.00000	5.00000

pq down to max leaves 119
Number of leaves in estimate: 119 s.
After split, getTotalAreaOfIntervalBand() =   0.041837
Computing time for pq split in estimate: 0.48 s.
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
100 training data and 50 validation data inserted.
Increment by : 10
Perform 5 iterations

Iteration 0......
Calling prioritySplitAndEstimate...
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
Calling prioritySplitAndEstimate...
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
Calling prioritySplitAndEstimate...
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
Computing time for MDE: 0.28 s.
  0.465737	22	  0.268294	6
Error computations output to results1.txt
```
The file `results1.txt` consists of 4 values:
 - The IAE of the histogram estimate using MDE
 - the number of leaves of histogram estimate using MDE
 - The IAE of the histogram estimate which gives the minimum IAE
 - the number of leaves of the histogram estimate with the minimum IAE 

### Simulations for Mapped Gaussian Densities
The script `RunMappedGaussianMDESimulations.sh`can be used to run repeated trials of the MDE. 
```%sh
$ cd examples/StatsSubPav/MinimumDistanceEstimation

$ vim RunMappedGaussianMDESimulations.sh
------Adjust parameters accordingly-----------------
#!/bin/bash
#File: RunMapedGaussianMDESimulations.sh

rm *.txt #be careful not to remove txt files that you want to keep

NUM_SIMS=3; #How many simulations
DIM=1 #dimensions
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate the num_check-th histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

for N in 100 #number of data points to generate
do
   #mkdir ${MYDIR}/n${N}critleaves${CRITLEAVES}
	for DATASEED in `seq 1 ${NUM_SIMS}`
		do 
		echo Simulation ${NUM_SIMS} for n = $N and CRITLEAVES = ${CRITLEAVES}
		./MappedGaussianMDE ${DATASEED} ${DIM} $N ${MAXLEAVESEST} ${CRITLEAVES} ${NUM_CHECKS} ${NUM_ITERS}  
		#mv *.txt ${MYDIR}/n${I}L${L}
		#append the results for each loop to an overall file #make sure in the correct folder depending on n and L
		cat results${DATASEED}.txt >> n${N}critleaves${CRITLEAVES}all.txt
		rm results${DATASEED}.txt
	done
done
----------------------------------

$ ./RunMappedGaussianMDESimulations.sh
```
An example of the output txt file for 3 simulations:
```%sh 
$ cat n100critleaves100all.txt 
  0.465737	22	  0.268294	6
  0.247111	14	  0.227635	11
  0.426323	17	  0.262360	8
```

### MDE for Mapped Rosenbrock Densities

### MDE for Uniform Densities

### Input .txt files
