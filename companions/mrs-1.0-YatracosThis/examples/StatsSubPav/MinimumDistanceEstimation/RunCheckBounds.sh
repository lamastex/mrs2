#!/bin/bash
#File: RunCheckBounds.sh

#rm *.txt
NUM_SIMS=5; #How many simulations
D=2
N=15000
HOLDOUTPERCENT=0.33
MAXLEAVESEST=10000 #maximum number of leaves in the function estimator
CRITLEAVES=1000 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

for DATASEED in `seq 1 ${NUM_SIMS}`
	do 
	echo Simulation ${NUM_SIMS} for n = $N and CRITLEAVES = ${CRITLEAVES}
	./CheckBounds $DATASEED $D $N $HOLDOUTPERCENT $MAXLEAVESEST $CRITLEAVES $NUM_CHECKS $NUM_ITERS
	cat theorem2_check${DATASEED}.txt >> Gaussian_d${DIM}_n${N}critleaves${CRITLEAVES}theorem2.txt
	cat theorem_values${DATASEED}.txt >> Gaussian_d${DIM}_n${N}critleaves${CRITLEAVES}theorem_values.txt
done


