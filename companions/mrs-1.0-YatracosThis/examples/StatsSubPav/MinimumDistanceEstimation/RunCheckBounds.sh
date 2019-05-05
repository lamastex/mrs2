#!/bin/bash
#File: RunCheckBounds.sh

NUM_SIMS=10; #How many simulations
DATASEED=1
D=1
N=150
HOLDOUTPERCENT=0.33
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=10 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

for DATASEED in `seq 1 ${NUM_SIMS}`
	do 
	echo Simulation ${NUM_SIMS} for n = $N and CRITLEAVES = ${CRITLEAVES}
	./CheckBounds $DATASEED $D $N $HOLDOUTPERCENT $MAXLEAVESEST $CRITLEAVES $NUM_CHECKS $NUM_ITERS
	cat theorem2_check${DATASEED}.txt >> n${N}critleaves${CRITLEAVES}all.txt
done


