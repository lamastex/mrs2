#!/bin/bash
#File: RunUniformMDESimulations.sh

rm *.txt #be careful not to remove txt files that you want to keep

NUM_SIMS=10; #How many simulations
HOLDOUTPERCENT=0.33
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

for N in 150 1500 15000 150000 1500000 15000000
do
	for DIM in 1 2 5 10 100 1000 
	do
		for DATASEED in `seq 1 ${NUM_SIMS}`
			do 
			echo Simulation ${NUM_SIMS} for n = $N and CRITLEAVES = ${CRITLEAVES}
			./UniformMDE ${DATASEED} ${DIM} $N ${HOLDOUTPERCENT} ${MAXLEAVESEST} ${CRITLEAVES} ${NUM_CHECKS} ${NUM_ITERS}  
			cat uniform_${DIM}d_${N}n_iaes_and_diff${DATASEED}.txt >> uniform_${DIM}d_${N}n_iaes_and_diff_${NUM_SIMS}sims.txt
			rm uniform_${DIM}d_${N}n_iaes_and_diff${DATASEED}.txt
		done
		echo Simulation results appended to uniform_${DIM}d_${N}n_iaes_and_diff_${NUM_SIMS}sims.txt
	done
done


