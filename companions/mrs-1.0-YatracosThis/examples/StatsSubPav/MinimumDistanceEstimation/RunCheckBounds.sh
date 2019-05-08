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

for DIM in 1 2 5 10 100 1000 #may need more CRITLEAVES for higher dims 
do
	for N in 150 1500 15000 150000 1500000 15000000
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

