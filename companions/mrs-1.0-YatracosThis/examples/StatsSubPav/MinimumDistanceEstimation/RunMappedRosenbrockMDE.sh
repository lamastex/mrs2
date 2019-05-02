#!/bin/bash

#File: RunMappedRosenbrockMDE.sh

DATASEED=1
D=2
N=150
HOLDOUTPERCENT=0.33
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

./MappedRosenbrockMDE $DATASEED $D $N $HOLDOUTPERCENT $MAXLEAVESEST $CRITLEAVES $NUM_CHECKS $NUM_ITERS
