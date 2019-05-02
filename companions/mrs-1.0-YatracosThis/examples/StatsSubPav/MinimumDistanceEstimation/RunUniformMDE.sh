#!/bin/bash
#File: RunUniformMDE.sh

DATASEED=1
D=5
N=15000
HOLDOUTPERCENT=0.33
MAXLEAVESEST=10000 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

./UniformMDE $DATASEED $D $N $HOLDOUTPERCENT $MAXLEAVESEST $CRITLEAVES $NUM_CHECKS $NUM_ITERS



