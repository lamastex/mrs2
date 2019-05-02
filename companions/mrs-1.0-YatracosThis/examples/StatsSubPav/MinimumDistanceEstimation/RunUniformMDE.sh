#!/bin/bash
#File: RunUniformMDE.sh

DATASEED=1
D=1
MIXSHAPE="1,1"
#MIXSHAPE="2,2,1"
#MIXSHAPE="3,3,2,1"
#MIXSHAPE="3,4,4, 2, 2, 3, 3"

N=10
HOLDOUTPERCENT=0.33
MAXLEAVESEST=100 #maximum number of leaves in the function estimator
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate num_checks histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

./UniformMDE $DATASEED $D $MIXSHAPE $N $HOLDOUTPERCENT $MAXLEAVESEST $CRITLEAVES $NUM_CHECKS $NUM_ITERS



