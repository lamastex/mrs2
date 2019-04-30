#!/bin/bash
#File: RunMDETest.sh

INPUTFILENAME="testdat.txt" #file name
CRITLEAVES=100 #split until this number of leaves in the PQ
NUM_CHECKS=10 #collate the num_check-th histogram
NUM_ITERS=5 #number of iterations for "zooming-in" 

./MDETest $INPUTFILENAME $CRITLEAVES $NUM_CHECKS $NUM_ITERS
