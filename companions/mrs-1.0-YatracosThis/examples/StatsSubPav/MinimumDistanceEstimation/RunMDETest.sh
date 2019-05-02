#!/bin/bash
#File: RunMDETest.sh

#INPUTFILENAME="simulated_rosenbrock_data1.txt" #file name
INPUTFILENAME="simulated_gaussian_data1.txt" #file name

HOLDOUTPERCENT=0.33  # % of data to be held out for validation
#the number of data points to be held out is obtained by multiplying
#the hold out percent with the total number of points, rounded to the nearest integer

CRITLEAVES=100 #split until this number of leaves in the PQ

NUM_CHECKS=10 #collate num_checks histogram

NUM_ITERS=5 #number of iterations for "zooming-in" 

./MDETest $INPUTFILENAME $HOLDOUTPERCENT $CRITLEAVES $NUM_CHECKS $NUM_ITERS
