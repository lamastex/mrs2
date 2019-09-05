#!/bin/bash
#File: RunDynamicTrajectory.sh
#Executes "DynamicTrajectory"

VOL=0.1 #note: can be put into the for loop if we know each individual craft size
STARTTIME=0;
TOTALTIMEBLOCK=4;
NUMBEROFFLIGHTS=2;

./DynamicTrajectory ${VOL} ${STARTTIME} ${TOTALTIMEBLOCK} ${NUMBEROFFLIGHTS}
