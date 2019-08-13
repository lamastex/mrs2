#!/bin/bash
#File: RunTrajectory.sh
#Executes "Trajectory"

#A txt file that contain the filenames of trajectory data
FILENAMES=trajectory_names

#For animal migration data:
#The birds have about a 2 metre wing span and are about 80 cm long. 
#They fly within about 50 m of the sea surface, and mostly within 2 m of it.
#So the volume of the object is 2/1000 * 80/100/1000 = 0.0000016
VOL=0.0000016;

./Trajectory ${FILENAMES}.txt ${VOL}
