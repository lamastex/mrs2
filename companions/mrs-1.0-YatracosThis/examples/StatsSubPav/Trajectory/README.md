
### Trajectory: arithmetic for enclosures on spatio-temporal data

Implementation of *Statistical regular pavings to analyze massive data of aircraft trajectories*, Gloria Teng, Kenneth Kuhn and Raazesh Sainudiin, [Journal of Aerospace Computing, Information, and Communication, Vol. 9, No. 1, pp. 14-25, doi: 10.2514/1.I010015, 2012](https://arc.aiaa.org/doi/abs/10.2514/1.I010015)

```%sh
$ make
```

Upon successful compilation, the folder will contain the following executables:

```%sh
Trajectory
DynamicTrajectory
```  

### 1. Single Trajectories and Collation of Trajectories
Input:
- A txt file that contain the filenames of positional data (1st col: longitude, 2nd col: latitude)
- vol of object

To execute program:

```%sh
./RunTrajectory.sh
```

Optional output:
The trajectories for each data set and a collated histogram of these trajectories. (uncomment the relevant lines in ``Trajectory.cpp``.)

Example:
```%sh
$ ./RunTrajectory.sh
Reading in file names for simulated data: 
0804LongLat.txt
0904LongLat.txt
1004LongLat.txt
1004LongLat.txt

Put all data in a container to get rootbox: 
0804LongLat.txt
End of reading data input file: 320 valid data points read in
0904LongLat.txt
Warning: adding to existing data - mixing datasets
End of reading data input file: 640 valid data points read in
1004LongLat.txt
Warning: adding to existing data - mixing datasets
End of reading data input file: 472 valid data points read in


A box is being made for the data.  The box is 
[1.666054E+002,1.708125E+002]  [-48.143752,-45.980408]  
Data has 2 dimensions.
getRootBoxVol
Vol: 1.6e-06	approxMinVol: 2.16987e-06
================1======================
Processing file 0904LongLat.txt
The file 0904LongLat.txt has 641 lines in it

Getting enclosure for this trajectory: 
End of reading data input file: 640 valid data points read in
Computing time : 0.07 s.

 Add into collator
================2======================
Processing file 1004LongLat.txt
The file 1004LongLat.txt has 473 lines in it

Getting enclosure for this trajectory: 
End of reading data input file: 472 valid data points read in
Computing time : 0.05 s.

 Add into collator
================3======================
Processing file 1004LongLat.txt
The file 1004LongLat.txt has 473 lines in it

Getting enclosure for this trajectory: 
End of reading data input file: 472 valid data points read in
Computing time : 0.05 s.

 Add into collator
The output of the accumulated AdaptiveHistograms has been written to coll.txt
```

### 2. Dynamic Trajectories	
Input:
- VOL: vol of object
- STARTTIME: starting time
- TOTALTIMEBLOCK: how many time blocks
- NUMBEROFFLIGHTS: the number of flights

To execute program:

```%sh
./RunDynamicTrajectory.sh
```

Optional output:
The collated histograms at each time block. (uncomment the relevant lines in ``DynamicTrajectory.cpp``.)

Example:
```%sh
$ ./RunDynamicTrajectory.sh 
Box is: [-10.000000, 10.000000]
[-10.000000, 10.000000]

getRootBoxVol
craftVol: 0.1	approxMinVol: 0.195312
Processing file Time0Flight1.txt

End of reading data input file: 4 valid data points read in
Processing file Time0Flight2.txt

End of reading data input file: 1 valid data points read in
The output of the accumulated AdaptiveHistograms has been written to spaceColl0.txt

Processing file Time1Flight1.txt

Error in AdaptiveHistogram::readRvectorsFromTxt: no data in input file Time1Flight1.txt
Processing file Time1Flight2.txt

End of reading data input file: 1 valid data points read in
The output of the accumulated AdaptiveHistograms has been written to spaceColl1.txt

Processing file Time2Flight1.txt

Error in AdaptiveHistogram::readRvectorsFromTxt: no data in input file Time2Flight1.txt
Processing file Time2Flight2.txt

End of reading data input file: 1 valid data points read in
The output of the accumulated AdaptiveHistograms has been written to spaceColl2.txt

Processing file Time3Flight1.txt

Error in AdaptiveHistogram::readRvectorsFromTxt: no data in input file Time3Flight1.txt
Processing file Time3Flight2.txt

Error in AdaptiveHistogram::readRvectorsFromTxt: no data in input file Time3Flight2.txt
The output of the accumulated AdaptiveHistograms has been written to spaceColl3.txt
```


