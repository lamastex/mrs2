
### Trajectory: arithmetic for enclosures on spatio-temporal data

*Implementation of* Statistical regular pavings to analyze massive data of aircraft trajectories, Gloria Teng, Kenneth Kuhn and Raazesh Sainudiin, [Journal of Aerospace Computing, Information, and Communication, Vol. 9, No. 1, pp. 14-25, doi: 10.2514/1.I010015, 2012](https://arc.aiaa.org/doi/abs/10.2514/1.I010015)

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
- A txt file that contain the filenames of positional data (1st col: longitude, 2nd col: latitude
- vol of object

To execute program:

```%sh
./RunTrajectory.sh
```

Optional output:
The trajectories for each data set and a collated histogram of these trajectories.

### 2. Dyanmic Trajectories	
	
