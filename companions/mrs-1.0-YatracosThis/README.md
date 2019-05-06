
## Minimum estimation density branch

## gat41's thesis branch from mrs-1

## To make files in path mrs2/companions/mrs-1.0-YatracosThis

```%sh
$ ./bootstrap

$ ./custom_config.sh

$ make
```

Upon successful compilation, the folder ``examples/StatsSubPav/MinimumDistanceEstimation`` will contain the following executables:

```%sh
CheckBounds
MappedGaussianMDE
MappedRosenbrockMDE
MDETest
UniformMDE
```  
Shell scripts have been set up for each executable. Details for each executable are given in the following.  

## 1. MDE for Densities Built using Piecewise Constant Functions
	
### Mapped Gaussian densities: 
```%sh
./RunMappedGaussianMDE.sh
./RunMappedGaussianMDESimulations.sh
```

### Mapped Rosenbrock densities
```%sh
./RunMappedRosenbrockMDE.sh
./RunMappedRosenbrockMDESimulations.sh
```

### Uniform densities
```%sh
./RunUniformMDE.sh
```

## 2. MDE for Input with Text Files
```%sh
./RunMDETest.sh
```


## 3. Check bounds of Theorems 2 and 3
```%sh
./RunCheckBounds.sh
```
