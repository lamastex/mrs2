## Project Structure

This started as Dillon Geroge's honours project at University of Canterbury in 2016.
The project had to come to a stop due to situations beyond control.
The project was supervised by Raazesh Sainudiin and Kourosh Neshatian on the basis of 
discussions between Luc Devroye and Raazesh Sainudiin during the end of 2015 during Luc's 
visit to Canterbury.

Raazesh Sainudiin and Dillon George, Sat Jul 30 23:28:40 CEST 2016.

The include directory contains all header files for the project.
With the source `cpp` files contained in the `src/` directory.
```
include/
└── DensityTree
    ├── PointUtils.h
    ├── bsp.h
    └── CGALTypeDefs.h
```
```
src/
├── bsp.cpp
├── pointUtils.cpp
```
## Build instructions
To build this first ensure that a copy of CGAL is installed.

Then it should be a matter of:

```
mkdir build/
cd build

cmake ../
make all #if CMake is set to generate Makefile
```


## Overview of Files and their intended use

### `include/DensityTree/CGALTypeDefs.hpp`
Contains typedefs for working with CGAL points. These typdefs are essentially the same as those used in
CGAL examples for the dD Geometry Kernel.

### `src/bsp.cpp` & `include/DensityTree/bsp.hpp`
Contains the definition of the space partitioning tree.

`BSP_node` represents a node in the tree. It has pointers to its left and right children, are null when no
children exist. It has a counter keeping track of the total number of points contained within that partition.
It stores a vector of those points enclosed within that partition.

The basic outline of how the tree is constructed:
	- Create the root node and insert it into a priority queue
	- While the largest (in number of points) node in the queue has over a specified threshold of points
	  within it:
		- Generate splitting hyperplane splitting the node region into left and right regions
		- insert those new left and right nodes into the priority queue

__TODO__: More detailed description of the actual tree construction algorithm

### `PointUtils.cpp/hpp`
Utilities for generating random points. Currently wraps the stdlib random number generators to generate
either uniform or normally distributed points. But the general method of generating points could easily be
extended other sources of randomness


### `hyperplane.cpp`

Simple application/example for generating random points and then constructing a bsp tree from them.
