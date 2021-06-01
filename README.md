# SBD-failures
Code to accompany "Failure of the SBD technique applied to complete and cluster synchronization of generic networks"

## Description
There are two ways to use the SBD library.
1. An executable called runSBD can be used to compute the simultaneous block diagonalizing transformation.
2. An API is available, in C, that can be integrated into larger programs.

Documentation for the API can be generated with doxygen by running 'make doc' 

### runSBD
The executable requires a working directory that contains a file named 'graph' that has the graph of interest written as a list of edges using 0-based indexing.
One can use either the graph's adjacency matrix or its Laplacian as the first matrix, and can also use either the orbit indicator matrices or the minimum balanced coloring matrices, as the remaining matrices for the SBD.

#### Usage Examples
Assume you have a directory named example1 that already contains a file named graph, and you want to compute the SBD of the Laplacian and the orbit indicator matrices and write both the P matrix and the B matrices to files

./runSBD -O -L -P -B ./example1

If instead you want to use the adjacency matrix, the minimal balanced coloring, and only write the transformation matrix

./runSBD -M -A -P ./example1


### API
The API allows additional freedom to compute the SBD of a set of M matrices generated in any way, for instance two or more different graphs.
Helper functions are available, such as a qgraph type to generate graphs with desired orbits of the automorphism group.

## Dependencies
The following libraries are used and variables to point to the appropriate headers/libraries should be set in the make.var file
1. igraph
2. arpack
3. lapack/blas
4. nauty

## Installation
The file 'make.var' is used to define the variables for the make build system. The C compiler and any compiling flags as well as the locations of the header files for igraph and nauty and the library locations for all of the dependencies.
