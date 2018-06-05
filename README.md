# nullspace-basis
Compute a sparse nullspace basis of a rectangular matrix

## Context
The nullspace basis needs to be sparse, well-conditioned, easy to apply. Those criteria are
difficult to match together, actually it can be rather ill-conditioned or dense. Paramount
among them is sparsity. The problem of finding a sparse nullspace basis is shown to be
NP-hard and even If the rectangular matrix is sparse that do not mean the sparsity of the nullbasis.

Here I developed a technique that suits those purposes, which is the sparse Gaussian elimination
approach that attempts to preserve sparsity while keeping rounding errors under control.
We compute the nullspace basis by performing LU on the matrix transpose of the rectangular matrix. 
This latter is called a skinny matrix, since its rows outnumbers its columns. As done above, there exists two
permutation matrices, P used for stability, and Q used for sparsity such that we can define the nullbasis Z :

![alt text](screenshot.png?raw=true "")

### How to use ?
This directory contains the Fortran nullspace basis computation algorithm using C routines in SuperLU.

To compile the examples, type:
```
make
```
To run the examples, type:
```
nullbasis < matC
```

### Input & Output 
The input is a sparse rectangular matrix in Harwell-Boeing format.
You can obtain L and U factor matrices using a MATLAB format (i, j, value) or a Harwell-Boeing format.

### Algorithm steps 

For a sparse rectangular matrix C. 
* LU factorization of tr(C)
* Bulding the inverse permutation
* Getting tr(L1) and tr(L2)
* Computing tr(L1)^-1 tr(L2)
* Storing the nullspace Z
