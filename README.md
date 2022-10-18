# Parallel Quicksort Algorithm Using the Message Passing Paradigm

This program will generate a random sequence of size n * 100 (n being the # of processors) and perform a parallel quicksort algorithm

### Compiling on an MPI cluster

		 mpicxx child.cpp -o c -std=c++14

### Running the program on the cluster

		mpirun -np 32 c
