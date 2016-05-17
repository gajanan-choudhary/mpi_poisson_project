# README #

### What is this repository for? ###

* This program solves the Poisson equation in Parallel. This project was done as a part of the course: Parallel Computing for Science and Engineering, The University of Texas at Austin, U.S.A.
* Version 0.1

### How do I get set up? ###

* The build directory is <source_dir>/CMAKE/
* To compile, go to <source_dir>/CMAKE/ and run :
     SERIAL:    "cmake -DUSE_MPI=OFF ..", followed by "make"
     PARALLEL:  "cmake -DUSE_MPI=ON  ..", followed by "make"
* The binary is built in the file <source_dir>/CMAKE/bin with the name "main"
* To run the binary, just type "./main <filename>" in the command line without the ".2dm" and "_part_#.2dm" at the end for serial and parallel runs, respectively

### Example parallel run with 4 nodes on a personal computer ###
* go to CMAKE/
* run "cmake -DUSE_MPI=ON .."
* run "make"
* go to CMAKE/bin/
* run "mpirun -np 4 ./main simple_square"
* Checkout the results


### Folder list ###
* main/     contains the main program
* include/  contains the basic header files
* structs/  contains the node, element, and model structs
* initio/   contains files for reading the input in serial or parallel
* solver/   contains the conjugate gradient solver function and helper functions
* testcase/ contains 2 serial programs "create_mesh" and "partition_mesh" for creating and partitioning a square domain
* CMAKE/    is the folder intended for compiling the program (you may make another folder)
* CMAKE/bin contains 2 sample input files for serial and parallel runs
