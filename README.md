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
* To run the binary, just type "./main <filename>" in the command line without the ".2dm" at the end

### Example parallel run with 4 nodes on a personal computer
* go to CMAKE/
* run "cmake -DUSE_MPI=ON .."
* run "make"
* go to CMAKE/bin/
* run "mpirun -np 4 ./main simple_square"
* Checkout the results
