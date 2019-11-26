# Parallel inverted dual time stepping method (PIDTS)

PIDTS is an open-source software for parallel-in-time implementation of the dual time stepping (DTS) method, focusing on fluid flow and heat transfer problems. The DTS method contains two loops: the outer loop sweeps all physical time steps and the inner loop iteratively solve the implicit equations at each time step. PIDTS is designed by inverting the outer loop and the inner loop within the framework of the DTS method. This change brings a concurrency of the solving at different time steps, which is the core idea behind PIDTS. PIDTS is a non-intrusive code with a flexible framework for DTS scenarios, i.e. users only need write simple wrappers for their spatial-discretization codes to provide PIDTS with the implicit iterative solver.
 
A user-defined wrapper of implicit iterative solver is required to integrate with PIDTS for a specific problem. To facilitate user's easy understanding,we add a steady/unsteady natural convection problem as an example, which are discretized with the Chebyshev pseudospectral method for spatial terms and the fully backward difference scheme for temporal terms. The employed implicit iterative solver is the fourth-order four-stage Runge-Kutta integrator in pseudo time. 
 
This is a scientific project. If you use PIDTS for publications or presentations in science, please support the project by citing our publications given in [references](REFERENCE.md). In case you have question regarding PIDTS, don't hesitate to contact us <wenqianchen2016@gmail.com>.

# Installation
The example code is compiled with gfortran 5.4, and has been tested for Centos 7 and Ubuntu 16.  In addition, the external librarie, namely MPI is used in PIDTS.
## Compiling the code
* Open a terminal 
* Change into the PIDTS directory
* Create a new subdirectory and use CMake to configure and compile the code.

             mkdir build; cd build
             cmake ../
             make
## Run the code
* Open a terminal 
* Change into the PIDTS directory
* Make sure there exist two directories (INPUT and OUTPUT), and creat a INPUT.txt in the directory INPUT.
* Execute command: " mpiexec -n np ./build/bin/PIDTS", where np is used to set the number of CPU cores for parallel computing.

Besides, We have built a batch file (batch.sh) to help test the parallel performance against various parameters. If you use it, just use chmod to add executable permissions by "chmod +x batch.sh" and then execute it by "./batch.sh".