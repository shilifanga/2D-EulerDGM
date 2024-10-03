This is an MPI code. You need the MPI environment. 
To run the MPI code, you need to do the following steps: 
(1) first step, you need to compile it and get an executable file "all".  The command is mpiifort -r8 all.f90 -o all
(2) Second step, to run the executable file "all".  The command is mpirun -np 16 ./all 

This code can solve the two-dimensional Euler equations, primarily using the Runge-Kutta discontinuous Galerkin method. To control numerical spurious oscillations, the jump filter is employed.

The jump filter in the discontinuous Galerkin method for hyperbolic conservation laws.
