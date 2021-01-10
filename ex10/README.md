Matteo Bortoletto

Quantum information and computing - Exercise 9

# Quantum Ising Model

This program computes ground state of a quantum anti-ferromagnetic Ising model using the Real Space Renormalization Group algorithm (RSRG).

## Contents
- Ex10-MatteoBortoletto.f90: Fortran code with main program;
- ising.f90: Fortran module which contains utilities for the Ising model, including the RSRG algorithm;
- Ex10-MatteoBortoletto-REPORT.pdf: report written in LaTeX;
- Ex10.x: executable file for the main program;
- compile.sh: bash script that compiles and executes the main program;
- debugger.f90: Fortran module for debugging;
- util.f90: Fortran module with some utilities;
- plot_gs: gnuplot script to plot the ground state for different values of N;

- Output files:
    - gs_N*.txt: ground state energy densities for different $\lambda$ (field strength);
    - error_gs_N2.txt: error wrt the theorical value $e=-1$ for $\lambda=0$ and $N=2$
  Output plots:
    - gs_N*.pdf: plot of the ground state for different values of N;
    - gs_diff.pdf: plot of the error wrt the mean field;
    - error_gs_N2.pdf: plot of the error wrt the theorical value $e=-1$ for $\lambda=0$ and $N=2$.


## Documentation
The following modules are used: 
- 'Ising', which contains the following functions/subroutines:
    - IsingHamiltonian: initializes the Hamiltonian;
    - DiagonalizeR: diagonalize the Hamiltonian using the LAPACK dsyev subroutine;
    - RSRG: Real Space Renormalization Group algorithm with fixed number of iterations;
    - RSRG1: Real Space Renormalization Group algorithm with convergence threshold;
    - idmatr: defines a identity matrix of size d;
    - KroneckerProd: computes the tensor product.

- 'debugger', which contains a module with a subroutine that can be used as checkpoint for debugging. It includes a control on a logical variable 'debug' and optional variable/array/message to be printed. Furthermore, it contains a logical variable 'end_program' which can stop the program in case of need. 

The main program (contained in Ex10-MatteoBortoletto.f90) is structured as follows. We use two nested `for` loops: the first running on different initial values of $N$ and the second on different values of $\lambda$. Inside the inner loop we initialize the Hamiltonian and we call the RSRG1 function to compute the ground state. For each value of $N$, we write the results in different output files.