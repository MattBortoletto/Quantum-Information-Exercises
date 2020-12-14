Matteo Bortoletto

Quantum information and computing - Exercise 9

# Quantum Ising Model

This program computes the Hamiltonian for a quantum anti-ferromagnetic Ising model, diagonalizes it and saves the first k eigenvalues.

## Contents
- Ex09-MatteoBortoletto.f90: Fortran code with main program;
- ising.f90: Fortran module which contains utilities for the Ising model;
- Ex09-MatteoBortoletto-REPORT.pdf: report written in LaTeX;
- Ex09.x: executable file for the main program;
- compile.sh: bash script that compiles and executes the main program;
- debugger.f90: Fortran module for debugging;
- util.f90: Fortran module with some utilities;
- plot_cpu.gnu: gnuplot script to plot the computation CPU time as a function of N;
- plot_gs: gnuplot script to plot the ground state for different values of N;
- run_plt.py: python script to run the gnuplot file 'plot_k_eig.gnu' multiple times.

- Output files:
    - eig_N.txt: first k eigenvalues for a certain value of N;
    - cpu_times: cpu times;
  Output plots:
    - gs.pdf: plot of the ground state for different values of N;
    - cpu_times.pdf: plot of the CPU times.

## Documentation
The following modules are used: 
- 'Ising', which contains the following functions/subroutines:
    - IsingHamiltonian: initializes the Hamiltonian
    - Diagonalize: diagonalize the Hamiltonian using the LAPACK zheev subroutine

- 'Utilities', which contains the following functions/subroutines:
    - str_i: converts an integer into string;
    - str_r_d: converts a real into string with decimal notation;
    - str_r_e: converts a real into string with exponential notation;
    - WriteEigenvalues: writes the computed eigenvalues in a text file;
    - WriteEigenvectors: writes the computed eigenvectors in a text file;
    - NormalizePsi: normalizes the psi and check the normalization.

- 'debugger', which contains a module with a subroutine that can be used as checkpoint for debugging. It includes a control on a logical variable 'debug' and optional variable/array/message to be printed. Furthermore, it contains a logical variable 'end_program' which can stop the program in case of need. 

The main program (contained in Ex09-MatteoBortoletto.f90) is structured as follows. We use two nested 'for' loops: the first running on different values of $N$ and the second on different values of $\lambda$. Inside the inner loop we initialize the Hamiltonian, we diagonalize it and we save the first $k$ eigenvalues. Then, for each $N$, we save the results in different output file which are named according to the value of $N$.