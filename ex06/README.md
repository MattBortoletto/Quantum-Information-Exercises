Matteo Bortoletto

Quantum information and computing - Exercise 6

# One-dimensional quantum harmonic oscillator

This program computes the (first k) eigenvalues and eigenvectors for the quantum harmonic oscillator in D=1.  

## Contents
- Ex06-MatteoBortoletto.f90: Fortran code;
- th_eigenf.f90: Fortran code that computes the theoretical eigenfunctions;
- Ex06-MatteoBortoletto-REPORT.pdf: report written in LaTeX
- Ex06.x: executable file for the Fortran code;
- plot_eigenval.gnu: gnuplot script to plot the eigenvalues;
- plot_eigenvect.gnu: gnuplot script to plot the eigenvectors;
- plot_prob.gnu: gnuplot script to plot the probability densities;
- compile.sh: bash script that compiles and executes the Fortran programs;
- run_plot.py: Python script that runs the gnuplot scripts to plot the data;
- debugger.f90: Fortran module for debugging;
- Output files (* = txt/pdf/log)): the notation is the following
    - data_L_Dx_omega_m_hbar.txt: text file with the numerical eigenfunctions;
    - data_t_k_L_Dx_omega_m_hbar.txt: text file with the first k theoretical eigenfunctions.
  where data = [en, ef, pr] = [energies, eigenfunctions, probability densities].
  Output plots:
  - eigenvalues.pdf: plot of the energies;
  - first_k_eigenfunc.pdf: plot of the first k numerical eigenfunctions;
  - first_k_eigenfunc_th.pdf: plot of the first k theoretical eigenfunctions;
  - first_k_prob.pdf: plot of the first k numerical probability densities;
  - first_k_prob_th.pdf: plot of the first k theoretical probability densities.

## Documentation
The debugger.f90 file contains a module with a subroutine that can be used as checkpoint for debugging. It includes a control on a logical variable 'debug' and optional variable/array/message to be printed. Furthermore, it contains a logical variable 'end_program' which can stop the program in case of need. 

The Ex06-MatteoBortoletto.f90 file contains two modules: 
- 'HarmonicOscillator1D' which contains the following functions/subroutines:
    - DiscretizedLapalacian: computes the discretized laplacian;
    - HarmonicPotential: computes the harmonic potential;
    - ComputeEigenvalues: computes the eigenvalues/eigenvectors using the LAPACK subroutine 'zheev';
    - ComputeProb: computes the probability densities corresponding to one eigenfunction.

- 'Utilities', which contains the following functions/subroutines:
    - str_i: converts an integer into string;
    - str_r_d: converts a real into string with decimal notation;
    - str_r_e: converts a real into string with exponential notation;
    - WriteEigenvalues: writes the computed eigenvalues in a text file;
    - WriteEigenvectors: writes the computed eigenvectors in a text file.

The main program ('harmonic_oscillator_1D') asks the user if he want to enter custom parameters or use the default ones. If he decides for the former, he will be asked to enter L, Dx, omega, m and hbar. According to these choices the results will be saved in different output files. 
The program th_eigenf.x can be used to compute the theoretical eigenfunctions. It works exactly as the main program.
Once the Fortran code is run, a Python script ('run_plot.py') will automatically run the gnuplot scripts to plot the results.