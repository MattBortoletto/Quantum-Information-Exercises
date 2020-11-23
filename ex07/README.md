Matteo Bortoletto

Quantum information and computing - Exercise 6

# Time-dependent Schrödinger equation

This program computes the time evolution for of the ground state for an Hamiltonian with shifted potential.  

## Contents
- Ex07-MatteoBortoletto.f90: Fortran code;
- time_evol.f90: Fortran code that computes the time evolution using the split-operator method;
- Ex07-MatteoBortoletto-REPORT.pdf: report written in LaTeX;
- Ex07.x: executable file for the Fortran code;
- plot.gnu: gnuplot script to plot the time evolution of the probability density;
- plot_mean.gnu: gnuplot script to plot the mean position of the probability density;
- plot_prob.gnu: gnuplot script to plot the probability densities;
- gif5.gnu: gnuplot script to generate a gif of the probability density time evolution; 
- compile.sh: bash script that compiles and executes the Fortran programs;
- run_plot.py: Python script that runs the gnuplot scripts to plot the data;
- debugger.f90: Fortran module for debugging;
- util.f90: Fortran module with some utilities;
- fftw3.f03: FFTW Fortran code;
- Output files:
    - pr_mean_time_evol_T.txt: mean position of the probability density
    - prob_time_evol_T.txt: probability time evolution
    - psi_imag_time_evol_T.txt: psi imaginary part time evolution
    - psi_real_time_evol_T.txt: psi real part time evolution
    - V_time_evol_T.txt: potential time evolution
    where T is the value of T. 
  Output plots:
    - evolution_T.gif: gif of the time evolution of the probability density and of the potential 
    - mean_T.pdf: plot of the mean position in time 
    - prob_T.pdf: plot of the probability time evolution

## Documentation
The Ex07-MatteoBortoletto.f90 file uses four modules: 
- 'HarmonicOscillator1D' which contains the following functions/subroutines:
    - DiscretizedLapalacian: computes the discretized laplacian;
    - HarmonicPotential: computes the harmonic potential;
    - ComputeEigenvalues: computes the eigenvalues/eigenvectors using the LAPACK subroutine 'zheev';
    - ComputeProbabilityDens: computes the probability densities corresponding to one eigenfunction;
    - ComputeThEnergy: computes the theoretical eigenvalues.

- 'Utilities', which contains the following functions/subroutines:
    - str_i: converts an integer into string;
    - str_r_d: converts a real into string with decimal notation;
    - str_r_e: converts a real into string with exponential notation;
    - WriteEigenvalues: writes the computed eigenvalues in a text file;
    - WriteEigenvectors: writes the computed eigenvectors in a text file;
    - NormalizePs: normalizes the psi and check the normalization.

- 'TimeEvolution', which contains the following functions/subroutines:
    - FourierTransform: computes the Fourier transform using FFTW
    - InverseFourierTransform: computes the inverse Fourier transform using FFTW
    - psiTimeEvol: computes the time evolution of the state

- 'debugger', which contains a module with a subroutine that can be used as checkpoint for debugging. It includes a control on a logical variable 'debug' and optional variable/array/message to be printed. Furthermore, it contains a logical variable 'end_program' which can stop the program in case of need. 

The main program ('time_evolution') asks the user to enter the value of T. Then the computation proceeds automatically.