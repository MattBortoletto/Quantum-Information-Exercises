Matteo Bortoletto

Quantum information and computing - Exercise 8

# Density matrices

This program deals with density matrices. It initializes a pure state for a composite system, computes the density matrix and computes the reduced density matrices in the case of N=2 subsystems.

## Contents
- ex08.f90: Fortran code with main program;
- ex08test.f90: Fortran test code with a N=2, D=2 system;
- ex08_eff_test.f90: Fortran code to test the efficiency;
- density_matr.f90: Fortran code which contains the module to deal with density matrices;
- Ex07-MatteoBortoletto-REPORT.pdf: report written in LaTeX;
- Ex08.x: executable file for the main program;
- Ex08Test.x: executable file for the test;
- Ex08eff.x: executable file for the efficiency test;
- cpu_plot.gnu: gnuplot script to plot the state initialization CPU time as a function of N both for separable and generic states;
- compile.sh: bash script that compiles and executes the main program;
- compile_test.sh: bash script that compiles and executes the test program;
- compile_eff.sh: bash script that compiles and executes the efficiency test program;
- debugger.f90: Fortran module for debugging;
- util.f90: Fortran module with some utilities;

- Output files:
    - rho.txt: density matrix;
    - rhoA.txt: reduced density matrix for subsystem A;
    - rhoB.txt: reduced density matrix for subsystem B;
    - cpu_times: initialization cpu times (N, sep_time, gen_time);
  Output plots:
    - cpu_times.pdf: plot of the CPU times.

## Documentation
The following modules are used: 
- 'DensityMatrices' which contains the following functions/subroutines:
    - InitPureSepState: initializes a pure separable state
    - InitPureGenState: initializes a pure general state
    - ComplexSquareMatrixTrace: computes the trace of a square complex matrix
    - BuildSeparableStateFromMatr: builds a separable state starting from the coefficient matrix;
    - BuildSeparableStateFromVect: builds a separable state starting from the coefficient vector;
    - ReducedDensityMatrixB: computes the reduced density matrix for subsystem B;
    - ReducedDensityMatrixA: computes the reduced density matrix for subsystem A;
    - WriteComplex16Matr: writes a complex*16 matrix in a text file.

- 'Utilities', which contains the following functions/subroutines:
    - str_i: converts an integer into string;
    - str_r_d: converts a real into string with decimal notation;
    - str_r_e: converts a real into string with exponential notation;
    - WriteEigenvalues: writes the computed eigenvalues in a text file;
    - WriteEigenvectors: writes the computed eigenvectors in a text file;
    - NormalizePsi: normalizes the psi and check the normalization.

- 'debugger', which contains a module with a subroutine that can be used as checkpoint for debugging. It includes a control on a logical variable 'debug' and optional variable/array/message to be printed. Furthermore, it contains a logical variable 'end_program' which can stop the program in case of need. 

The main program (ex08.f90) initializes both a separable and generic state, computes the density matrix and check that its trace is 1. The program in ex08test.f90 considers a state psi = (|00> + |01>) / sqrt(2) and computes the density matrix and the two reduced density matrices for the subsystems. The program in ex08_eff_test.f90 initializes new states for D=2 and varying N until there is an error, measuring the CPU time. 