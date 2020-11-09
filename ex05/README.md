Matteo Bortoletto

Quantum information and computing - Exercise 5

# Random matrix program 

This program performs matrix diagonalization using the LAPACK cheev subroutine and then it computes the distribution of the (normalized) spacings between eigenvalues for different matrix samples. 

## Contents
- Ex05-MatteoBortoletto.f90: Fortran code;
- Ex05.x: executable file for the Fortran code;
- herm_fit.gnu: gnuplot script to fit the spacings distribution of complex hermitian matrices;
- diag_fit.gnu: gnuplot script to fit the spacings distribution of real diagonal matrices;
- compile.sh: bash script that compiles and executes the Fortran program;
- run_fit.py: Python script that runs the gnuplot scripts to fit the data;
- debugger.f90: Fortran module for debugging;
- Output files (* = txt/pdf/log)):
    - diag_glo.*: text file with data, pdf with the plot and log file with the fit results for global-normalized spacings for real diagonal matrices;
    - diag_loc_#.*: text file with data, pdf with the plot and log file with the fit results for local-normalized spacings for real diagonal matrices (# = 1, 5, 10, 50, 100);
    - herm_glo.*: text file with data, pdf with the plot and log file with the fit results for global-normalized spacings for complex hermitian matrices;
    - diag_loc_#.*: text file with data, pdf with the plot and log file with the fit results for local-normalized spacings for complex hermitian matrices (# = 1, 5, 10, 50, 100).

## Documentation
The debugger.f90 file contains a module with a subroutine that can be used as checkpoint for debugging. It includes a control on a logical variable 'debug' and optional variable/array/message to be printed. Furthermore, it contains a logical variable 'end_program' which can stop the program in case of need. 

The Ex05-MatteoBortoletto.f90 file contains one module, 'rand_matrix' which contains the following functions/subroutines
- MatrixInit: initializes the matrix of the chosen type. The matrices are filled with random numbers uniformed distributed between -1 and 1;
- ComputeEigenvalues: computes the eigenvalues of a (complex) matrix using the LAPACK cheev subroutines;
- ComputeSpacings: computes the spacings between the eigenvalues of a matrix;
- ComputeNormSpacings: normalizes the spacings using the global mean;
- ComputeNormSpacingsLocal: normalizes the spacings using the local mean;
- ComputePDF: given the spacings, builds an histogram and normalizes it using the area of the bins;
- ComputeR: computes < r >;
- str: converts an integer into string.

The main program ('RandomMatrix') asks the user to enter the dimension of the matrix, its type and if he/she wants to consider global-normalized or local-normalized spacings. Finally, he/she will be asked to enter the number of samples matrices to generate. According to these choices -- and to the locality level -- the results will be saved in different output files.
Once the Fortran code is run, a Python script ('run_fit.py') will automatically run two different gnuplot scripts ('herm_fit.gnu' for complex hermitian matrices, 'diag_fit.gnu' for diagonal real matrices) which plots the distribution P(s) of the spacings vs. the spacings s. 