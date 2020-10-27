EXERCISE 4 - Matteo Bortoletto

CONTENTS
* Ex04-MatteoBortoletto.f90: code 
* Ex04.x: exacutable file 
* Ex04-MatteoBortoletto-REPORT: report written in LaTeX
* ex04run.py: Python script that compiles and runs the Fortran code   
* ex04fit_run.py: Python script that runs the gnuplot script to fit the CPU times vs. matrix dimension 
* ex04plot.gnu: gnuplot script to plot the CPU times vs. matrix dimension for one method
* ex04fit.gnu: gnuplot script to fit and plot the CPU times vs. matrix dimension for one method
* ex04fit_all.gnu: gnuplot script to fit and plot the CPU times vs. matrix dimension all the three methods
* Documentation.txt: text file which contains the program reference guide 
* latex_plots: folder which contains the code to produce LaTeX plots 

OUTPUT FILES    
* ex04run.py: 
    - matrix_dimensions.txt: text file with the dimension of the two matrices to multiply
    - not-optimized.txt: CPU times for the non-optimized matrix multiplication function 
    - optimized.txt: CPU times for the optimized matrix multiplication function 
    - matmul.txt: CPU times for the intrisic matrix multiplication function 
* ex04fit_run.py:
    - CPU_time_fit_not-optimized.png: plot with fit of the CPU times for the non-optimized 
                                      matrix multiplication function
    - CPU_time_fit_optimized.png: plot with fit of the CPU times for the optimized matrix 
                                  multiplication function
    - CPU_time_fit_matmul.png: plot with fit of the CPU times for the intrinsic matrix 
                               multiplication function
    - fit logs:
        + not-optimized_fit.log
        + optimized_fit.log 
        + matmul_fit.log 