import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess


N_min = 0
while (N_min <= 0):
    try:
        N_min = int(input("Please enter the minimum dimension N_min: "))
    except ValueError:
        print("Error! The input must be an integer.")

N_max = 0
while (N_max <= 0 or N_max < N_min):
    try:
        N_max = int(input("Please enter the maximum dimension N_max: "))
    except ValueError:
        print("Error! The input must be an integer.")
    if N_max < N_min:
        print("Error! N_max must be greater than N_min.")

iterations = 50
##### for generating equally spaced points in log-log. 
##### however it is better to have less points in the beginning
##### so we use a linspace
# dimensions = np.logspace(np.log10(N_min), np.log10(N_max), iterations, dtype=int)
#####
dimensions = np.linspace(N_min, N_max, iterations, dtype=int)

fileSource = "Ex04-MatteoBortoletto.f90"
fileExec = "Ex04.x"
dimensionsFile = "matrix_dimensions.txt"
outputFiles = ['not-optimized.txt', 'optimized.txt', 'matmul.txt']

# if there already exists a file with the dimensions of the matrices, delete it
if os.path.exists(dimensionsFile):
    os.remove(dimensionsFile)

# if there already exist files with the results, delete them
for f in outputFiles:
    if os.path.exists(f):
        os.remove(f)

# if there already exists the executable file, delete it 
if os.path.exists(fileExec):
    os.remove(fileExec)

# compile the Fortran code 
subprocess.call(["gfortran", fileSource, "-o", fileExec, "-O3"]) 

# store the dimensions in the 'dimensionsFile' and launch the Fortran program
for d in dimensions:
    with open(dimensionsFile, "w+") as inputfile:
        print("Computing for N =", str(d))
        inputfile.write(str(d) + '\n') 
    subprocess.run("./" + fileExec)

# if there already exists the plot, delete it
plotImg = "CPU_time_plt.png"
if os.path.exists(plotImg):
    os.remove(plotImg)

# execute the gnuplot script
plotFile = "ex04plot.gnu"
subprocess.call(["gnuplot", plotFile])