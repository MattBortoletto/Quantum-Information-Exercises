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
    # if N_max <= 0:
    #     print("The input must be an integer.")
    if N_max < N_min:
        print("N_max must be greater than N_min.")

iterations = 20
dimensions = np.logspace(np.log10(N_min), np.log10(N_max), iterations, dtype=np.int16)
print('dimensions:', dimensions)

fileSource = "Ex04-MatteoBortoletto.f90"
fileExec = "Ex04.x"
dimensionsFile = "matrix_dimensions.txt"
outputFiles = ['not-optimized.txt', 'optimized.txt', 'matmul.txt']

# if there exists already a file with the dimensions of the matrices, delete it
if os.path.exists(dimensionsFile):
    os.remove(dimensionsFile)

# if there exist already files with the results, delete them
for f in outputFiles:
    if os.path.exists(f):
        os.remove(f)

# if there exists already the executable file, delete it 
if os.path.exists(fileExec):
    os.remove(fileExec)

# compile the Fortran code 
subprocess.call(["gfortran", fileSource, "-o", fileExec, "-O3"]) 

# store the dimensions in the 'dimensionsFile' and launch the Fortran program
for d in dimensions:
    with open(dimensionsFile, "w+") as inputfile:
        print("I'm printing", str(d))
        inputfile.write(str(d) + '\n' + str(d) + '\n' + str(d) + '\n' + str(d) + '\n')
    """subprocess.call(["./", fileExec])    # does not work ('Permission denied' error)"""
    subprocess.run("./" + fileExec)