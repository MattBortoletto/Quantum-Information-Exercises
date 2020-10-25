import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess


#asking the user to insert the dimensions
N_min = 0
while (N_min<=0):
    # catching eventual non-integer values.
    try:
        N_min = int(input("Please enter the minimum size: "))
    except ValueError:
        print("Invalid value: integer needed.")
    # if the integer is non valid
    if (N_min<=0):
        print("Invalid dimension: less or equal to 0.")

N_max = 0
while (N_max<=0 or N_max<N_min):
    # catching eventual non-integer values.
    try:
        N_max = int(input("Please enter the maximum size: "))
    except ValueError:
        print("Invalid value: integer needed.")
    # if the integer is non valid
    if (N_max<=0 or N_max<N_min):
        print("Invalid dimension: N_max should be positive and greater than N_min=",N_min)

# the number of steps from N_min to N_man
num_intermediate = 30

# creating the sizes; logarithmic
dims = np.logspace(np.log10(N_min), np.log10(N_max), num_intermediate, dtype=np.int16)

print(dims)





# iterations = 5
# n_runs = [i for i in range(N_min, N_max)]
# print(n_runs)

# fileSource = "Ex04-MatteoBortoletto.f90"
# fileExec = "Ex04.x"

# subprocess.call(["gfortran", "-o", fileExec, fileSource, "-O3"]) 
# subprocess.call(["./", fileExec])                          
