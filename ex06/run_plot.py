import os
import subprocess

k = 0
while (k <= 0):
    try:
        k = int(input("Please enter the number of eigenfunctions you want to plot: "))
    except ValueError:
        print("Error! The input must be an integer.")

# fileSource = "Ex06-MatteoBortoletto.f90"
# fileExec = "Ex06.x"

# # compile the Fortran code 
# subprocess.call(["gfortran", fileSource, "-o", fileExec, "-O3", "-framework Accelerate"]) 

# if there already exists the plot, delete it
# eigenval_plot = "eigenvalues.pdf"
# eigenvect_plot = "first_"+str(k)+"_eigenfunc.pdf"
# if os.path.exists(eigenval_plot):
#     os.remove(eigenvect_plot)

# execute the gnuplot script
eigenval_plotFile = "plot_eigenval.gnu"
eigenvect_plotFile = "plot_eigenvect.gnu"
prob_plotFile = "plot_prob.gnu"
command_eigenval = ["gnuplot", eigenval_plotFile]
command_eigenvect = ["gnuplot", "-e", "k="+str(k), eigenvect_plotFile]
command_prob = ["gnuplot", "-e", "k="+str(k), prob_plotFile]
commands = [command_eigenval, command_eigenvect, command_prob]
for i in commands:
    subprocess.Popen(i)