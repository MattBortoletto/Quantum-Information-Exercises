import os
import subprocess

plotFile = "plot_k_eig.gnu"

command2 = ["gnuplot", "-e", "filename='eig_2.txt'", "-e", "N=2", "-e", "k=4", plotFile]
command3 = ["gnuplot", "-e", "filename='eig_3.txt'", "-e", "N=3", "-e", "k=4", plotFile]
command4 = ["gnuplot", "-e", "filename='eig_4.txt'", "-e", "N=4", "-e", "k=4", plotFile]
command5 = ["gnuplot", "-e", "filename='eig_5.txt'", "-e", "N=5", "-e", "k=4", plotFile]
command6 = ["gnuplot", "-e", "filename='eig_6.txt'", "-e", "N=6", "-e", "k=4", plotFile]
command7 = ["gnuplot", "-e", "filename='eig_7.txt'", "-e", "N=7", "-e", "k=4", plotFile]
command8 = ["gnuplot", "-e", "filename='eig_8.txt'", "-e", "N=8", "-e", "k=4", plotFile]
command9 = ["gnuplot", "-e", "filename='eig_9.txt'", "-e", "N=9", "-e", "k=4", plotFile]
command10 = ["gnuplot", "-e", "filename='eig_10.txt'", "-e", "N=10", "-e", "k=4", plotFile]

commands = [command2, command3, command4, command5, command6, 
            command7, command8, command9, command10]

for i in commands:
    subprocess.Popen(i)