import os
import subprocess

print("Running...")

# gnuplot script
plotFile = "plot_test.gnu"

# write the command for each file
command_0020 = ["gnuplot", "-e", "filename='prob_time_evol_0.020.txt'", 
                "-e", "T='0.020'", plotFile]
command_0001 = ["gnuplot", "-e", "filename='prob_time_evol_0.001.txt'", 
                "-e", "T='0.001'", plotFile]
command_0100 = ["gnuplot", "-e", "filename='prob_time_evol_0.100.txt'", 
                "-e", "T='0.100'", plotFile]

commands = [command_0020, command_0001, command_0100] 

# run the commands
for i in commands:
    subprocess.Popen(i)

print("End.")

#gnuplot -e "filename='prob_time_evol_0.020.txt'" -e "T='0.020'" plot_test.gnu