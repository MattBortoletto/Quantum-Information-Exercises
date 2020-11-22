import os
import subprocess

print("Running...")

# gnuplot script
plotFile = "plot_test.gnu"

# write the command for each file
command_1 = ["gnuplot", "-e", "filename='prob_time_evol_0.00010.txt'", 
                "-e", "T='0.0001'", plotFile]
command_2 = ["gnuplot", "-e", "filename='prob_time_evol_0.00050.txt'", 
                "-e", "T='0.0005'", plotFile]
command_3 = ["gnuplot", "-e", "filename='prob_time_evol_0.00100.txt'", 
                "-e", "T='0.001'", plotFile]
command_4 = ["gnuplot", "-e", "filename='prob_time_evol_0.00500.txt'", 
                "-e", "T='0.005'", plotFile]

commands = [command_1, command_2, command_3, command_4] 

# run the commands
for i in commands:
    subprocess.Popen(i)

print("End.")

#gnuplot -e "filename='prob_time_evol_0.020.txt'" -e "T='0.020'" plot_test.gnu