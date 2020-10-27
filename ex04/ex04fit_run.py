import os
import subprocess


# if there already exists the plot, delete it
# plotImg = "CPU_time_fit.png"
# if os.path.exists(plotImg):
#     os.remove(plotImg)

# try:
#     # Change the current working directory    
#     os.chdir("/plots")
#     print("Directory changed.")
# except OSError:
#     print("Can't change the current working directory.")

# execute the gnuplot script
plotFile = "../ex04fit.gnu"
# write the command for each file
command_not_opt = ["gnuplot", "-e", "filename='not-optimized.txt'", "-e", "name='not-optimized'", plotFile]
command_opt = ["gnuplot", "-e", "filename='optimized.txt'", "-e", "name='optimized'", plotFile]
command_matmul = ["gnuplot", "-e", "filename='matmul.txt'", "-e", "name='matmul'", plotFile]
commands = [command_not_opt, command_opt, command_matmul]

for i in commands:
    subprocess.Popen(i)

