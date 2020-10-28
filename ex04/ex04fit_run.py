import os
import subprocess

logFiles = ["not-optimized_fit.log", "optimized_fit.log", "matmul_fit.log"]

# if the log files already exist, delete them
for logfile in logFiles:
    if os.path.exists(logfile):
        os.remove(logfile)

# gnuplot script
plotFile = "ex04fit.gnu"
# write the command for each file
command_not_opt = ["gnuplot", "-e", "filename='not-optimized.txt'", "-e", "name='not-optimized'", plotFile]
command_opt = ["gnuplot", "-e", "filename='optimized.txt'", "-e", "name='optimized'", plotFile]
command_matmul = ["gnuplot", "-e", "filename='matmul.txt'", "-e", "name='matmul'", plotFile]
commands = [command_not_opt, command_opt, command_matmul]
# run the commands
for i in commands:
    subprocess.Popen(i)

