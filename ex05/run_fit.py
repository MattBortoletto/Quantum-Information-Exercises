import os
import subprocess

N = 0
while (N <= 0):
    try:
        N = int(input("Please enter N (read from the txt output files): "))
    except ValueError:
        print("Error! The input must be an integer.")

logFiles = ["herm_spacings.log", "herm_loc_aver_"+str(N)+".log", 
            "diag_spacings.log", "diag_loc_aver_"+str(N)+".log"]

# if the log files already exist, delete them
for logfile in logFiles:
    if os.path.exists(logfile):
        os.remove(logfile)

# gnuplot script
hermFitFile = "herm_fit.gnu"
diagFitFile = "diag_fit.gnu"
# write the command for each file
command_herm_glo = ["gnuplot", "-e", "filename='herm_spacings.txt'", 
                    "-e", "name='herm_glo'", hermFitFile]
command_herm_loc = ["gnuplot", "-e", "filename='herm_loc_aver_"+str(N)+".txt'", 
                    "-e", "name='herm_loc'", hermFitFile]
command_diag_glo = ["gnuplot", "-e", "filename='diag_spacings.txt'", 
                    "-e", "name='diag_glo'", diagFitFile]
command_diag_loc = ["gnuplot", "-e", "filename='diag_loc_aver_"+str(N)+".txt'",
                    "-e", "name='diag_loc'", diagFitFile]
commands = [command_herm_glo, command_herm_loc, command_diag_glo, command_diag_loc]
# run the commands
for i in commands:
    subprocess.Popen(i)
