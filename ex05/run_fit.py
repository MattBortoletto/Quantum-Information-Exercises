import os
import subprocess

logFiles = ["herm_spacings.log", "herm_loc_aver.log", 
            "diag_spacings.log", "diag_loc_aver.log"]

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
command_herm_loc = ["gnuplot", "-e", "filename='herm_loc_aver.txt'", 
                    "-e", "name='herm_loc'", hermFitFile]
command_diag_glo = ["gnuplot", "-e", "filename='diag_spacings.txt'", 
                    "-e", "name='diag_glo'", diagFitFile]
command_diag_loc = ["gnuplot", "-e", "filename='diag_loc_aver.txt'",
                    "-e", "name='diag_loc'", diagFitFile]
commands = [command_herm_glo, command_herm_loc, command_diag_glo, command_diag_loc]
# run the commands
for i in commands:
    subprocess.Popen(i)
