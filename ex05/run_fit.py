import os
import subprocess

logFiles = ["herm_glo_spacings.log", "herm_loc_spacings.log", 
            "diag_glo_spacings.log", "diag_loc_spacings.log"]

# if the log files already exist, delete them
for logfile in logFiles:
    if os.path.exists(logfile):
        os.remove(logfile)

# gnuplot script
plotFile = "fit.gnu"
# write the command for each file
command_herm_glo = ["gnuplot", "-e", "filename='herm_glo_spacings.txt'", 
                    "-e", "name='herm_glo'", plotFile]
command_herm_loc = ["gnuplot", "-e", "filename='herm_loc_spacings.txt'", 
                    "-e", "name='herm_loc'", plotFile]
command_diag_glo = ["gnuplot", "-e", "filename='diag_glo_spacings.txt'", 
                    "-e", "name='diag_glo'", plotFile]
command_diag_loc = ["gnuplot", "-e", "filename='diag_loc_spacings.txt'",
                    "-e", "name='diag_loc'", plotFile]
commands = [command_herm_glo, command_herm_loc, command_diag_glo, command_diag_loc]
# run the commands
for i in commands:
    subprocess.Popen(i)
