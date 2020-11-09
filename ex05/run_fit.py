import os
import subprocess


logFiles = ["herm_loc_aver_100.log", "herm_loc_aver_50.log", "herm_glo.log",  
            "herm_loc_aver_10.log", "herm_loc_aver_5.log", "diag_glo.log", 
            "diag_loc_aver_100.log", "diag_loc_aver_50.log","diag_loc_aver_10.log", 
            "diag_loc_aver_5.log"]

# if the log files already exist, delete them
for logfile in logFiles:
    if os.path.exists(logfile):
        os.remove(logfile)

# gnuplot script
hermFitFile = "herm_fit.gnu"
diagFitFile = "diag_fit.gnu"

# write the command for each file
command_herm_glo = ["gnuplot", "-e", "filename='herm_glo.txt'", 
                    "-e", "name='herm_glo'", hermFitFile]
command_herm_loc_5 = ["gnuplot", "-e", "filename='herm_loc_aver_5.txt'", 
                      "-e", "name='herm_loc_5'", hermFitFile]
command_herm_loc_10 = ["gnuplot", "-e", "filename='herm_loc_aver_10.txt'", 
                       "-e", "name='herm_loc_10'", hermFitFile]
command_herm_loc_50 = ["gnuplot", "-e", "filename='herm_loc_aver_50.txt'", 
                       "-e", "name='herm_loc_50'", hermFitFile]
command_herm_loc_100 = ["gnuplot", "-e", "filename='herm_loc_aver_100.txt'", 
                        "-e", "name='herm_loc_100'", hermFitFile]

command_diag_glo = ["gnuplot", "-e", "filename='diag_glo.txt'", 
                    "-e", "name='diag_glo'", diagFitFile]
command_diag_loc_5 = ["gnuplot", "-e", "filename='diag_loc_aver_5.txt'", 
                      "-e", "name='diag_loc_5'", diagFitFile]
command_diag_loc_10 = ["gnuplot", "-e", "filename='diag_loc_aver_10.txt'", 
                       "-e", "name='diag_loc_10'", diagFitFile]
command_diag_loc_50 = ["gnuplot", "-e", "filename='diag_loc_aver_50.txt'", 
                       "-e", "name='diag_loc_50'", diagFitFile]
command_diag_loc_100 = ["gnuplot", "-e", "filename='diag_loc_aver_100.txt'", 
                        "-e", "name='diag_loc_100'", diagFitFile]

commands = [command_herm_glo, 
            command_herm_loc_5, 
            command_herm_loc_10, 
            command_herm_loc_50, 
            command_herm_loc_100, 
            command_diag_glo, 
            command_diag_loc_5, 
            command_diag_loc_10, 
            command_diag_loc_50, 
            command_diag_loc_100] 

# run the commands
for i in commands:
    subprocess.Popen(i)
