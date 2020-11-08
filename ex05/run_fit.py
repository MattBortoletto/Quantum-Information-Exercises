import os
import subprocess

N = [100, 50, 10, 5, 1]

logFiles = ["herm_spacings.log", "herm_loc_aver_100.log", "herm_loc_aver_50.log", 
            "herm_loc_aver_10.log", "herm_loc_aver_5.log", "herm_loc_aver_1.log",
            "diag_spacings.log", "diag_loc_aver_100.log", "diag_loc_aver_50.log",
            "diag_loc_aver_10.log", "diag_loc_aver_5.log", "diag_loc_aver_1.log",]

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

command_herm_loc = [0]*5
for i in range(0, 4):
    print(N[i])
    command_herm_loc[i] = ["gnuplot", "-e", "filename='herm_loc_aver_"+str(N[i])+".txt'", 
                           "-e", "name='herm_loc"+str(N[i])+"'", hermFitFile]

command_diag_glo = ["gnuplot", "-e", "filename='diag_spacings.txt'", 
                    "-e", "name='diag_glo'", diagFitFile]

command_diag_loc = [0]*5
for i in range(0, 4):
    command_diag_loc[i] = ["gnuplot", "-e", "filename='diag_loc_aver_"+str(N[i])+".txt'",
                           "-e", "name='diag_loc"+str(N[i])+"'", diagFitFile] 

commands = [command_herm_glo, 
            command_herm_loc[0], 
            command_herm_loc[1], 
            command_herm_loc[2], 
            command_herm_loc[3], 
            command_herm_loc[4],
            command_diag_glo, 
            command_diag_loc[0], 
            command_diag_loc[1], 
            command_diag_loc[2], 
            command_diag_loc[3], command_diag_loc[4]]

# run the commands
for i in commands:
    subprocess.Popen(i)
