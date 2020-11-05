# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" -e "name='foo'" ex04fit.gnu
# -------------------------------------------------------------
# where 'foo.txt' is the file you want to fit and plot and 'name' 
# is an identificative name for the dataset (e.g. 'optimized') 
# which will be used in the output file name. 

set terminal png size 1024, 768 font "Verdana, 25"
set output sprintf("%s%s%s", "CPU_time_fit_", name, ".png")
set title "Matrix multiplication - CPU times" font "Verdana, 27"
set xlabel "Matrix dimension"
set ylabel "CPU time (s)"
set grid
set logscale 
set key bottom right box height 1.7
set format y '%2.0t*10^{%T}'

fileLog = name."_fit.log"
set fit logfile fileLog

a = 1e-9
b = 3
f(x)= a*x**b
fit f(x) filename using 1:2 via a, b 

plot f(x) title 'f(x)=ax^b' lw 2, filename using 1:2 title name with points lw 2