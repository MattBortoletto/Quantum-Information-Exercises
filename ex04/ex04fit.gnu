# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" -e "name='foo'" ex04fit.gnu
# -------------------------------------------------------------
# where 'foo.txt' is the file you want to fit and plot and 'name' 
# is an identificative name for the dataset (e.g. 'optimized') 
# which will be used in the output file name. 

set terminal png size 1024, 768 font "Verdana, 18"
set output sprintf("%s%s%s", "CPU_time_fit_", name, ".png")
set title "Matrix multiplication - CPU times" font "Verdana, 20"
set xlabel "Matrix dimension"
set ylabel "CPU time (s)"
set grid
set logscale y

fileLog = name."_fit.log"
set fit logfile fileLog

f(x) = a + b*x + c*x**2 + d*x**3
fit f(x) filename using 1:2 via a, b, c, d

plot f(x) title 'f(x)=a+b*x+c*x**2+d*x**3' lw 2, filename using 1:2 title filename with linespoints lw 2