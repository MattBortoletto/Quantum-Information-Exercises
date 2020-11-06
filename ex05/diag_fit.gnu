# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" -e "name='foo'" plot.gnu
# -------------------------------------------------------------
# where 'foo.txt' is the file you want to fit and plot and 'name' 
# is an identificative name for the dataset (e.g. 'optimized') 
# which will be used in the output file name. 

set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output sprintf("%s%s", name, ".pdf")
set title "Spacings distribution" font "Latin Modern Math, 27"
set xlabel "s"
set ylabel "P(s)"
set grid
set key font ",18" box height 1 width 1

fileLog = name."_fit.log"
set fit logfile fileLog

a = 1.8
b = 2 
beta = 0.8 
alpha = 0.01
p(x) = a*exp(-b*x**beta)*x**alpha

fit p(x) filename using 1:2 via a, b, beta, alpha  

plot filename using 1:2 with p title "data", p(x) title "fit" lw 2