# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "k=n" plot.gnu
# -------------------------------------------------------------
# where 'foo.txt' is the file you want to plot and plot and n
# is the number of eigenfunctions you want to plot 

set terminal pdf size 8, 5 font "Latin Modern Math, 25"
set output sprintf("%s%d%s", "first_", k, "_eigenfunc.pdf")
set title "First ".(k)." eigenfunctions" font ", 27"
set xlabel "x"
set ylabel "{/Symbol Y}(x)"
set grid
set key outside

plot for [i=2:k+1] 'ef_250_0.4E-02_0.1E+01_1.000_1.000.txt' using 1:i with lines title '{/Symbol Y}_{'.(i-1).'}'
plot for [i=2:k+1] 'ef_t_3_250_0.8E-02_0.2E+02_1.000_1.000.txt' using 1:i with lines title '{/Symbol Y}_{th,'.(i-1).'}' 
# lc 0*i dashtype i