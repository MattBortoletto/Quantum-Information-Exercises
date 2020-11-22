# $ gnuplot -e "filename='pr_mean_time_evol_0.10000.txt'" -e "T='0.1'" plot_mean.gnu

set terminal pdf size 8, 5 font "Latin Modern Math, 25"
set output sprintf("%s%s%s", "mean_", T, ".pdf")

set xlabel "t"
set ylabel "Mean position"
set grid 
set autoscale xy

plot filename with lines notitle