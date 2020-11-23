# $ gnuplot -e "filename='pr_mean_time_evol_0.10000.txt'" -e "T='0.1'" plot_mean.gnu

set terminal pdf size 7, 5 font "Latin Modern Math, 25"

set xlabel "t"
set ylabel "Mean position"
set grid 
set autoscale xy

omega = 1
TT = 1
dt = 0.01

a = 0.00097
f(x) = a*x 
fit f(x) filename via a 
set output sprintf("%s%s%s", "mean_", T, "_detr.pdf")
h(x) = - sin(x*dt*omega) / (omega*TT)
plot filename using 1:($2-f($1)) title "Detrended series" w lines lt 2 lw 2, h(x) title "Theory" w l lt 1 dt 2 lw 2


set output sprintf("%s%s%s", "mean_", T, ".pdf")
g(x) = x*dt / TT - sin(x*dt*omega) / (omega*TT) 
plot filename with lines lt 2 lw 2 title "Mean position", g(x) title "Theory" w l lt 1 dt 2 lw 2

