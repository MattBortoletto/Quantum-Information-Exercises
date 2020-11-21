set terminal pdf size 8, 5 font "Latin Modern Math, 25"
set output 'prob.pdf'

set title "Probability density time evolution" font ", 27"
set xlabel "x"
set ylabel "|{/Symbol Y}(x)|^2"
set grid
set autoscale xy
set key outside 

stats 'prob_time_evol_0.020.txt' using 2 nooutput

plot for [i=1:300:30] 'prob_time_evol_0.020.txt' u 1:(column(i+1)) with lines title "t=".i