set terminal png font "Latin Modern Math, "
#set encoding utf8
set output 'prob.png'

stats 'prob_time_evol.txt' using 2 nooutput

set xlabel 'x' 
set ylabel '|{/Symbol Y}(x)|^2' 

plot for [i=1:300:10] 'prob_time_evol.txt' u 1:(column(i+1)) with lines notitle