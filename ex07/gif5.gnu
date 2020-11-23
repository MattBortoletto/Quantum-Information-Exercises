set encoding utf8
set term gif animate font "Latin Modern Math, 25" size 960,690 delay 1 
set output "evolution_5.gif"

#setting the number of frames
tmax = 500

idx=1
prob_filename = 'prob_time_evol_5.00000.txt'
pot_filename = 'V_time_evol_5.00000.txt'

set xrange [-5:5]
set grid
set yrange[-0.1-(0.5*(2*idx-1)):1]
set xlabel 'x' 
set ylabel 'P(x)'

do for [time = 0:tmax] {
    set title sprintf("t=%i", time)
    plot prob_filename u 1:(column(time+2)) with lines title 'P(x)', \
        pot_filename u 1:(column(time+2))-(0.5*(2*idx-1)) with lines title 'V(x)'
}
