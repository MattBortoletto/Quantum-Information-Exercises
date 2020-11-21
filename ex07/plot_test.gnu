# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" -e "T='foo'" plot.gnu
# -------------------------------------------------------------
# where 'foo.txt' is the file you want to fit and plot and 'foo' 
# is T.

set terminal png font "Latin Modern Math, "
#set encoding utf8
set output sprintf("%s%s%s", "prob_", T, ".png")

stats filename using 2 nooutput

set xlabel 'x' 
set ylabel 'P(x)'
set cblabel 't'

set pm3d
set palette defined (0 '#052a36',  1 '#73cdeb')
set cbrange [0:100]

plot for [i=1:100:5] filename u 1:(column(i+1)) with lines notitle palette cb (i)