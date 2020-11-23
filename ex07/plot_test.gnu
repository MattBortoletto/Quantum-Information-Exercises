# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" -e "T='foo'" plot_test.gnu
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

tmax = 9000
step = 10

set pm3d
#set palette defined (0 '#052a36',  1 '#73cdeb')
#set palette functions sqrt(gray), gray**3, sin(gray*2*pi) 
set cbrange [0:tmax]

set xrange [-3:4]


plot for [i=1:tmax:step] filename u 1:(column(i+1)) with lines notitle palette cb (i)