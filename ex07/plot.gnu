# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" -e "T='foo'" plot_test.gnu
# -------------------------------------------------------------
# where 'foo.txt' is the file you want to fit and plot and 'foo' 
# is T.

set encoding utf8
set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output sprintf("%s%s%s", "prob_", T, ".pdf")

set pm3d

set xlabel "x"
set ylabel "|{/Symbol Y}(x)|^2"
set cblabel 't'
set grid

stats filename using 2 nooutput

tmax = 90
step = 5

#set palette defined (0 '#052a36',  1 '#73cdeb')
#set palette functions sqrt(gray), gray**3, sin(gray*2*pi) 

set cbrange [0:tmax]
set xrange [-3:4]

plot for [i=1:tmax:step] filename u 1:(column(i+1)) with lines notitle palette cb (i)