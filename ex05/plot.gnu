# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" plot.gnu
# -------------------------------------------------------------
# where 'foo.txt' is the file you want to fit and plot.

set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output "hist.pdf"
set title "Spacings distribution" font "Latin Modern Math, 27"
set xlabel "s"
set ylabel "P(s)"
set grid
set key font ",18" box height 1 width 1

plot filename using 1:2 with lp lw 2 title "data"