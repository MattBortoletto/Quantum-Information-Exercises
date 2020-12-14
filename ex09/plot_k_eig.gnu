# this gnuplot script must be run using the following command:
# ------------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" -e "N=n" -e "k=k" plot_k_eig.gnu
# ------------------------------------------------------------------
# where 'foo.txt' is the file you want to plot and plot, n is the
# number of subsystems and k is the number of eigenvalues you want 
# to plot 

set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output "eig_N".N.".pdf"

set encoding utf8

set xlabel "λ" 
set ylabel "e"
set grid

set key left bottom 

plot for [i=2:k+1] filename using 1:i with lines title 'e_{'.(i-2).'}'
