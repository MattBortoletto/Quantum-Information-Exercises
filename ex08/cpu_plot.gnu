set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output "cpu_times.pdf"

set xlabel "log(N)"
set ylabel "log(CPU time)"
set grid
set key top left
set format y "10^{%T}"
set logscale

plot "cpu_times.txt" u 1:2 w lp title "separable", "cpu_times.txt" u 1:3 w lp title "general"