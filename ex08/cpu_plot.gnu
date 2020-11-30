set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output "cpu_times.pdf"

set xlabel "N"
set ylabel "CPU time"
set grid
set key top left
set logscale 
set format y "10^{%T}"

plot "cpu_times.txt" u 1:2 w lp title "separable", "cpu_times.txt" u 1:3 w lp title "general"