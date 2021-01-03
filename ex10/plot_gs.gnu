set terminal pdf size 9, 5 font "Latin Modern Math, 25"
set encoding utf8

#set title "Ground state density for different N"
set xlabel "λ" 
set ylabel "e"
set grid
#set key outside 

f(x) = 0 <= x && x <= 2 ? -1 - x**2 / 4 : -x               

set output "gs_N2.pdf"
plot "gs_N2.txt" u 1:2 w l title "N=2", \
     f(x) title "Mean field" dt 2

set output "gs_N3.pdf"
plot "gs_N2.txt" u 1:2 w l title "N=2", \
     f(x) title "Mean field" dt 2

set output "gs_N4.pdf"
plot "gs_N2.txt" u 1:2 w l title "N=2", \
     f(x) title "Mean field" dt 2
