set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set encoding utf8

#set title "Ground state density for different N"
set xlabel "λ" 
set ylabel "e"
set grid
#set key outside 

f(x) = 0 <= x && x <= 2 ? -1 - x**2 / 4 : -x               

set output "gs_N2.pdf"
plot "gs_N2.txt" u 1:2 w l title "RSRG", \
     f(x) title "MF" dt 2

set output "gs_N3.pdf"
plot "gs_N3.txt" u 1:2 w l title "RSRG", \
     f(x) title "MF" dt 2

set output "gs_N4.pdf"
plot "gs_N4.txt" u 1:2 w l title "RSRG", \
     f(x) title "MF" dt 2

set output "gs_N5.pdf"
plot "gs_N5.txt" u 1:2 w l title "RSRG", \
     f(x) title "MF" dt 2

g(x, y) = abs(f(x) - y)
set output "gs_diff.pdf"
set key at graph 0.25,0.98
set ylabel "|e_{MF} - e_{RSRG}|"
plot "gs_N2.txt" u 1:(g($1, $2)) w lp title "N=2", \
     "gs_N3.txt" u 1:(g($1, $2)) w lp title "N=3", \
     "gs_N4.txt" u 1:(g($1, $2)) w lp title "N=4", \
     "gs_N5.txt" u 1:(g($1, $2)) w lp title "N=5"