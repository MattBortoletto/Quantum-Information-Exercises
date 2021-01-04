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
     "eig_2.txt" u 1:2 w l title "Exact", \
     f(x) title "MF" dt 2

set output "gs_N3.pdf"
plot "gs_N3.txt" u 1:2 w l title "RSRG", \
     "eig_3.txt" u 1:2 w l title "Exact", \
     f(x) title "MF" dt 2

set output "gs_N4.pdf"
plot "gs_N4.txt" u 1:2 w l title "RSRG", \
     "eig_4.txt" u 1:2 w l title "Exact", \
     f(x) title "MF" dt 2

set output "gs_N5.pdf"
plot "gs_N5.txt" u 1:2 w l title "RSRG", \
     "eig_5.txt" u 1:2 w l title "Exact", \
     f(x) title "MF" dt 2
