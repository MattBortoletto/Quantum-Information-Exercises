set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output "gs.pdf"

set encoding utf8

#set title "Ground state density for different N"
set xlabel "λ" 
set ylabel "e"
set grid

plot "eig_2.txt" u 1:2 w l title "N=2", \
     "eig_3.txt" u 1:2 w l title "N=3", \
     "eig_4.txt" u 1:2 w l title "N=4", \
     "eig_5.txt" u 1:2 w l title "N=5", \
