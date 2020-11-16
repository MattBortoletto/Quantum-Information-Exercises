# set terminal pdf size 7, 5 font "Latin Modern Math, 25"
# #set output sprintf("%s%s", name, ".pdf")
# set output "eigenvalues.pdf"
# set title "Eigenvalues" font ", 27"
# set xlabel "n"
# set ylabel "E_n"
# set grid

# set xrange [0:100]

# hbar = 1
# omega = 5
# f(x) = hbar*omega*(x+1/2)
# plot "en_500_0.1E-01_0.5E+01_1.000_1.000.txt" title "numerical solution" w l, f(x) title "theoretical values"

set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output "eigenvalues.pdf"

set multiplot

set title "Eigenvalues" font ", 27"
set xlabel "n"
set ylabel "E_n"
set grid
hbar = 1
omega = 5
f(x) = hbar*omega*(x+1/2)
plot "en_500_0.1E-01_0.5E+01_1.000_1.000.txt" title "numerical solution" w l, f(x) title "theoretical values"


set origin .2, .5
set size .35, .35
clear
unset title 
unset key
unset grid
unset object
unset arrow
unset xlabel 
unset ylabel 
set xtics 20
set ytics 200
set bmargin 1
set tmargin 1
set lmargin 3
set rmargin 1
set xrange [0:100]
set yrange [0:700]
set grid 
plot "en_500_0.1E-01_0.5E+01_1.000_1.000.txt" title "numerical solution" w l, f(x) title "theoretical values"
unset multiplot

