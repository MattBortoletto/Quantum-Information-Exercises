set terminal pdf size 7, 5 font "Latin Modern Math, 25"
#set output sprintf("%s%s", name, ".pdf")
set output "eigenvalues.pdf"
set title "Eigenvalues" font ", 27"
set xlabel "n"
set ylabel "E_n"
set grid

set xrange [0:100]

hbar = 1
omega = 5
f(x) = hbar*omega*(x+1/2)
plot "en_500_0.1E-01_0.5E+01_1.000_1.000.txt" title "numerical solution" w l, f(x) title "theoretical values"
