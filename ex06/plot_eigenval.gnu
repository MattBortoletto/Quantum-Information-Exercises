set terminal pdf size 7, 5 font "Latin Modern Math, 25"
#set output sprintf("%s%s", name, ".pdf")
set output "eigenvalues.pdf"
set title "Eigenvalues" font ", 27"
set xlabel "n"
set ylabel "E_n"
set grid
unset key 

plot "en_250_0.8E-02_0.2E+02_1.000_1.000.txt" 
#plot "en_t_500_0.1E-03_0.1E+05_1.000_1.000.txt" 
