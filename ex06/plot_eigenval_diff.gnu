set terminal pdf size 8, 5 font "Latin Modern Math, 25"
set output "eigenvalues_error.pdf"
set title "Eigenvalues error, {/Symbol w}=5, hbar=m=1" font ", 27"
set xlabel "n"
set ylabel "Error"
set grid
set logscale
set key outside 

set format y "10^{%T}";

plot 'en_err_500_0.1E+00_0.5E+01_1.000_1.000.txt' with lines title 'dx=1e-01', \
     'en_err_500_0.1E-01_0.5E+01_1.000_1.000.txt' with lines title 'dx=1e-02', \
     'en_err_500_0.1E-02_0.5E+01_1.000_1.000.txt' with lines title 'dx=1e-03', \
     'en_err_500_0.1E-03_0.5E+01_1.000_1.000.txt' with lines title 'dx=1e-04', \
     'en_err_500_0.1E-04_0.5E+01_1.000_1.000.txt' with lines title 'dx=1e-05'