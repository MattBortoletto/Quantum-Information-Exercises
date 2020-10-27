set terminal epslatex standalone color colortext 12
set output "CPU_time_fit.tex"
set title "Matrix multiplication - CPU times" 
set xlabel "Matrix dimension" 
set ylabel "CPU time (s)" 
set grid
set logscale y
set key bottom right horizontal width -3 box
set format y '$%2.0t\times10^{%T}$'

array methods_files[3]
methods_files[1] = "../not-optimized.txt"
methods_files[2] = "../optimized.txt"
methods_files[3] = "../matmul.txt"

logFile = "CPU_time_fit.log"
set fit logfile logFile

f(x) = a + b*x + c*x**2 + d*x**3
fit f(x) methods_files[1] using 1:2 via a, b, c, d

g(x) = i + j*x + k*x**2 + l*x**3
fit g(x) methods_files[2] using 1:2 via i, j, k, l

h(x) = p + q*x + r*x**2 + s*x**3
fit h(x) methods_files[3] using 1:2 via p, q, r ,s 

array label_list[3]
label_list[1] = "not optimized"
label_list[2] = "optimized"
label_list[3] = "matmul"

plot for [i=1:3] methods_files[i] using 1:2 title label_list[i] with linespoints lw 2, f(x) lw 2, g(x) lw 2, h(x) lw 2
