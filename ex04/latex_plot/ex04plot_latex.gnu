set terminal epslatex standalone color colortext 12
set output 'CPU_time_plt.tex'
set title "Matrix multiplication - CPU times"
set xlabel "Matrix dimension"
set ylabel "CPU time (s)"
set grid
set logscale y
set key bottom right box
set format y '$%2.0t\times10^{%T}$'

array methods_files[3]
methods_files[1] = "../not-optimized.txt"
methods_files[2] = "../optimized.txt"
methods_files[3] = "../matmul.txt"

array label_list[3]
label_list[1] = "not optimized"
label_list[2] = "optimized"
label_list[3] = "matmul"

plot for [i=1:3] methods_files[i] using 1:2 title label_list[i] with linespoints lw 2
