set terminal png size 1024, 768 font "Verdana, 18"
set output "CPU_time_plt.png"
set title "Matrix multiplication - CPU times" font "Verdana, 20"
set xlabel "Matrix dimension"
set ylabel "CPU time (s)"
set grid
set logscale y
set key bottom right box height 1.5

array methods_files[3]
methods_files[1] = "not-optimized.txt"
methods_files[2] = "optimized.txt"
methods_files[3] = "matmul.txt"

array label_list[3]
label_list[1] = "not optimized"
label_list[2] = "optimized"
label_list[3] = "matmul"

plot for [i=1:3] methods_files[i] using 1:2 title label_list[i] with linespoints lw 2
