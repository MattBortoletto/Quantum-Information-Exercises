set terminal epslatex standalone color colortext 12
set output "CPU_time_fit.tex"
set title "Matrix multiplication - CPU times" 
set xlabel "Matrix dimension" 
set ylabel "CPU time (s)" 
set grid

set logscale y
set key bottom right horizontal width -3 box
set format y '$%2.0t\times10^{%T}$'
set yrange [1e-7:100]

#set key top left width -2 box


array methods_files[3]
methods_files[1] = "../not-optimized.txt"
methods_files[2] = "../optimized.txt"
methods_files[3] = "../matmul.txt"

logFile = "CPU_time_fit.log"
set fit logfile logFile

# a = 0.1
# b = 1e-10
# c = 3
# f(x) = a + b*x**c #+ c*x**2 + d*x**3
# fit f(x) methods_files[1] using 1:2 via a, b, c #, d

# i = 0.1
# j = 1e-10
# k = 3
# g(x) = i + j*x**k #+ k*x**2 + l*x**3
# fit g(x) methods_files[2] using 1:2 via i, j, k #, l

# p = 0.1
# q = 1e-10
# r = 3
# h(x) = p + q*x**r #+ r*x**2 + s*x**3
# fit h(x) methods_files[3] using 1:2 via p, q, r # ,s 

a = 1e-9
b = 3
f(x)= a*x**b
fit f(x) methods_files[1] using 1:2 via a, b 

i = 1e-9
j = 3
g(x)= i*x**j
fit g(x) methods_files[2] using 1:2 via i, j

p = 1e-9
q = 3
h(x)= p*x**q
fit h(x) methods_files[3] using 1:2 via p, q 

array label_list[3]
label_list[1] = "not optimized"
label_list[2] = "optimized"
label_list[3] = "matmul"

plot for [i=1:3] methods_files[i] using 1:2 title label_list[i] with linespoints lw 2, f(x) lw 2 title "$f(x) = ax^b$", g(x) lw 2 title "$g(x) = ix^j$", h(x) lw 2 title "$h(x) = px^q$"
