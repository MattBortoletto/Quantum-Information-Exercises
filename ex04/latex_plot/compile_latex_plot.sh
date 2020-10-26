gnuplot ex04plot_latex.gnu
latex CPU_time_plt.tex
dvips -o CPU_time_plt.ps CPU_time_plt.dvi

rm CPU_time_plt.aux
rm CPU_time_plt.dvi
rm CPU_time_plt.log
rm CPU_time_plt-inc.eps