gnuplot ex04plot_all_latex.gnu
latex CPU_time_fit.tex
dvips -o CPU_time_fit.ps CPU_time_fit.dvi

rm CPU_time_fit.aux
rm CPU_time_fit.dvi
rm CPU_time_fit-inc.eps