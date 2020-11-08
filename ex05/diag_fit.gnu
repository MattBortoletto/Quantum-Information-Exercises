# this gnuplot script must be run using the following command:
# -------------------------------------------------------------
# $ gnuplot -e "filename='foo.txt'" -e "name='foo'" plot.gnu
# -------------------------------------------------------------
# where 'foo.txt' is the file you want to fit and plot and 'name' 
# is an identificative name for the dataset (e.g. 'optimized') 
# which will be used in the output file name. 

set terminal pdf size 7, 5 font "Latin Modern Math, 25"
set output sprintf("%s%s", name, ".pdf")
set title "Spacings distribution - Real diagonal matrix" font "Latin Modern Math, 27"
set xlabel "s"
set ylabel "P(s)"
set grid
set key box height 1 width 1

fileLog = name."_fit.log"
set fit logfile fileLog

p(x) = a * (x**alpha) * exp(-b*(x**beta))

fit p(x) filename via a, b, alpha, beta 

set encoding iso_8859_1
stats filename using 2 nooutput
ylabel_position = STATS_max
stats filename using 1 nooutput
xlabel_position = STATS_max
set label sprintf("a = %1.3f \261 %1.3f \
                   \n{/Symbol a} = %1.3f \261 %1.3f \
                   \nb = %1.3f \261 %1.3f \
                   \n{/Symbol b} = %1.3f \261 %1.3f", \
                  a, a_err, alpha, alpha_err, b, b_err,beta, beta_err) \
                  at 0.68*xlabel_position, 0.7*ylabel_position

plot filename using 1:2 with p title "data", p(x) title "fit" lw 2