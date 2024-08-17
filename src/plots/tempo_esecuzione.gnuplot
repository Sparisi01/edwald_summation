set terminal png size 1600,1200 enhanced font ',30' lw 2
set output './png/tempo_esecuzione.png'
set datafile separator ";"
set title font ',36'

set yrange [0:80000]
set xrange [1000:11000]

set format y "%.0t×10^{%S}"
set monochrome

set key spacing 2 left reverse Left box

# set key spacing 2 reverse Left

set xtics 0,2000
set mxtics 4


set title "CPU time as a function of N for optimized and non-optimized Ewald"
set ylabel 'CPU time [ms]' font ',36'
set xlabel 'N particles' font ',36'
set grid

f(x) = a+b*(x**c)

a = -62.9995
b = 0.00449944
c = 1.58543

g(x) = d+h*(x**l)

d = 298.972 # 82
h = 0.00173 # 0.00011
l = 1.974 # 0.007

plot '../data/Save_fit_corretto/execution_times_file.csv' using 1:4 w p pt 6 ps 1.5 t "Optimized", f(x) notitle, '../data/Save_fit_corretto/execution_times_file_square.csv' using 1:4 w p pt 8 ps 1.5 t "Not optimized", g(x) notitle




