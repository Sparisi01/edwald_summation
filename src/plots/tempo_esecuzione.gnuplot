set terminal png size 1600,1200 enhanced font 'FiraCode,30' lw 2
set output './png/tempo_esecuzione.png'
set datafile separator ";"
set title font ',36'

#set xrange [1000:10000]


# set key spacing 2 reverse Left

set title ""
set ylabel 'CPU time [ms]' font ',36'
set xlabel 'N' font ',36'
set grid

f(x) = a+b*(x**c)

a = -62.9995
b = 0.00449944
c = 1.58543

plot '../data/execution_times_file.csv' using 1:4 w lp pt 2 ps 1.5, f(x), '../data/execution_times_file_square.csv' using 1:4 w lp pt 3 ps 1.5




