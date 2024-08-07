set terminal png size 1600,1200 enhanced font 'FiraCode,30' lw 2
set output './png/tempo_esecuzione.png'
set datafile separator ";"
set title font ',36'

set xrange [0:15000]

# set monochrome

# set key spacing 2 reverse Left

set title ""
set ylabel 'CPU time [ms]' font ',36'
set xlabel 'N' font ',36'
set grid

f(x) = 0.0115*(x**(3./2)) +100
g(x) = 0.0001*(x**(2)) + 2300

plot '../data/execution_times_file.csv' using 1:2 w lp pt 4 ps 2.5, f(x), g(x)




