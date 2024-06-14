set terminal png size 1600,1200 enhanced font 'FiraCode,30' lw 2
set output './png/errore_real_space.png'
set datafile separator ";"
set title font ',36'

set monochrome
set key spacing 2 reverse Left

set xrange [0:7]

set logscale y

set format y "10^{%T}"

set title "Real space relative-error"
set ylabel 'Relative error' font ',36'
set xlabel 'αL' font ',36'
set grid
a = -1.2766682221E+00
f(x) = exp(-(x*x))/100
plot '../data/comparison_real_rec.csv' using 1:(abs(($5 - a)/a)) w lp pt 4 ps 2.5 t "Real"




