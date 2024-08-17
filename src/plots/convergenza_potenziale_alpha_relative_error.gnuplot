set terminal png size 1600,1200 enhanced font ',30' lw 2
set output './png/convergenza_range.png'
set datafile separator ";"
set title font ',36'

set monochrome
set key spacing 2 center right outside reverse Left title "αL"
set xrange [1:16]
set format y "10^{%T}"

set xtics 1,2
set mxtics 2

set logscale y


set title "Energy convergence as a function of the long range cut-off"
set ylabel 'σ_{rel}' font ',36'
set xlabel 'Long range cut-off in units of 2π/L' font ',36'
set grid
plot "../data/convergenza_range/range_variabile_coulomb.csv" using 1:3 w lp pt 5 ps 2.5 t "No Ewd",\
"../data/convergenza_range/range_variabile_edw_2.csv" using 1:3 w lp pt 2 ps 2.5 t "2",\
"../data/convergenza_range/range_variabile_edw_3_5.csv" using 1:3 w lp pt 6 ps 2.5 t "3.5",\
"../data/convergenza_range/range_variabile_edw_5.csv" using 1:3 w lp pt 10 ps 2.5 lw 1 t "5",\
"../data/convergenza_range/range_variabile_edw_7.csv" using 1:3 w lp pt 8 ps 2.5 lw 1 t "7",\
"../data/convergenza_range/range_variabile_edw_9.csv" using 1:3 w lp pt 4 ps 2.5 lw 1 t "9"