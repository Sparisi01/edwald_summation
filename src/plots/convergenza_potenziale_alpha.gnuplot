set terminal png size 1600,1200 enhanced font 'FiraCode,30' lw 2
set output './png/convergenza_range.png'
set datafile separator ";"
set title font ',36'

set monochrome
set key spacing 2 center right outside reverse Left title "αL"
set xrange [0.64:15]
set format y "10^{%T}"

set xtics 0,2
set mxtics 2

set logscale y
set logscale x

set title "Convergence on K-range"
set ylabel 'σ_{rel}' font ',36'
set xlabel 'Range' font ',36'
set grid
plot "../data/convergenza_range/range_variabile_coulomb.csv" using 1:3 w lp pt 5 ps 2.5 t "0",\
"../data/convergenza_range/range_variabile_edw_1.csv" using 1:3 w lp pt 12 ps 2.5 t "1",\
"../data/convergenza_range/range_variabile_edw_2.csv" using 1:3 w lp pt 2 ps 2.5 t "2",\
"../data/convergenza_range/range_variabile_edw_3_5.csv" using 1:3 w lp pt 6 ps 2.5 t "3.5",\
"../data/convergenza_range/range_variabile_edw_5.csv" using 1:3 w lp pt 10 ps 2.5 lw 1 t "5",\
"../data/convergenza_range/range_variabile_edw_7.csv" using 1:3 w lp pt 8 ps 2.5 lw 1 t "7",\
"../data/convergenza_range/range_variabile_edw_9.csv" using 1:3 w lp pt 4 ps 2.5 lw 1 t "9"