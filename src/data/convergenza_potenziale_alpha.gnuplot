set terminal png size 1600,1200 enhanced font 'FiraCode,32' lw 4
set output '../plots/convergenza_range.png'
set datafile separator ";"
set title font ',36'
set monochrome
set key spacing 2
set key bottom left

set format y "10^{%T}"

set logscale y

set title "Errore relativo in funzione del range"
set ylabel 'σ rel'
set xlabel 'Range'

plot "range_variabile_coulomb.csv" using 1:3 w lp pt 4 ps 4 t "Coulomb",\
"range_variabile_edw_3_5.csv" using 1:3 w lp pt 6 ps 4 t "α = 3.5",\
"range_variabile_edw_7.csv" using 1:3 w lp pt 8 ps 4 t "α = 7"
