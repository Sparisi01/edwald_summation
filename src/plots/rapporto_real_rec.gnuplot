set terminal png size 1600,1200 enhanced font ',30' lw 2
set output './png/rapporto_real_rec.png'
set datafile separator ";"
set title font ',36'

set monochrome
set key spacing 2 left reverse Left box

set xrange [1.5:15]
set yrange [-2:9]

set xtics 0,1
set mxtics 2

set ytics -3,1
set mytics 2

set lmargin 7

set title "Energy contributions as functions of alpha"
set ylabel 'Energy per particle [g cm^{2} s^{-2}]' font ',36'
set xlabel 'αL' font ',36'
set grid

plot '../data/comparison_real_rec.csv' using 1:($2) w lp pt 4 ps 2.5 t "Short",\
'../data/comparison_real_rec.csv' using 1:($3) w lp pt 6 ps 2.5 t "Long",\
'../data/comparison_real_rec.csv' using 1:($4) w lp pt 8 ps 2.5 t "Self",\
'../data/comparison_real_rec.csv' using 1:($2+$3-$4) w lp pt 13 ps 2.5 t "Short + Long - Self"



