set terminal png size 1600,1200 enhanced font 'FiraCode,30' lw 2
set output './png/rapporto_real_rec.png'
set datafile separator ";"
set title font ',36'

set monochrome
set key spacing 2 left reverse Left
set xrange [0:10]
set yrange [-3:6]

set xtics 0,1
set mxtics 2

set ytics -3,1
set mytics 2

set title "Contribute of Real-Reciprocal-Self energy"
set ylabel 'Energy per particle' font ',36'
set xlabel 'αL' font ',36'
set grid

plot '../data/comparison_real_rec.csv' using 1:($2) w lp pt 4 ps 2.5 t "E_{R}",\
'../data/comparison_real_rec.csv' using 1:($3) w lp pt 6 ps 2.5 t "E_{K}",\
'../data/comparison_real_rec.csv' using 1:($4) w lp pt 8 ps 2.5 t "E_{self}",\
'../data/comparison_real_rec.csv' using 1:($2+$3-$4) w lp pt 12 ps 2.5 t "E_{tot}"



