set terminal png size 1600,1200 enhanced font 'FiraCode,30' lw 2
set output 'rapporto_real_rec.png'
set datafile separator ";"
set title font ',36'


set monochrome
set key spacing 2 left reverse Left
set xrange [0:10]
set yrange [-23:60]



set xtics 0,2
set mxtics 2

set title "Contribute of Real-Reciprocal-Self therm"
set ylabel 'Energy' font ',36'
set xlabel 'αL' font ',36'
set grid

plot '../data/comparison_real_rec2.csv' using 1:($2) w lp pt 4 ps 2.5 t "Real",\
'../data/comparison_real_rec2.csv' using 1:($3) w lp pt 6 ps 2.5 t "Reciprocal",\
'../data/comparison_real_rec2.csv' using 1:($4) w lp pt 8 ps 2.5 t "Self",\
'../data/comparison_real_rec2.csv' using 1:($2+$3-$4) w lp pt 12 ps 2.5 t "Total"


