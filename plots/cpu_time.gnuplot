
set datafile separator ";"
#usa la prima riga come legend
set key autotitle columnhead
unset key

set key
set grid

#set logscale x
#set logscale y

set yrange [0:0.8]


set title ""

set xlabel 'N particles'
set ylabel 'CPU time per verlet step [s]'

plot "../output/cpu_time_comparison.csv" using 1:2 w l t "Edwald Tabulated", \
"" using 1:3 w l t "Coulomb range 5"