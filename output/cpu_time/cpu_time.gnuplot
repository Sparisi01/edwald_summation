
#set datafile separator ";"
#usa la prima riga come legend
set key autotitle columnhead
unset key

set key
set grid

#set logscale x
#set logscale y

set yrange [0:0.12]
set xrange [0:1000]


a = 0.000000115
f(x) = a*(x**(2))
g(x) = a*(x**(3./2))

set title ""

set xlabel 'N particles'
set ylabel 'CPU time per verlet step [s]'

plot "./cpu_time_edwald.dat" using 1:2 w lp pt 5 t "Edwald Tabulated", \
"./cpu_time_edwald_cutoff.dat" using 1:2 w lp pt 7 t "Edwald Cutoff", \
"./cpu_time_edwald_parallelized_5.dat" using 1:2 w lp pt 9 t "Edwald Cutoff", \
f(x), g(x)
#"" using 1:3 w l t "Coulomb range 5"