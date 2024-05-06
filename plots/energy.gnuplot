set logscale y
set yrange [1e-1:1e6]
set xrange [0:1]

plot "../output/energy.dat" using 1:2 t "Energy" w l, \
"" using 1:4 t "Potental" w l, \
"" using 1:3 t "Kinetic" w l, \
"" using 1:5 t "Temperature" w l