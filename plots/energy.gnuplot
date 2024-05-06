set multiplot layout 2,1

set xrange [0:6]

plot "../output/energy.dat" using 1:2 t "Energy" w l, \
"" using 1:4 t "Potental" w l, \
"" using 1:3 t "Kinetic" w l

plot "" using 1:5 t "Temperature" w l