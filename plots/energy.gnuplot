set yrange [-0.0001:0.0002]

plot "../output/energy.dat" using 1:2 t "Energy" w l, \
"" using 1:4 t "Kinetic" w l, \
"" using 1:3 t "Potential" w l