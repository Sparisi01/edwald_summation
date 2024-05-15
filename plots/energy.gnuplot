set terminal png size 1600,1200 enhanced font ',28' lw 2
set output './png/energy.png'

set multiplot layout 2,1

set title "Relative error on long distance force due to trilinear approximation" 
set title font ",32"

plot "../output/energy.dat" using 1:2 t "Energy" w l \
# "" using 1:4 t "Potental" w l, \
# "" using 1:3 t "Kinetic" w l

unset title
unset yrange

plot "" using 1:5 t "Temperature" w l