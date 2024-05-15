
set terminal png size 1600,1200 enhanced font ',28' lw 2
set output './png/table_error.png'

set multiplot layout 2,1

set title "Relative error on long distance force due to trilinear approximation" 
set title font ",32"

set key horizontal spacing 3
set key bottom
#set key box lt -1 lw 1

set yrange [-0.015:0.005]

set ylabel 'Relative error'
set xlabel 'Magniture of rᵢⱼ'

plot "../output/table_error_file.dat" using 1:2 pt 7 ps 0.5 lc 1 t "x component",\
"" using 1:3 pt 7 ps 0.5 lc 2 t "y component",\
"" using 1:4 pt 7 ps 0.5 lc 4 t "z component"

unset yrange
plot "../output/table_error_file.dat" using 1:5

