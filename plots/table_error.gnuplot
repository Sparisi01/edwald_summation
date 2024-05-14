
set terminal png size 1600,1200 enhanced font ',28' lw 2
set output './png/table_error.png'

set title "Relative error on long distance force due to trilinear approximation" 
set title font ",32"

set key horizontal spacing 3
set key bottom
#set key box lt -1 lw 1

set yrange [-0.04:0.015]

set ylabel 'Relative error'
set xlabel 'Magniture of rᵢⱼ'

# NO CORRECTION:
# Relative error x: -2.670E-03 ± 3.024E-03
# Relative error y: -2.643E-03 ± 3.001E-03
# Relative error z: -2.619E-03 ± 3.079E-0

# WITH CORRECTION:
# Relative error x: -4.148E-06 ± 3.032E-03
# Relative error y : -6.957E-06 ± 3.008E-03
# Relative error z : -1.372E-05 ± 3.087E-03

plot "../output/table_error_file.dat" using 1:2 pt 7 ps 1 lc 1 t "x component",\
"" using 1:3 pt 7 ps 1 lc 2 t "y component",\
"" using 1:4 pt 7 ps 1 lc 4 t "z component"
