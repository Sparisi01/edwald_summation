
set title "Relative error on long distance force due to trilinear approximation"

set key horizontal spacing 3
set key bottom
#set key box lt -1 lw 1

set yrange [-0.02:0.01]

set ylabel 'Relative error'
set xlabel 'Magniture of rᵢⱼ'

# Relative error x: -2.670E-03 ± 3.024E-03
# Relative error y: -2.643E-03 ± 3.001E-03
# Relative error z: -2.619E-03 ± 3.079E-03

plot "../output/table_error_file.dat" using 1:2 pt 7 ps 0.5 t "x component",\
"" using 1:3 pt 7 ps 0.5 t "y component",\
"" using 1:4 pt 7 ps 0.5 t "z component"
