
set multiplot layout 2,1

plot "../output/error_table_reciprocal_128_lin.dat" using 1:3 t "Valore Tabulato" w l, \
"" using 1:4 t "Valore Esatto" w l

set ylabel 'Eₜ/Eᵣ'
set xlabel '|rᵢⱼ|'

plot "../output/error_table_reciprocal_128_lin.dat" using 1:2 notitle w l