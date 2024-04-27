
set xrange [-1.5:1.5]
set yrange [-1.5:1.5]

set xlabel 'x'
set ylabel 'y'

plot "../output/end_pos.dat" using 1:2 w p pt 7 t "End position", "../output/start_pos.dat" using 1:2 w p pt 7 t "Start position"