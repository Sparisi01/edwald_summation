set xrange [-1.3:1.3]
set yrange [-1.3:1.3]

set xlabel 'x'
set ylabel 'y'

plot "../output/end_pos.dat" using 1:2 w p pt 7 ps 1 t "End position", \
"../output/start_pos.dat" using 1:2:3 w p pt 7 ps 1 t "Start position"