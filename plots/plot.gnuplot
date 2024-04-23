set multiplot layout 3,1
set xlabel 'x'
set ylabel 'y'
plot "../output/output.dat" using 2:3 with lines, "" using 8:9 with lines
set xlabel 'y'
set ylabel 'z'
plot "../output/output.dat" using 3:4 with lines, "" using 9:10 with lines
set xlabel 'x'
set ylabel 'z'
plot "../output/output.dat" using 2:4 with lines, "" using 8:10 with lines