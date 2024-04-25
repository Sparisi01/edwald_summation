set terminal pdfcairo font ",14" lw 2
set output 'plot.pdf' 

set title "Titolo Grafico" font ",16"

set xlabel 'Posizione x →'
set ylabel 'Velocità x'

set grid

plot '../output/output.dat' using 2:5 with lines