set terminal png size 1600,1200 enhanced font ',30' lw 2
set output './png/tempo_esecuzione.png'
set datafile separator ";"
set title font ',36'

set yrange [0:80000]
set xrange [1000:11000]

set format y "%.0t×10^{%S}"
set monochrome

set key spacing 2 left reverse Left box

# set key spacing 2 reverse Left

set xtics 0,2000
set mxtics 4



set title "CPU time as a function of N for optimized and non-optimized Ewald"
set ylabel 'CPU time [ms]' font ',36'
set xlabel 'N particles' font ',36'
set grid

#f(x) = a+b*(x**c)

# Final set of parameters            Asymptotic Standard Error
# =======================            ==========================
# a               = -71.3822         +/- 35.05        (49.1%)
# b               = 0.0347844        +/- 0.005355     (15.4%)
# c               = 1.56976          +/- 0.04609      (2.936%)

f(x) = a+(b*x)**c
a = -71.38
b = 0.03478
c = 1.56976

g(x) = d+h*(x**l)

d = 298.972 # 82
h = 0.00173 # 0.00011
l = 1.974 # 0.007

# Final set of parameters            Asymptotic Standard Error
# =======================            ==========================
# a               = 89.2679          +/- 87.11        (97.59%)
# b               = 0.0147261        +/- 0.001678     (11.4%)
# c               = 1.8863           +/- 0.04271      (2.264%)

q(x) = m+(n*x)**o

m = 89.26 # 82
n = 0.014 # 0.00011
o = 1.8863 # 0.007

plot '../data/Save_fit_corretto/execution_times_file.csv' using 1:4 w p pt 6 ps 1.5 t "Optimized", f(x) notitle, q(x), '../data/Save_fit_corretto/execution_times_file_square.csv' using 1:4 w p pt 8 ps 1.5 t "Not optimized", g(x) notitle lw 1




