set terminal png size 1600,1200 enhanced font ',30' lw 2
set output './png/tempo_esecuzione.png'
set datafile separator ";"
set title font ',36'

set yrange [0:40000]
set xrange [500:15000]

#set format y "%.0t×10^{%S}"
set monochrome

set key spacing 2 right reverse Left box

# set key spacing 2 reverse Left

set xtics 0,2000
set mxtics 4

set mytics 5

set title "CPU time as a function of N for optimized and non-optimized Ewald"
set ylabel 'CPU time [ms]' font ',36'
set xlabel 'N particles' font ',36'
set grid


# Final set of parameters            Asymptotic Standard Error
# =======================            ==========================
# a               = -1.60379         +/- 58.34        (3638%)
# b               = 0.0226072        +/- 0.003164     (14%)
# c               = 1.5872           +/- 0.03921      (2.471%)

# Final set of parameters            Asymptotic Standard Error
# =======================            ==========================
# a               = -85.8332         +/- 111.3        (129.7%)
# b               = 0.0368587        +/- 0.008153     (22.12%)
# c               = 1.49683          +/- 0.05362      (3.582%)

f(x) = a+(b*x)**c
a = -85.8332
b = 0.0368587
c = 1.49683

# =======================            ==========================
# a               = 298.972          +/- 81.98        (27.42%)
# b               = 0.0399024        +/- 0.0008585    (2.152%)
# c               = 1.97399          +/- 0.007466     (0.3782%)

g(x) = d+(h*x)**l

d = 298.972 # 82
h = 0.0399024 # 0.00011
l = 1.974 # 0.007

# Final set of parameters            Asymptotic Standard Error
# =======================            ==========================
# a               = 89.2679          +/- 87.11        (97.59%)
# b               = 0.0147261        +/- 0.001678     (11.4%)
# c               = 1.8863           +/- 0.04271      (2.264%)

# =======================            ==========================
# a               = -83.7766         +/- 45.54        (54.36%)
# b               = 0.0387571        +/- 0.007685     (19.83%)
# c               = 1.53661          +/- 0.05661      (3.684%)

q(x) = m+(n*x)**o

m = 89.26 # 82
n = 0.0147261 # 0.00011
o = 1.8863 # 0.007

plot '../data/execution_times_file.csv' every 1 using 1:4 w p pt 4 ps 1.5 t "Optimized Verlet", f(x) notitle, '../data/Save_fit_corretto/execution_times_file.csv' every 2 using 1:4 w p pt 6 ps 1.5 t "Optimized No-Verlet", q(x) notitle, '../data/Save_fit_corretto/execution_times_file_square.csv' every 2 using 1:4 w p lw 1 pt 8 ps 1.5 t "Not Optimized", g(x) notitle lw 1




