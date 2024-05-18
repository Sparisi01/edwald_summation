set datafile separator ";"
#usa la prima riga come legend
set key autotitle columnhead
unset key
set key

plot "./convergenza_coulomb.csv" using 1:2 w lp pt 5 t "Yukawa λ=L/2"