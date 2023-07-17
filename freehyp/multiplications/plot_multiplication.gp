set term qt persist size 2560,1294 font ',20' lw 2
set xrange [5.8:7.5]
set yrange [-20:30]
set xtics 0.05
set mxtics 5
set xzeroaxis
set mytics 5
set grid
set grid mxtics
set grid mytics
set xlabel "R / Angstrom"
set ylabel "Re or Im (see key) of psi_hypfin * psi_n" noenhanced
set key font "Monospace" noenhanced
FILE = ARG1
plot FILE u 1:($2*$4) w l title FILE
