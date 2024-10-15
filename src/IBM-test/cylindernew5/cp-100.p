set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "cp-100.pdf"
set size ratio 0.8
set yrange [-1.5:1.6]

a = 1
b = 2

set xlabel "{/Symbol q} (degrees)" font ",26"
set ylabel "C_p" font ",26"

set key top center font ",24"

set xtics 45
set mxtics 2
set mytics 2

set style line 1 lt 1 lc rgb "black" lw 3
set style line 2 lt 1 lc rgb "red" lw 3 dt 2
set style line 3 lt 1 lc rgb "blue" lw 3 dt 3

plot '../cylindernew4/surface-100new2.txt' u (abs(abs($1)-180)):2 t "Present" w l ls 1, \
        '../cylinderl11v2/surface-100new.txt' u (abs(abs($1)-180)):2 t "Embedded boundary" w l ls 3, \
            '../cylindernew4/uhlmann_cp.dat' u 1:2 t "Uhlmann 2005" w l ls 2, \
                '../cylindernew4/park_cp_100.dat' u 1:2 t "Park et al. 1998" lc rgb "black" pt 4

set output
