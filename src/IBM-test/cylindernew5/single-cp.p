set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "cp-100.pdf"
set size ratio 0.8

a = 1
b = 2

set xlabel "Theta (degrees)"
set ylabel "Cp" rotate by 90

set key top center

set style line 1 lt 1 lc rgb "black" lw 3
set style line 2 lt 1 lc rgb "black" lw 3 dt 2

plot '../cylindernew4/surface-100new2.txt' u (abs(abs($1)-180)):2 t "Present" w l ls 1, \
        '../cylindernew4/uhlmann_cp.dat' u 1:2 t "Uhlmann 2005" w l ls 2

set output
