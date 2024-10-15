set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "cp.pdf"
set size ratio 0.8

a = 1
b = 2

set xlabel "Theta (degrees)"
set ylabel "Cp"

set key top center
unset key

set arrow from 80,1 to 40,-0.5 lw 2 lc rgb "black" head filled size screen 0.03,15,45
set label "Re" at 38,-0.65 font "Arial,22"

set arrow from 160,-1.9 to 155,0 lw 2 lc rgb "black" head filled size screen 0.03,15,45
set label "Re" at 150,0.2 font "Arial,22"

set style line 1 lt 1 lc rgb "black" lw 3

plot 'surface-2new.txt' u (abs(abs($1)-180)):2 t "Re = 2" w l ls 1, \
       'surface-5new.txt' u (abs($1-180)):2 t "Re = 5" w l ls 1, \
         'surface-10new.txt' u (abs($1-180)):2 t "Re = 10" w l ls 1, \
           '../cylindernew4/surface-20new.txt' u (abs($1-180)):2 t "Re = 20" w l ls 1, \
             '../cylindernew4/surface-40new.txt' u (abs($1-180)):2 t "Re = 40" w l ls 1, \

set output
