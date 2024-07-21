set xrange [0:0.16]
set yrange [-2.5:2.5]

tt = 0.15625
U0 = 1
offset = 0

set title "VOF = 1" font ",16"
set xlabel "u.x" font ",15"
set ylabel "y" font ",15"
set ytics 0.5 font ",12"
set xtics 0.025 font ",12"
set grid

f(x) = sqrt(1 - (x/tt))*(5 + 2 * offset)/2

set key right top

 plot 'vprof-1' u 3:2 t "IBM lvl 7", \
        'vprofref-1' u 3:2 t "Analytical" w l lc rgb "blue", \

