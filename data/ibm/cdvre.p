set xrange [0:50]
set yrange [0:15]

set xlabel "Re" font ",15"
set ylabel "Cd" font ",15"
set xtics 10 font ",12"
set ytics 5 font ",12"
set pointsize 1.5



plot 'experimental-cd.txt' u 1:2 t "Experiment" pt 2 lc rgb "black", \
       'cdvre.txt' u 1:2 t "IBM" pt 5 lc rgb "blue"
