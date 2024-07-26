set term pngcairo size 700,600 enhanced
set output 'cp.png'
set size ratio 0.8
set xrange []
set yrange [0:20]

set title "Re = 1" font ",16"
set xlabel "time" font ",15"
set ylabel "Cd" font ",15"
set tics font ",12"
set pointsize 1.5

set key top center

a = 2
b = 10
lwide = 2

plot 'c11v3/log' u a:($3==b?$8:NaN) w l lw lwide t "Single Projection", \
        'c11dpv1/log' u a:($3==b?$8:NaN) w l lw lwide t "Double Projection v1", \
                'c11dpv2/log' u a:($3==b?$8:NaN) w l lw lwide t "Double Projection v2", \
                         'ce11/log' u a:($3==b?$8:NaN) w l lw lwide t "embed"                   


