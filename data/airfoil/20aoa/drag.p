set term pngcairo size 700,600 enhanced
set output 'c1-1.png'
set size ratio 0.8
set xrange [0:15]
set yrange [0:2.5]

set title "AOA = 20, RE = 2000" font ",16"
set xlabel "time" font ",15"
set ylabel "Cl" font ",15"
set tics font ",12"
set pointsize 1.5

set key top center

a = 2
b = 9
lwide = 2

plot 'log' u a:($3==0?$9:NaN) w l lw lwide t "sharp IBM", \
        'log' u a:($3==1?$9:NaN) w l lw lwide t "SPM", \
                'embed/log' u a:b w l lw lwide t "embed"

