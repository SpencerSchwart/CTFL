set term pngcairo size 700,600 enhanced
set output 'cl.png'
set size ratio 0.8
set xrange [14:15]
set yrange [-4:4]

set title "Re = 200  KC = 10" font ",16"
set xlabel "t/T" font ",15"
set ylabel "CD,CL" font ",15"
set tics font ",12"
set pointsize 1.5

set key top center

a = 2
b = 9
T = 5 # period
lwide = 2

plot 'log' u ($2/T):($3==1?$9:NaN) w l linecolor rgb "blue" lw lwide t "CD v2", \
        'fig15_cycle14' u ($1+14):($2/2) w l linecolor rgb "red" lw lwide t "paper", \
                'logv1' u ($2/T):($3==1?$9:NaN) w l linecolor rgb "green" lw lwide t "CD v1"
#        'log' u ($2/T):($3==1?$8:NaN) w l linecolor rgb "red" lw lwide t "CL v2"

