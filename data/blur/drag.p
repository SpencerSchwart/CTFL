set term pngcairo size 700,600 enhanced
set output 'u0.png'
set size ratio 0.8
set xrange [200:400]
set yrange [0:1.2]

set title "x/H=0" font ",16"
set xlabel "tUg/H" font ",15"
set ylabel "C" font ",15"
set tics font ",12"
set pointsize 1.5

set key top right

a = 1
b = 4  # vof = 4, u.y = 5
lwide = 2

plot '1-data/0-1-xh0' u ($1*160):4 w l lw lwide lt 1 dt 1 lc "red" t 'M=80', \
        '2-data/0-2-xh0' u ($1*160):4 w l lw lwide lt 1 dt 8 lc "green" t 'M=320', \
                '3-data/0-3-xh0' u ($1*160):4 w l lw lwide lt 1 dt 2 lc "blue" t 'M=1280'

#plot '1-data/0-1-xh2' u ($1*160):($5/160) w l lw lwide lt 1 dt 1 lc "red" t 'M=80', \
#        '2-data/0-2-xh2' u ($1*160):($5/160) w l lw lwide lt 1 dt 8 lc "green" t 'M=320', \
#                '3-data/0-3-xh2' u ($1*160):($5/160) w l lw lwide lt 1 dt 2 lc "blue" t 'M=1280'
