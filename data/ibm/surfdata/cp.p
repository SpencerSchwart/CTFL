set term pngcairo size 700,600 enhanced
set output 'cp.png'
set size ratio 0.8
set xrange []
set yrange []

set title "Re = 40" font ",16"
set xlabel "Theta (degrees)" font ",15"
set ylabel "U.y" font ",15"
set tics font ",12"
set pointsize 1.5

set key top center
unset key

a = 1
b = 5
wide = 2

set style line 1 lt 1 lc rgb "black" lw 1 dt 4
set style line 2 lt 1 lc rgb "blue" lw wide dt 5
set style line 3 lt 1 lc rgb "red" lw wide
set style line 4 lt 1 lc rgb "green" lw wide
set style line 5 lt 1 lc rgb "orange" lw wide


plot 'embed/surfaceVelocity1-40' u (abs($3)):b t 'embed' ls 1, \
            'ibm/surfaceVelocity1-40' u (abs($3)):b t 'IBM w/ z' ls 3, \
                'ibmnz/surfaceVelocity1-40' u (abs($3)):b t 'IBM w/o z' ls 5

