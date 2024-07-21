set term pngcairo size 800,800 enhanced
set output 'cl-v-time.png'
set size ratio 0.8
set xrange [2:4]
set yrange [0.4:0.6]

set title "" font ",16"
set xlabel "y" font ",15"
set ylabel "u.x" font ",15"
set tics font ",12"
set pointsize 1.5

a = 3

plot 'ibmv2/vprofx3-10' u a:4 t 'new IBM' w l lc "green", \
   'c11v3/vprofx3-25' u a:4 t 'old IBM' w l lc "blue", \
    'ce11/vprofx3-15' u a:4 t 'embed' w l lc "red"





