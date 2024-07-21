set term pngcairo size 800,800 enhanced
set output 'cl-v-time.png'
set size ratio 0.8
set xrange [0:15]
set yrange [0:]

set title "CL vs Time" font ",16"
set xlabel "Time" font ",15"
set ylabel "CL" font ",15"
set tics font ",12"
set pointsize 1.5

a= 2

plot 'logaf' u a:($3==0?$9:NaN) t 'SPM' w l lc "green", \
   'logaf' u a:($3==1?$9:NaN) t 'no SPM' w l lc "blue", \
    'logafe' u a:9 t 'embed' w l lc "red"





