set xrange [-.5:5]
set yrange [-5:5]

# zoom
set yrange [-1:1]
set xrange [-0.1:0.4]

tt = 0.078125
U0 = 5

set title "VOF = 0.5" font ",16"
set xlabel "u.x" font ",15"
set ylabel "y" font ",15"
set tics t font ",12"



f(y) = (y*(5 + tt * 1.75)/U0) - (tt * 1.75)
g(x) = -tt * (1.75)

set key center top

plot 'vprof3' u 3:2 t "IBM" ps 3, \
        f(x) t "Analytical" lw 3, \
          g(x) t "Interface" lw 1 lc rgb "black"
