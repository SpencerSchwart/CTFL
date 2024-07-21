set xrange [-.1:1]
set yrange [-5:5]

# zoom
set yrange [-0.5:1]
set xrange [-0.1:0.2]

tt = 0.15625
U0 = 1
offset = 0.75;



set title "VOF = 0.25" font ",16"
set xlabel "u.x" font ",15"
set ylabel "y" font ",15"
set ytics tt font ",12"
set xtics 1 font ",12"
set grid


f(y) = (y*(5 + tt * offset)/U0) - (tt * offset)
g(x) = -tt * (offset)

set key center top

plot 'vprof-0.25' u 3:2 t "IBM" ps 3, \
        f(x) t "Analytical" lw 3, \
          g(x) t "Interface" lw 1 lc rgb "black"
