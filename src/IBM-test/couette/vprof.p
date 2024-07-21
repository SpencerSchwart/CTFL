set xrange [0:]
set yrange [0:5]

set title "U = 5"

set xlabel "u.x"
set ylabel "y"

f(y) = y

set key center top

plot 'vprof' u 3:2 t "IBM", \
        f(x) t "Equation"
