set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "cycle.pdf"
set size ratio 0.8

a = 1
b = 2

set xlabel "y_c(t)/D"
set ylabel "C_D"

set xrange [-.24:.24]
set yrange [1.17:1.5]

set key top center
unset key

set style line 1 lt 1 lc rgb "black" lw 3

plot '../cylinder-oscv4/cycle.dat' u 12:(-$14) w l ls 1

set output
