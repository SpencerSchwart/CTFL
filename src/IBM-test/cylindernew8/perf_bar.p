set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "perf-bar.pdf"
set size ratio 0.8

set style line 2 lc rgb 'black' lt 1 lw 1
set style data histogram
set style histogram cluster gap 1
set style fill pattern border -1
set boxwidth 0.9
set xtics format ""
set grid ytics

set xrange [0:9]

set xlabel "Re"
set ylabel "t_i (seconds)"

set key top left

plot "perf.dat" using 2:xtic(1) title "N_i=1" ls 2,   \
     "perf.dat" using 3 title "N_i=5" ls 2,   \
     "perf.dat" using 4 t "N_i=10" ls 2,    \
