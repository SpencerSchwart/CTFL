set terminal pdfcairo enhanced size 8,6 font "Arial,20"
set output "cd.pdf"
set size ratio 0.8
set yrange [0:6]

a = 1
b = 2

set xlabel "Re" font ",26"
set ylabel "C_D" font ",26"

set key top center font ",24"

set xtics 0,50,200
set mxtics
set mytics 2

set style line 1 lt 1 lc rgb "black" lw 3
set style line 2 lt 1 lc rgb "red" lw 3 dt 2
set style line 3 lt 1 lc rgb "blue" lw 3 dt 3

plot 'cd5.dat' u 1:2 t "Tritton (1959)"  lc rgb "black" lw 2 pt 1 ps 2, \
        'cd4.dat' u 1:2 t "Le et al. (2006)" lc rgb "dark-green" lw 2 pt 2 ps 2, \
            'cd3.dat' u 1:2 t "Park et al. (1998)" lc rgb "blue" pt 9 ps 2, \
                'cd2.dat' u 1:2 t "Calhoun (2001)" pt 13 ps 2, \
                    'cd1.dat' u 1:2 t "Present" lc rgb "red" lw 2 pt 5 ps 2



set output
