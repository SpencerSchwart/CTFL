set term pngcairo size 1600,1000 enhanced
set output 'x2.png'
set size ratio 0.8
set xrange [9.7:10.3 ]

# set title "Re = 1" font ",16"
set xlabel "x" font ",15"
set ylabel "U.x" font ",15"
set tics font ",12"
set pointsize 1.5

a = 3 # x = 2 y = 3
b = 4
wide = 3 # line width

set style line 1 lt 1 lc rgb "black" lw wide dt 4
set style line 2 lt 1 lc rgb "blue" lw wide dt 5
set style line 3 lt 1 lc rgb "red" lw wide
set style line 4 lt 1 lc rgb "green" lw wide
set style line 5 lt 1 lc rgb "orange" lw wide


set macros

NOXTICS = "set format x ''; unset xlabel"
NOYTICS = "unset ylabel"
XTICS = "set format x '%g'; set xlabel 'y' font ',15'"
YTICS = "set format y '%g'; set ylabel 'u.x' font ',15'"

TMARGIN = "set tmargin at screen 0.95; set bmargin at screen 0.60"
BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.20"
LMARGIN = "set lmargin at screen 0.053; set rmargin at screen 0.35"
CMARGIN = "set lmargin at screen 0.37; set rmargin at screen 0.75"
RMARGIN = "set lmargin at screen 0.69; set rmargin at screen 0.99"
POS = "at graph 0.45,0.9 font ',12'"

set multiplot layout 2,3 columnsfirst
set key top center
unset key

@TMARGIN; @LMARGIN; @NOXTICS; @YTICS
set label 1 'Re = 1' @POS
plot 'embed/vprofx2-1' u a:b t 'embed' ls 1, \
        'ibm/vprofx2-1' u a:b t 'old IBM' ls 5, \
            'ibmnz/vprofx2-1' u a:b t 'bilinear IBM' ls 3, \
                'embed/vprofx2-1' u a:b notitle w l ls 1 lw 1

@BMARGIN; @LMARGIN; @XTICS; @YTICS
set label 1 'Re = 2' @POS
plot 'embed/vprofx2-2' u a:b t 'embed'  ls 1, \
        'ibm/vprofx2-2' u a:b t 'old IBM' ls 5, \
            'ibmnz/vprofx2-2' u a:b t 'bilinear IBM' ls 3, \
                'embed/vprofx2-2' u a:b notitle w l ls 1 lw 1

@TMARGIN; @CMARGIN; @NOXTICS; @NOYTICS
set label 1 'Re = 5' @POS
plot 'embed/vprofx2-5' u a:b t 'embed' ls 1, \
        'ibm/vprofx2-5' u a:b t 'old IBM' ls 5, \
            'ibmnz/vprofx2-5' u a:b t 'bilinear IBM' ls 3, \
                'embed/vprofx2-5' u a:b notitle w l ls 1 lw 1

@BMARGIN; @CMARGIN; @XTICS; @NOYTICS
set label 1 'Re = 10' @POS
plot 'embed/vprofx2-10' u a:b t 'embed' ls 1, \
        'ibm/vprofx2-10' u a:b t 'old IBM' ls 5, \
            'ibmnz/vprofx2-10' u a:b t 'bilinear IBM' ls 3, \
                'embed/vprofx2-10' u a:b notitle w l ls 1 lw 1

@TMARGIN; @RMARGIN; @NOXTICS; @NOYTICS
set label 1 'Re = 20' @POS
plot 'embed/vprofx2-20' u a:b t 'embed' ls 1, \
        'ibm/vprofx2-20' u a:b t 'old IBM' ls 5, \
            'ibmnz/vprofx2-20' u a:b t 'bilinear IBM' ls 3, \
                'embed/vprofx2-20' u a:b notitle w l ls 1 lw 1

@BMARGIN; @RMARGIN; @XTICS; @NOYTICS
set label 1 'Re = 40' @POS
plot 'embed/vprofx2-40' u a:b t 'embed' ls 1, \
        'ibm/vprofx2-40' u a:b t 'old IBM' ls 5, \
            'ibmnz/vprofx2-40' u a:b t 'bilinear IBM' ls 3, \
                'embed/vprofx2-40' u a:b notitle w l ls 1 lw 1

unset multiplot

