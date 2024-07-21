set term pngcairo size 1600,1000 enhanced
set output 'y1.png'
set size ratio 0.8
# set xrange [0:15]
 set xrange [4.5:5.5]
# set yrange [:0.7]

# set title "Re = 1" font ",16"
set xlabel "x" font ",15"
set ylabel "U.x" font ",15"
set tics font ",12"
set pointsize 1.5


set key top center
# set tmargin 0
# set bmargin 0
# set rmargin 3
# set lmargin 8

a = 2 # x = 2 y = 3
b = 4
wide = 2 # line width

set style line 1 lt 1 lc rgb "black" lw wide dt 4
set style line 2 lt 1 lc rgb "blue" lw wide dt 5
set style line 3 lt 1 lc rgb "red" lw wide
set style line 4 lt 1 lc rgb "green" lw wide
set style line 5 lt 1 lc rgb "orange" lw wide

set macros

NOXTICS = "set format x ''; unset xlabel"
NOYTICS = "unset ylabel"
XTICS = "set format x '%g'; set xlabel 'x' font ',15'"
YTICS = "set format y '%g'; set ylabel 'u.x' font ',15'"

TMARGIN = "set tmargin at screen 0.95; set bmargin at screen 0.60"
BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.20"
LMARGIN = "set lmargin at screen 0.053; set rmargin at screen 0.35"
CMARGIN = "set lmargin at screen 0.37; set rmargin at screen 0.75"
RMARGIN = "set lmargin at screen 0.69; set rmargin at screen 0.99"
POS = "at graph 0.1,0.1 font ',12'"

set multiplot layout 2,3 columnsfirst
# unset key

@TMARGIN; @LMARGIN; @NOXTICS; @YTICS
set label 1 'Re = 1' @POS
plot 'ce11/vprofy1-10' u a:b t 'embed' w l ls 1, \
        'c11/vprofy1-10' u a:b t 'old IBM' w l ls 2, \
            'ibmv4/vprofy1-10' u a:b t 'non-iterative IBM w/z v2' w l ls 5, \
                'ibmv5/vprofy1-10' u a:b t 'non-iterative IBM w/z' w l ls 4

@BMARGIN; @LMARGIN; @XTICS; @YTICS
set label 1 'Re = 2' @POS
plot 'ce11/vprofy1-11' u a:b t 'embed' w l ls 1, \
        'c11/vprofy1-11' u a:b t 'old IBM' w l ls 2, \
            'ibmv4/vprofy1-11' u a:b t 'non-iterative IBM w/z v2' w l ls 5, \
                'ibmv5/vprofy1-11' u a:b t 'non-iterative IBM w/z' w l ls 4

@TMARGIN; @CMARGIN; @NOXTICS; @NOYTICS
set label 1 'Re = 5' @POS
plot 'ce11/vprofy1-12' u a:b t 'embed' w l ls 1, \
        'c11/vprofy1-12' u a:b t 'old IBM' w l ls 2, \
            'ibmv4/vprofy1-12' u a:b t 'non-iterative IBM w/z v2' w l ls 5, \
                'ibmv5/vprofy1-12' u a:b t 'non-iterative IBM w/z' w l ls 4

@BMARGIN; @CMARGIN; @XTICS; @NOYTICS
set label 1 'Re = 10' @POS
plot 'ce11/vprofy1-13' u a:b t 'embed' w l ls 1, \
        'c11/vprofy1-13' u a:b t 'old IBM' w l ls 2, \
            'ibmv4/vprofy1-13' u a:b t 'non-iterative IBM w/z v2' w l ls 5, \
                'ibmv5/vprofy1-13' u a:b t 'non-iterative IBM w/z' w l ls 4

@TMARGIN; @RMARGIN; @NOXTICS; @NOYTICS
set label 1 'Re = 20' @POS
plot 'ce11/vprofy1-14' u a:b t 'embed' w l ls 1, \
        'c11/vprofy1-14' u a:b t 'old IBM' w l ls 2, \
            'ibmv4/vprofy1-14' u a:b t 'non-iterative IBM w/z v2' w l ls 5, \
                'ibmv5/vprofy1-14' u a:b t 'non-iterative IBM w/z' w l ls 4

@BMARGIN; @RMARGIN; @XTICS; @NOYTICS
set label 1 'Re = 40' @POS
plot 'ce11/vprofy1-15' u a:b t 'embed' w l ls 1, \
        'c11/vprofy1-15' u a:b t 'old IBM' w l ls 2, \
            'ibmv4/vprofy1-15' u a:b t 'non-iterative IBM w/z v2' w l ls 5, \
                'ibmv5/vprofy1-15' u a:b t 'non-iterative IBM w/z' w l ls 4

unset multiplot

