set term pngcairo size 1600,1000 enhanced
set output 'y1.png'
set size ratio 0.8
# set xrange [2.5:3.5]
# set xrange [0:6]
#set xrange [0:15]
set xrange [4.5:5.5]
# set yrange [:0.6]

# set title "Re = 1" font ",16"
set xlabel "y" font ",15"
set ylabel "U.x" font ",15"
set tics font ",12"
set pointsize 1.5


# set key top center
# set tmargin 0
# set bmargin 0
# set rmargin 3
# set lmargin 8

a = 2
b = 4
wide = 2 # line width

set macros

NOXTICS = "set format x ''; unset xlabel"
NOYTICS = "set format y '%g'; unset ylabel"
XTICS = "set format x '%g'; set xlabel 'y' font ',15'"
YTICS = "set format y '%g'; set ylabel 'u.x' font ',15'"

TMARGIN = "set tmargin at screen 0.95; set bmargin at screen 0.60"
BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.20"
LMARGIN = "set lmargin at screen 0.05; set rmargin at screen 0.35"
CMARGIN = "set lmargin at screen 0.35; set rmargin at screen 0.65"
RMARGIN = "set lmargin at screen 0.65; set rmargin at screen 0.95"
POS = "at graph 0.4,0.95 font ',12'"

set multiplot layout 2,3 columnsfirst
unset key

@TMARGIN; @LMARGIN; @NOXTICS; @YTICS
set label 1 'Re = 1' @POS
plot 'c11v3/vprofy1-10' u a:b t 'Single P' w l lw wide lc "orange", \
   'c11dpv1/vprofy1-10' u a:b t 'Double Pv1' w l lw wide lc "blue", \
    'c11dpv2/vprofy1-10' u a:b t 'Double Pv2' w l lw wide lc "green", \
      'ce11/vprofy1-10' u a:b t 'embed' w l lw wide lc "red", \

@BMARGIN; @LMARGIN; @XTICS; @YTICS
set label 1 'Re = 2' @POS
plot 'c11v3/vprofy1-11' u a:b t 'Single P' w l lw wide lc "orange", \
   'c11dpv1/vprofy1-11' u a:b t 'Double Pv1' w l lw wide lc "blue", \
    'c11dpv2/vprofy1-11' u a:b t 'Double Pv2' w l lw wide lc "green", \
      'ce11/vprofy1-11' u a:b t 'embed' w l lw wide lc "red"

@TMARGIN; @CMARGIN; @NOXTICS; @NOYTICS
set label 1 'Re = 5' @POS
plot 'c11v3/vprofy1-12' u a:b t 'Single P' w l lw wide lc "orange", \
   'c11dpv1/vprofy1-12' u a:b t 'Double P' w l lw wide lc "blue", \
    'c11dpv2/vprofy1-12' u a:b t 'no SPM' w l lw wide lc "green", \
      'ce11/vprofy1-12' u a:b t 'embed' w l lw wide lc "red"

@BMARGIN; @CMARGIN; @XTICS; @NOYTICS
set label 1 'Re = 10' @POS
plot 'c11v3/vprofy1-13' u a:b t 'SPM v1' w l lw wide lc "orange", \
   'c11dpv1/vprofy1-13' u a:b t 'SPM v2' w l lw wide lc "blue", \
    'c11dpv2/vprofy1-13' u a:b t 'no SPM' w l lw wide lc "green", \
      'ce11/vprofy1-13' u a:b t 'embed' w l lw wide lc "red"

@TMARGIN; @RMARGIN; @NOXTICS; @NOYTICS
set label 1 'Re = 20' @POS
plot 'c11v3/vprofy1-14' u a:b t 'SPM v1' w l lw wide lc "orange", \
   'c11dpv1/vprofy1-14' u a:b t 'SPM v2' w l lw wide lc "blue", \
    'c11dpv2/vprofy1-14' u a:b t 'no SPM' w l lw wide lc "green", \
      'ce11/vprofy1-14' u a:b t 'embed' w l lw wide lc "red"

@BMARGIN; @RMARGIN; @XTICS; @NOYTICS
set label 1 'Re = 40' @POS
plot 'c11v3/vprofy1-15' u a:b t 'SPM v1' w l lw wide lc "orange", \
   'c11dpv1/vprofy1-15' u a:b t 'SPM v2' w l lw wide lc "blue", \
    'c11dpv2/vprofy1-15' u a:b t 'no SPM' w l lw wide lc "green", \
      'ce11/vprofy1-15' u a:b t 'embed' w l lw wide lc "red"

unset multiplot

