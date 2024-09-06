set xlabel "Theta (degrees)" font ",13"
set ylabel "u.x" font ",13"

# set xrange [0:360]
set xrange [-180:180]

set key top center

a = 3
b = 4
c = b - 3

set title "Re = 20" font ",15"

# plot 'surfaceVelocity2-185' u a:b t "IBM lvl 11", \
#        '../cylinderE10/surfaceVelocity2-185' u a:b t "embed", \
#            '../cylindernewlvl10/surfaceVelocity2-185' u a:b t "IBM lvl 10"

#plot '../cylinderl11/surface-100' u 1:7 t "embed" lc rgb "black", \
#        'surface-100' u 1:b t "IBM HD", \
#            '../cylindernew1/surface-100' u 1:b t "IBM FD", \
#                '../cylindernew2/surface-100' u 1:b t "IBM FD (before Projection)"

plot '../cylinderl11/surfaceVelocity1-20' u a:b t "embed" lc rgb "black", \
        'surfaceVelocity1-20' u a:b t "IBM HD", \
            '../cylindernew1/surfaceVelocity1-20' u a:b t "IBM FD", \
                '../cylindernew2/surfaceVelocity1-20' u a:b t "IBM FD (before Projection)"
            
