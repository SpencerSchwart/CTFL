set xlabel "Theta (degrees)"
set ylabel "u.y"

a = 3
b = 5  # u.x=4, u.y=5, p=6
wide = 2

set style line 1 lt 1 lc rgb "violet" lw wide
set style line 2 lt 1 lc rgb "blue" lw wide
set style line 3 lt 1 lc rgb "cyan" lw wide 
set style line 4 lt 1 lc rgb "green" lw wide
set style line 5 lt 1 lc rgb "orange" lw wide
set style line 6 lt 1 lc rgb "black" lw wide dt 4
set style line 7 lt 1 lc rgb "red" lw wide dt 5

plot '../ibmcylinder/ring2-40' u a:b t 'Accel' ls 1, \
        '../ibmcylinder1/ring2-40' u a:b t 'Adv+Accel' ls 2, \
            '../ibmcylinder2/ring2-40' u a:b t 'endtime' ls 3, \
                '../ibmcylinder3/ring2-40' u a:b t 'Adv+endtime' ls 4, \
                        '../ibmcylinder4/ring2-40' u a:b t 'Adv' ls 5, \
                            '../../IBM-test/cylindernew5/ring2-40' u a:b t 'DF IBM' ls 6, \
                                '../../IBM-test/cylinderE10/ring2-40' u a:b t 'Embed'ls 7
