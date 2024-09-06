a = 22
b = 2 + a
c = 2 + b
d = 2 + c

plot [1:] 'log' u 2:a w l t "embed tech CD", \
    'log' u 2:b w l t "2 point CD", \
        'log' u 2:c w l t "3 point CD", \
            'log' u 2:d w l t "extrapolate CD"
#                '../cylinderE10/log' u 2:12 w l t 'embed CD'
