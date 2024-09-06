#include "navier-stokes/centered.h"
#include "../immersed-new.h"
#include "view.h"

#define L0 1
#define LEVEL 8

int maxlevel = 8;
// int Re = 100;

double t_end = 10;
double U0 = 1;
double x_max = 0.7481;
double x_min = 0.74415;

coord vc = {0, 0};
coord ci = {0.75, 0};  // coordinate of wall

scalar vof[];
face vector sf[];
// face vector muv[];

u.n[left]   = dirichlet (U0);
u.t[left]   = dirichlet (0);
p[left]     = neumann (0);
pf[left]    = neumann (0);

u.n[right]  = dirichlet (0);
u.t[right]  = dirichlet (0);
p[right]    = neumann (0);
pf[right]   = neumann (0);

u.n[top]    = neumann (0);
p[top]      = dirichlet (0);
pf[top]     = dirichlet (0);

u.n[bottom] = neumann (0);
p[bottom]   = dirichlet (0);
pf[bottom]  = dirichlet (0);


int main()
{
    size(L0);
    init_grid (1 << (LEVEL - 2));
    // mu = muv;
    DT = 0.01;

    run();
}


event tracer_advection (i++)
{
    solid (vof, sf, x >= ci.x);
    foreach()
        vof[] = x < x_max && x > x_min? 0.49: vof[];
}


event init (t = 0)
{
    solid (vof, sf, x >= ci.x);
    foreach()
        u.x[] = (1 - vof[])*U0;
}

/*
event properties (i++)
{
    foreach_face()
        muv.x[] = fm.x[]*(U0)*(L0)/(Re);
    boundary ((scalar *) ({muv}));
}
*/

event logfile (i++, t <= t_end)
{
    fprintf (stderr, "%d %g\n", i, t);
}

event adapt (i++)
{
    adapt_wavelet ({vof,u}, (double[]){1e-5,3e-3,3e-3},
                   maxlevel = LEVEL, minlevel = LEVEL - 3);
}

