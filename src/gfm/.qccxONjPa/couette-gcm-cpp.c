@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#include "common.h"
#include "grid/quadtree.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "couette-gcm.c"
#include "navier-stokes/centered.h"
#include "../ibm-gfm.h"
#include "view.h"

#define L0 10
int maxlevel = 7;

const double U0 = 5;
double Re = 50;
double t_end = 100;
double H = 5;
coord vc = {0, 0};

scalar vof[];
face vector sf[];
face vector muv[];

u.n[left] = dirichlet (y > 5? U0: neumann (0));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (U0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);


int main() {
  size(L0);
  init_grid (1 << (8));
  mu = muv;
  TOLERANCE = 1.e-5;

  run();
}


event moving_wall (i++) {
  solid (vof, sf, y < (H));
}

event init (t = 0) {
  foreach()
    u.x[] = vof[]*1;
}

event propertires (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(H)/(Re);
  boundary ((scalar *) {muv});
}

event profile (t = t_end) {
  double delta = L0/pow(2,maxlevel);
  char name[80];

  sprintf (name, "vprof");
  FILE * fv = fopen (name, "w");
  for (double i = 0; i <= L0; i += delta)
    foreach_point (L0/2, i)
      fprintf (fv, "%g %g %g\n", x, y, u.x[]);
  fclose (fv);
}

event logfile (i++, t <= t_end) {
  fprintf (stderr, "%d %g\n", i, t);
}

event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = maxlevel, minlevel = maxlevel - 2);
}



#endif
