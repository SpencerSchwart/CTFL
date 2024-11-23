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
//#include "navier-stokes/centered.h"
#include "../ibm-gfm.h"
#include "../my-centered.h"
#include "view.h"


#define L0 10
int maxlevel = 7;

const double U0 = 5;
double Re = 10;
double t_end = 100;
double H = 10;
// coord vc = {0, 0};
const double delta = 10./pow(2,7);
const double int_frac = 0.75;
double y_min = 4.90;
double y_max = 5.02;
double tt = delta * (1 - int_frac);

// scalar vof[];
// face vector sf[];
face vector muv[];

u.n[left] = neumann (0);
// u.n[left] = dirichlet (1);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (U0);
p[top] = neumann (0);
pf[top] = neumann (0);

// u.n[bottom] = dirichlet (0);
// u.t[bottom] = dirichlet (0);
// u.t[bottom] = neumann (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

int main() {
  size(L0);
  init_grid (1 << (maxlevel));
  mu = muv;
  // TOLERANCE = 1.e-5;
  stokes = true;
  DT = 0.1;

  run();
}

bid solid;
u.t[solid]=dirichlet(0.);

scalar cy[];
face vector cyf[];
bid cylinder;
u.t[solid] = dirichlet(0.);

event wall (i=0) {
  solid (vof, sf, y > (5));
  // solid (cy, cyf, - sq(x-5-t*(0.01)) - sq(y-7.5) + sq(1));
  //mask (cy[] == 1? cylinder: none);
  //mask (y < 4.90? solid : none);
  foreach()
    vof[] = fabs(y) < y_max && fabs(y) > y_min? int_frac: vof[];
  boundary ({vof, sf});
}

face vector mf[];
face vector alphamf[];
scalar mc[];
scalar rhot[];

event propertires (i++) {
  foreach()
    mc[] = cm[];
    rhot[] = rho;
  foreach_face() {
      mf.x[] = fm.x[];
      alphamf.x[] = alpha.x[];
      muv.x[] = fm.x[]*(U0)*(H)/(Re);
      // muv.x[] = alpha.x[]*1;
  }
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
  fprintf (stderr, "%d %g %d %d %d %d\n", i, t,
           mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
}

/*
event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = maxlevel, minlevel = maxlevel - 2);
}
*/


#endif
