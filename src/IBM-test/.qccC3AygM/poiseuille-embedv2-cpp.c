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
#define BGHOSTS 2
#include "common.h"
#include "grid/quadtree.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "poiseuille-embedv2.c"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define L0 20
int maxlevel = 7;

const double U0 = 1;
double t_end = 100;
double H = 5;
const double delta = 20/pow(2,7);
const double int_frac = 1;
double y_min = 2.5;
double y_max = 2.73;
double tt = delta * (1 - int_frac);

face vector muc;

u.n[left] = neumann (0);
p[left]   = dirichlet (1);
pf[left]   = dirichlet (1);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

int main() {
  size(L0);
  init_grid (1 << (maxlevel));
  origin (0, -(L0/2));
  TOLERANCE = 1.e-5;
  DT = 0.1;
  mu = muc;

  run();
}


vector uref[];
event init (t = 0) {
  solid (cs, fs, fabs(y) < (H/2));
  // foreach()
   // cs[] = fabs(y) < y_max && fabs(y) > y_min? (1 - int_frac): cs[];
  u.t[embed] = dirichlet(0);
  u.n[embed] = dirichlet(0);
  double U_max = (sq(H)*1)/(8*1*L0);
  foreach() {
    u.x[] = (cs[])*U0;
    uref.x[] = fabs(y) <= y_min? U_max*(1 - sq((2*y/(H + 2*tt)))): 0;
  }
}

vector fmm[]
event properties (i++) {
  
  foreach_face() {
    muc.x[] = fm.x[];
    fmm.x[] = fm.x[];
  }
  boundary ((scalar*){muc});
}



event profile (t = t_end) {
  char name[80];

  sprintf (name, "vprof-%g", int_frac);
  FILE * fv = fopen (name, "w");
  for (double i = -L0/2; i <= L0/2; i += delta)
    foreach_point (L0 - (delta), i) {
      fprintf (fv, "%g %g %g\n", x, y, u.x[]);
    }
  fclose (fv);

  sprintf (name, "vprofref-%g", int_frac);
  FILE * fref = fopen (name, "w");
  for (double i = -L0/2; i <= L0/2; i += delta)
    foreach_point (L0 - (delta), i)
      fprintf (fref, "%g %g %g\n", x, y, uref.x[]);
  fclose (fref);
}


scalar e[];
scalar perr[];

event logfile (i++, t <= t_end) {
  trash ({e});
  
  double U_max = (sq(H)*1)/(8*1*L0);
  foreach() {
    double uref = U_max*(1 - sq((2*y/(H + 2*tt))));
    e[] = cs[] >= (1 - int_frac)? fabs(u.x[] - uref): 0;
    perr[] = cs[] >= (1 - int_frac)? fabs(u.x[] - uref)/uref: 0;
  }
  fprintf (stderr, "%d %g\n", i, t);
}



#endif
