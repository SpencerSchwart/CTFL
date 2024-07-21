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
#line 1 "couettev2.c"
#include "navier-stokes/centered.h"
#include "../immersed.h"
#include "view.h"

#define L0 10
int maxlevel = 6;

const double U0 = 5;
double Re = 50;
double t_end = 100;
double H = 5;
const double delta = 10/pow(2,6);
double tt = 3*delta/4;
double y_min = 0;
double y_max = -0.23;

coord vc = {0, 0};

scalar vof[];
face vector sf[];
face vector muv[];

u.n[left] = dirichlet (y < y_max? 0: neumann (0));
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
  init_grid (1 << (maxlevel));
  origin (0, -5);
  mu = muv;
  // TOLERANCE = 1.e-5;

  run();
}


event moving_wall (i++) {
  solid (vof, sf, y < y_min);
  foreach() 
    vof[] = y < y_max && y > y_min? 0.5: vof[];
}

event init (t = 0) {
  foreach()
    u.x[] = 1 - vof[];
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(5)/(Re);
  boundary ((scalar *) {muv});
}

event profile (t = t_end) {
  char name[80];

  sprintf (name, "vprof4");
  FILE * fv = fopen (name, "w");
  for (double i = -L0/2; i <= L0/2; i += delta)
    foreach_point (L0/2, i) {
      coord ucb = {0};
      if (vof[] > 0 && vof[] < 1) {
        coord pc, ucb;
        coord m = {x, y};
        coord n = interface_normal (point, vof);
        double alpha = plane_alpha (vof[], n);
        double area = plane_area_center(n, alpha, &pc);

        foreach_dimension()
          pc.x = m.x + pc.x*Delta; 

        bilinear_interpolation (point, u, pc, &ucb);
      } else if (vof[] == 1 && vof[0,1] == 0){
        double xc = x;
        double yc = y + delta/2;
        coord pc = {xc, yc}, ucb;
        bilinear_interpolation (point, u, pc, &ucb);
    }
        fprintf (fv, "%g %g %g %g\n", x, y, u.x[], ucb.x);
    }
  fclose (fv);
}

event logfile (i++, t <= t_end) {

  int counter = 0;
  coord usum = {0};
  coord pe = {0};
  double max = -HUGE;
  scalar e[];
  foreach() {
    if (vof[] > 0 && vof[] < 1) {
      counter++;
      coord pc, uc, ucb;
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

     foreach_dimension()
        pc.x = m.x + pc.x*Delta; 

     bilinear_interpolation (point, u, pc, &ucb);

      foreach_dimension() 
        usum.x += ucb.x;
    } else if (vof[] == 1 && vof[0,1] == 0){
      counter++;
      double xc = x;
      double yc = y + delta/2;
      coord pc = {xc, yc}, ucb;
      bilinear_interpolation (point, u, pc, &ucb);

      foreach_dimension() 
        usum.x += ucb.x;
    }
    double ref = U0*y/(H + tt);
    e[] = vof[] < 1? u.x[] - ref: 0;
    if (e[] > max) {
      max = e[];
      pe.x = x;
      pe.y = y;
    }
  }
  norm n = normf (e);

  double uavg_x = counter > 0? usum.x/counter: 0;
  double uavg_y = counter > 0? usum.y/counter: 0;
  fprintf (stderr, "%d %g %g %g %g %g %g %g %g %g\n", i, t, uavg_x, uavg_y, n.avg, n.rms, n.max, max, pe.x, pe.y);
}

/*
event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = maxlevel, minlevel = maxlevel - 2);
}
*/


#endif
