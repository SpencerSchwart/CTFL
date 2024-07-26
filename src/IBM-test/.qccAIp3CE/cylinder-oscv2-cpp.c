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
#line 1 "cylinder-oscv2.c"
#include "navier-stokes/centered.h"
#include "../immersed.h" 
#include "view.h"

#define L0 40.
#define D 0.50
#define LEVEL 12

int maxlevel = 12;
int Re;
double U0 =  1.; // inlet velocity
double t_end = 50;
coord ci = {5, 20}; // initial coordinates of cylinder
coord vc = {0, 0}; // velocity of cylinder
coord xc = {0, 0};
const double A = 0.2*D;
double freq = 0.3;

scalar vof[];
face vector sf[];
face vector muv[];

u.n[left] = dirichlet ((U0));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = neumann (0);
u.t[top] = dirichlet (U0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = neumann (0);
u.t[bottom] = dirichlet (U0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);


int main() {
  size(L0);
  init_grid (2 << (6));
  mu = muv;
  TOLERANCE = 1.e-6; 
  // DT = 0.003;
  // CFL = 0.6;

  Re = 185;
  run();
}


event moving_cylinder (i++) {
  xc.y = A*sin(2*M_PI*freq*t);
  vc.y = 2*M_PI*freq*A*cos(2*M_PI*freq*t); 
  solid (vof, sf, - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2));
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}

scalar e[];
vector Uc[];
vector ub[];
double avgCD = 0, avgCL = 0;
int count = 0;

event logfile (i++){

  xc.y = A*sin(2*M_PI*freq*t);
  vc.y = 2*M_PI*freq*A*cos(2*M_PI*freq*t); 

  coord Fp = {0,0};
  immersed_forcev2 (vof, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*(D));
  avgCD += t > 20? CD: 0;
  double CL = (Fp.y)/(0.5*sq(U0)*(D));
  avgCL += t > 20? CL: 0;
  count += t > 20? 1:0;

  char name[80];
  sprintf (name, "avg_v.txt");
  FILE * fpv = fopen(name, "a");
  coord total = {0};

  int counter = 0;
  coord usum = {0};
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

      coord sum = {0};
        foreach_neighbor() {
          double delta_u = delta_func (x, y, pc.x, pc.y, Delta);
          foreach_dimension()
	    sum.x += u.x[]*delta_u*dv();
        }
      bilinear_interpolation (point, u, pc, &ucb);

      foreach_dimension() {
        uc.x = sum.x;
	Uc.x[] = ucb.x;
	usum.x += ucb.x;
      }
      e[] = fabs(u.y[] - vc.y) / vc.y;
    }
    else
      foreach_dimension() {
        Uc.x[] = 0;
        ub.x[] = 0;
      }
  }

  norm n = normf (e);
  double uavg_x = usum.x/counter;
  double uavg_y = usum.y/counter;

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
	   i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax,
	   CD, avgCD/(count + 1e-6), CL, avgCL/(count + 1e-6),
	   uavg_x, uavg_y, xc.y, vc.y, n.avg, n.max);
}

event profile (t = t_end) {
  double delta = 20/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "surface-%d", Re);
  FILE * fv5 = fopen (name, "w");
  foreach() {
    if (vof[] > 0 && vof[] < 1) {
      coord pc, ucb, fcb;
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      bilinear_interpolation (point, u, pc, &ucb);
      bilinear_interpolation (point, fc, pc, &fcb);
      double pcb = scalar_bilinear_interpolation (point, p, pc);


      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI) + 180;
      double umag = sqrt(sq(ucb.x) + sq(ucb.y));
      double fmag = sqrt(sq(fcb.x) + sq(fcb.y));
      double cp = pcb / (0.5 * sq(U0));
      fprintf (fv5, "%g %g %g %g\n", fabs(theta), umag, fmag, cp);
    }
  }
  fflush (fv5);
  fclose (fv5);

}

scalar omega[];

event movie (t += 0.1) {
  vorticity (u, omega);
  char name[80];
  sprintf (name, "vort-%d-%g.png", Re, t);
  view (fov = 3, tx = 0 ty = 0, bg = {1,1,1,});
  isoline ("omega", min = -10, max = 10);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  squares ("omega", map = cool_warm, min = -10, max = 10);
  save (name);
}


event snapshot (t = t_end) {
  vorticity (u, omega);

  char name[80];
  sprintf (name, "vort-%d.png", Re);
  view (fov = 2, tx = -(L0/ci.x), ty = -(L0/ci.y),
	width = 800, height = 400); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  save (name);
}


event adapt (i++) {
  vorticity (u, omega);
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}


event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}

#endif
