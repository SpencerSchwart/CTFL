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
#line 1 "cylindernew.c"
#include "navier-stokes/centered.h"
#include "../immersed.h" 
#include "view.h"

#define L0 15.
#define D 0.5
#define LEVEL 10

int maxlevel = 10;
double Re;
double U0 =  1.; // inlet velocity
double t_end = 20;
double xi = 4.8266;
double yi = 3.0835;
coord ci = {5, 3}; // initial coordinates of cylinder
coord vc = {0, 0}; // velocity of cylinder
int j;


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
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = neumann (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);


int main() {
  size(L0);
  init_grid (2 << (6));
  mu = muv;
  TOLERANCE = 1.e-5; 

  j = 10;
  Re = 40.;
  run();
}


scalar ref[];
face vector rf[];

event moving_cylinder (i++) {
  solid (vof, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
  solid (ref, rf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));

  // vof.refine = fraction_refine;
  
}

event init (t = 0) {
  
  mask(y > 6 ? top: y < -6 ? bottom : none);
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}


event logfile (i++){
  coord Fp = {0,0}, Fmu = {0,0};
  immersed_force (vof, sf, &Fp);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(D));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(D));

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL);
	   
}


event snapshot (t = t_end) {
  scalar omega[];
  vorticity (u, omega);
  solid (ref, rf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));

  char name[80];
  sprintf (name, "vort-%d", j);
  FILE * fp1 = fopen (name, "w");
  view (fov = 2, tx = -0.375, ty = -0.20,
	width = 800, height = 400); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("ref", "rf", filled = 1, lw = 5);
  save (fp = fp1);
}

event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
  // refine (vof[] > 0 && vof[] < 1 && level < LEVEL);
}

event movie (t += 0.01; t <= t_end) {
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", j, s.real, i, s.speed);
  fflush (fp);
  return 1;
}

#endif
