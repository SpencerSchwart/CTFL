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
#line 1 "cylinderE10.c"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define L0 15.
#define D 1.
#define LEVEL 9

int Re;
double U0 =  1.; // inlet velocity
double t_end = 20;
coord ci = {L0/4, L0/2}; // initial coordinates of cylinder

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

double SDF (double x, double y) {
   return - sq(x - ci.x) - sq(y - ci.y) + sq(D/2);
}

int main() {
  size(L0);
  init_grid (1 << (LEVEL - 3));
  mu = muv;
  TOLERANCE = 1.e-6 [*]; 
  CFL = 0.8;

  Re = 40;
  run();
}


event init (t = 0) {
  refine (level <= LEVEL*(1. - sqrt(fabs(sq(x-ci.x) + sq(y-ci.y) - sq(D/2.)))/2.));
  solid (cs, fs, sq(x - ci.x) + sq(y - ci.y) - sq(D/2));
  foreach()
    u.x[] = cs[] ? U0 : 0.;
  
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
 // boundary ((scalar *) {muv});
}


event logfile (i++; t <= t_end){
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(D));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(D));

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g\n",
	   i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL);
}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}

event snapshot (t = t_end) {
  FILE * file = fopen("data.txt", "w");
  output_field (fp = file);
  fclose(file);
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}

#endif
