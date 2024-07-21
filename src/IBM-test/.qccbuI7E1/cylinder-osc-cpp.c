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
#line 1 "cylinder-osc.c"
#include "navier-stokes/centered.h"
#include "../immersed.h" // IBM
#include "curvature.h"
#include "view.h"

#define LEVEL 9
double RE, KC, Amp, freq = 1 [*];
const double U0 =  1 [1,-1]; // inlet velocity
const double D = 0.5 [1];
double t_end = 62.5 [0,1];

coord ci = {0}; // initial coordinates of cylinder
coord vc = {0}; // velocity of cylinder
coord xc = {0}; // displacement of cylinder

scalar airfoil[];
face vector sf[];
face vector muv[];

u.n[left] = neumann (0);
u.t[left] = neumann (0);
p[left]   = dirichlet (0);
pf[left]  = dirichlet (0);

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
  L0 = 12.8*D;
  ci.x = L0/2;
  ci.y = L0/2;
  size(L0);
  init_grid (2 << (LEVEL-4));
  mu = muv;
  TOLERANCE = 1.e-5 [*]; 
  DT = 0.0005;

  RE = 100;
  KC = 5;
  run();
}

event moving_cylinder (i++) {
  Amp = KC*D/(2*M_PI);
  freq = 0.4 [*];
  vc.y = -Amp*2*M_PI*freq*cos(2*M_PI*freq*t);
  xc.y = -Amp*sin(2*M_PI*freq*t);
  solid (airfoil, sf, - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2));
}

event init (t = 0) {
 //  mask(y > 10*D ? top: y < -10*D ? bottom : none);

  foreach()
    foreach_dimension()
      u.x[] = 0.;
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(RE);
   boundary ((scalar *) {muv});
}

int T = 0;
event period (t += 1/freq) {
  T++;
}

event logfile (i++; t <= t_end){
  coord Fp = {0};
  immersed_force (airfoil, sf, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*D);
  double CL = (Fp.y)/(0.5*sq(U0)*D);
  
  fprintf (stderr, "%d %g %d %d %d %d %g %g %g %g %g %d\n",
	   i, t, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, xc.y, vc.y, dt, T);
	   
}

event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}



#endif
