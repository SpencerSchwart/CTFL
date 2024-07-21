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
#line 1 "cylinderO.c"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "../immersed.h" // IBM
#include "curvature.h"
#include "view.h"

const double D = 0.5 [1];
#define LEVEL 9
#define RE 16.7
const double ST = 0.625;

const double U0 =  1 [1,-1];
double pulsation = ST*U0/(D/2);
double A0 = (0.8*D);
double t_end = 20 [0,1];
coord ci = {5, 0}; // initial coordinates of cylinder
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
  size(20*D);
  init_grid (2 << (LEVEL-4));
  mu = muv;
  TOLERANCE = 1.e-5 [*]; 

  DT = 0.0005;
  run();
}

event moving_cylinder (i++) {
  vc.x = A0*pulsation*sin(pulsation*t); // velocity of cylinder
  xc.x = - A0*cos(pulsation*t);
  solid (airfoil, sf, - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2));
}

event init (t = 0) {
  mask(y > 10*D ? top: y < -10*D ? bottom : none);

  foreach()
    foreach_dimension()
      u.x[] = 0.;
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(RE);
   boundary ((scalar *) {muv});
}


event logfile (i++; t <= t_end){
  coord Fp = {0,0};
  immersed_force (airfoil, sf, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*pi*sq(D/2));
  double CL = (Fp.y)/(0.5*sq(U0)*pi*sq(D/2));

  fprintf (stderr, "%d %g %d %d %d %d %g %g %g %g %g\n",
	   i, t, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, vc.x, vc.y, dt);
	   
}


event movie (t += 0.05) {
  scalar omega[];
  vorticity (u, omega);
  char name[80];
  sprintf (name, "snap-%g.png", t);
  FILE * fp1 = fopen(name, "w");
  view (fov = 10, tx = -0.5, width = 600, height = 675);
  draw_vof ("airfoil", "sf", lw = 3);
  squares ("omega", map = cool_warm);
  cells();
  mirror ({0,1}) {
    draw_vof("airfoil", "sf", lw = 3);
    squares ("p", map = cool_warm);
    cells();
  }
  save (fp = fp1);
  fclose (fp1);
}


event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 4);
}


#endif
