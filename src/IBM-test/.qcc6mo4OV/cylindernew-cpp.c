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
  Re = 1.;
  run();

  j = 11;
  Re = 2.;
  run();

  j = 12;
  Re = 5.;
  run();

  j = 13;
  Re = 10.;
  run();

  j = 14;
  Re = 20.;
  run();

  j = 15;
  Re = 40.;
  run();
}


event moving_cylinder (i++) {
  solid (vof, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
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
  coord Fp = {0,0};
  immersed_forcev2 (vof, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*(D));
  double CL = (Fp.y)/(0.5*sq(U0)*(D));

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL);
}


event profile (t = t_end) {
  int k = 0;
  double delta = 15/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "vprofx1-%d", j); // x = 4.8125
  FILE * fv = fopen(name, "w");
  for(double i = 0; i <= 6; i += delta) {
    foreach_point (4.8125, i) {
      if (vof[] > 0 && vof[] < 1)
        k = 2.;
      else if (vof[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv, "%d %g %g %g %g %g\n", k, x, y, u.x[], u.y[], p[]); 
    }
  }
  fflush (fv);
  fclose (fv);

  sprintf (name, "vprofx2-%d", j); // x = 5
  FILE * fv1 = fopen(name, "w");
  for(double i = 0; i <= 6; i += delta) {
    foreach_point (5, i) {
      if (vof[] > 0 && vof[] < 1)
        k = 2.;
      else if (vof[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv1, "%d %g %g %g %g %g\n", k, x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv1);
  fclose (fv1);

  sprintf (name, "vprofx3-%d", j); // x = 10
  FILE * fv2 = fopen(name, "w");
  for(double i = 0; i <= 6; i += delta) {
    foreach_point (10, i) {
      if (vof[] > 0 && vof[] < 1)
        k = 2.;
      else if (vof[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv2, "%d %g %g %g %g %g\n", k, x, y, u.x[], u.y[], p[]); 
    }
  }
  fflush (fv2);
  fclose (fv2);

  sprintf (name, "vprofy1-%d", j); // y = 3
  FILE * fv3 = fopen(name, "w");
  for(double i = 0; i <= 15; i += delta) {
    foreach_point (i, 3) {
      if (vof[] > 0 && vof[] < 1)
        k = 2.;
      else if (vof[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv3, "%d %g %g %g %g %g\n", k, x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv3);
  fclose (fv3);

  sprintf (name, "vprofy2-%d", j); // y = 3.1875
  FILE * fv4 = fopen(name, "w");
  for(double i = 0; i <= 15; i += delta) {
    foreach_point (i, 3.1875) {
      if (vof[] > 0 && vof[] < 1)
        k = 2.;
      else if (vof[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv4, "%d %g %g %g %g %g\n", k, x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv4);
  fclose (fv4);
 
}


event snapshot (t = t_end) {
  scalar omega[];
  vorticity (u, omega);

  char name[80];
  sprintf (name, "vort-%d.png", j);
  FILE * fp1 = fopen (name, "w");
  view (fov = 2, tx = -0.375, ty = -0.20,
	width = 800, height = 400); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  save (fp = fp1);
}


event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
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
