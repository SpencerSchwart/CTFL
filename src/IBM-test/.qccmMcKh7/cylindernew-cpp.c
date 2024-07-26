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

#define L0 20.48
#define D 0.5
#define LEVEL 10

int maxlevel = 10;
int Re;
double U0 =  1.; // inlet velocity
double t_end = 50;
double tf_start = 5;
coord ci = {5, 10}; // initial coordinates of cylinder
coord vc = {0, 0}; // velocity of cylinder

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

  Re = 185;
  run();

  /*
  Re = 2;
  run();

  Re = 5;
  run();

  Re = 10;
  run();

  Re = 20;
  run();

  Re = 40;
  run();
  
  Re = 50;
  run();
  */
}


event moving_cylinder (i++) {
  solid (vof, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
}


event init (t = 0) { 
  // mask(y > 6 ? top: y < -6 ? bottom : none);
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}

vector Uc[];
vector ub[];
event logfile (i++){
  coord Fp = {0,0};
  immersed_forcev2 (vof, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*(D));
  double CL = (Fp.y)/(0.5*sq(U0)*(D));

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
      
      // fprintf(fpv, "%d %g %g %g %g %g %g %g %g %g %g\n", i, x, y, u.x[], u.y[], pc.x, pc.y, uc.x, uc.y, ub.x[], ub.y[]);
    }
    else
      foreach_dimension()
        Uc.x[] = 0;
        ub.x[] = 0;
  }

 fclose(fpv);

 double uavg_x = usum.x/counter;
 double uavg_y = usum.y/counter;
 fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g\n",
	   i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, uavg_x, uavg_y);


}

event profile (t = t_end) {
  int k = 0;
  double delta = 15/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "vprofx1-%d", Re); // x = 4.8125
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

  sprintf (name, "vprofx2-%d", Re); // x = 5
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

  sprintf (name, "vprofx3-%d", Re); // x = 10
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

  sprintf (name, "vprofy1-%d", Re); // y = 3
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

  sprintf (name, "vprofy2-%d", Re); // y = 3.1875
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

  sprintf (name, "surface-%d", Re);
  FILE * fv5 = fopen (name, "w");
  foreach() {
    if (vof[] > 0 && vof[] < 1) {
      coord pc, ucb;
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      bilinear_interpolation (point, u, pc, &ucb);
      double pcb = scalar_bilinear_interpolation (point, p, pc);


      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      double mag = sqrt(sq(ucb.x) + sq(ucb.y));
      double cp = pcb / (0.5 * sq(U0));
      fprintf (fv5, "%g %g %g\n", fabs(theta), mag, cp);
    }
  }
  fflush (fv5);
  fclose (fv5);

}


event snapshot (t = t_end) {
  scalar omega[];
  vorticity (u, omega);

  char name[80];
  sprintf (name, "vort-%d.png", Re);
  view (fov = 2, tx = -0.275, ty = -0.495,
	width = 800, height = 400); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  save (name);
}


event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}

int count = 0;
event frequency (i++; t >= tf_start) {
  char name[80];
  sprintf (name, "freq.dat");
  FILE * fp = fopen (name, "a")
  foreach_point(7.500, ci.y) {
    fprintf (fp, "", count, x, y, u.x[], u.y[], p[]);
    count++;
  }
}

event movie (t += 0.01; t <= t_end) {
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}

#endif
