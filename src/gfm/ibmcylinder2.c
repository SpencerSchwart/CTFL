#include "navier-stokes/centered.h"
#include "../ibm-gfm.h" 
#include "view.h"

#define L0 15
#define D 0.5
#define LEVEL 10

int Re;
int maxlevel = 10;
double U0 =  1.;             // inlet velocity
double t_end = 50;
double tf_start = 25;
coord ci = {L0/4, L0/2};     // initial coordinates of cylinder
coord vc = {0, 0};           // velocity of cylinder

scalar vof[];
face vector sf[];
face vector muv[];

u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = neumann (0);

u.n[bottom] = neumann (0);

int main() {
  size(L0);
  init_grid (1 << (LEVEL - 2));
  mu = muv;
  TOLERANCE = 1.e-6; 
  CFL = 0.8;

  Re = 1;
  run();

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

  Re = 80;
  run();

  Re = 100;
  run();

  Re = 200;
  run();

  Re = 300;
  run();
}

event init (t = 0) {
  solid (vof, sf, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq(D/2));
  refine (vof[] < 1 && vof[] > 0 && level < LEVEL);
  solid (vof, sf, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq(D/2));
}

/*
event moving_cylinder (i++) {
  solid (vof, sf, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq(D/2));
  // foreach()
  //  u.x[] = (1 - vof[])*u.x[];
}
*/

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}

event logfile (i++) {
  coord Fp, Fmu;

  immersed_force (vof, p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x) / (0.5*sq(U0)*(D));
  double CL = (Fp.y + Fmu.y) / (0.5*sq(U0)*(D));

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g %g\n",
          i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, Fp.x, Fp.y, Fmu.x, Fmu.y); // 11
}

int count = 0;
event frequency (i++) {
  if (t >= tf_start && Re >= 50) {
    char name[80];
    sprintf (name, "freq-7.5-%d.dat", Re);
    FILE * fp = fopen (name, "a");
    foreach_point(7.5000, 7.5000) {
      fprintf (fp, "%d %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], p[]);
    }
    fclose (fp);

    sprintf (name, "freq-10-%d.dat", Re);
    FILE * fp1 = fopen (name, "a");
    foreach_point(10.0000, 7.5000) {
      fprintf (fp1, "%d %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], p[]);
    }
    fclose (fp1);
    count++;
  }
}


event profile (t = t_end) {
  int k = 0;
  double delta = L0/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "vprofx1-%d", Re); // x = 4.8125
  FILE * fv = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (ci.x + (20*delta), i) {
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
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (ci.x, i) {
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
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (12, i) {
      if (vof[] > 0 && vof[] < 1)
        k = 2.;
      else if (vof[] == 1)
	k = 1.;
      else
	k = 0.;
      fprintf (fv2, "%d %g %g %g %g %g %g\n", k, x, y, u.x[], u.y[], p[], Delta); 
    }
  }
  fflush (fv2);
  fclose (fv2);

  sprintf (name, "vprofy1-%d", Re); // y = L0/2
  FILE * fv3 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (i, ci.y) {
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
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (i, ci.y+(delta*15)) {
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

  sprintf (name, "right-boundary-%d.dat", Re);
  FILE * fpv = fopen (name, "w");
  foreach_boundary(right)
    fprintf (fpv, "%g %g %g %g %g %g\n", x, y, Delta, u.x[], u.y[], p[]);
  fflush (fpv);
  fclose (fpv);

  scalar ring[];
  fraction (ring, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq((0.5+0.0146)/2));
  sprintf (name, "ring1-%d", Re);
  FILE * fv6 = fopen (name, "w");
  foreach()
    if (ring[] > 0 && ring[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv6, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv6);
  fclose (fv6);

  fraction (ring, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq((0.5+0.0293)/2));
  sprintf (name, "ring2-%d", Re);
  FILE * fv7 = fopen (name, "w");
  foreach()
    if (ring[] > 0 && ring[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv7, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv7);
  fclose (fv7);

  fraction (ring, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq((0.5+0.0439)/2));
  sprintf (name, "ring3-%d", Re);
  FILE * fv8 = fopen (name, "w");
  foreach()
    if (ring[] > 0 && ring[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv8, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv8);
  fclose (fv8);
}

event snapshot (t = t_end) {
  scalar omega[];
  vorticity (u, omega);

  char name[80];
  view (fov = 2, tx = -0.25, ty = -0.465,
        width = 3000, height = 1500); 
  sprintf (name, "%g-pressure-%d.png", t, Re);
  clear();
  draw_vof ("vof", "sf", lw = 5, lc = {0,0,0});
  squares ("p", min = -0.5, max = 0.5, map = blue_white_red);
  save (name);

  sprintf (name, "%g-vort-%d.png", t, Re);
  clear();
  squares ("omega", min = -3, max = 3, map = blue_white_red);
  draw_vof ("vof", "sf", lw = 5, lc = {0,0,0});
  save (name);

  sprintf (name, "%g-pressureiso-%d.png", t, Re);
  clear();
  draw_vof ("vof", "sf", lw = 5, lc = {0,0,0});
  isoline ("p", n = 20, min = -0.5, max = 0.5, lc = {0,0,0});
  squares ("p", min = -0.5, max = 0.5, map = blue_white_red);
  save (name);

  sprintf (name, "%g-vortiso-%d.png", t, Re);
  clear();
  isoline ("omega", n = 15, min = -3, max = 3, lc = {0,0,0});
  squares ("omega", min = -3, max = 3, map = blue_white_red);
  draw_vof ("vof", "sf", lw = 5,  lc = {0,0,0});
  save (name);
}

event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fprintf (fp, "%d\t%d\t%d\t%d\n", mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);
  fflush (fp);
  return 1;
}

