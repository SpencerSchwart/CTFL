#include "navier-stokes/centered.h"
// #include "../centered_immersed.h"
#include "../immersed.h" 
#include "view.h"

#define L0 20
#define D 0.5
#define LEVEL 11

int maxlevel = 11;
int Re;
double U0 =  1.; // inlet velocity
double t_end = 40;
double tf_start = 40;
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
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = neumann (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);


int main() {
  size(L0);
  init_grid (1 << (LEVEL-2));
  mu = muv;
  TOLERANCE = 1.e-6; 
  DT = 0.01;

  Re = 40;
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


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}

vector Uc[];
vector ub[];
double avgCD = 0;
double avgCL = 0;
int cdcount = 0;

event logfile (i++){
  coord desiredF = {0};
  coord markerF = {0};
  coord interpolateF = {0};
  coord gridF = {0};

  coord Fctotal = {0};
  immersed_forcev2 (vof, &desiredF, &markerF, &interpolateF, &gridF);
  // double CD = (Fp.x)/(1);
  // avgCD += t > tf_start? CD : 0;
  // double CL = (Fp.y)/(1);
  // avgCL += t > tf_start? CL : 0;
  // cdcount += t > tf_start? 1: 0;

  char name[80];
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
    }
    else
      foreach_dimension()
        Uc.x[] = 0;
        ub.x[] = 0;
    if(Fc.x[]) {
      foreach_dimension()
	Fctotal.x += Fd.x[];
    }
  }

  double area = interface_area (vof);

 double uavg_x = usum.x/counter;
 double uavg_y = usum.y/counter;
 fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, 
	   uavg_x, uavg_y, Fctotal.x, Fctotal.y, area,
       desiredF.x, desiredF.y, markerF.x, markerF.y, 
       interpolateF.x, interpolateF.y, gridF.x, gridF.y);


}

event profile (t = t_end) {
  int k = 0;
  double delta = L0/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "vprofx1-%d", Re); // x = 4.8125
  FILE * fv = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
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
  for(double i = 0; i <= L0; i += delta) {
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
  for(double i = 0; i <= L0; i += delta) {
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
      fprintf (fv5, "%g %g %g %g %g %g %g\n", theta, mag, cp, Fc.x[], Fc.y[], Fd.x[], Fd.y[]);
    }
  }
  fflush (fv5);
  fclose (fv5);

  scalar boundaryVelocity[];
  fraction (boundaryVelocity, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq(0.51758/2));
  sprintf (name, "surfaceVelocity1-%d", Re);
  FILE * fv6 = fopen (name, "w");
  foreach()
    if (boundaryVelocity[] > 0 && boundaryVelocity[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv6, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv6);
  fclose (fv6);

  fraction (boundaryVelocity, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq(0.53516/2));
  sprintf (name, "surfaceVelocity2-%d", Re);
  FILE * fv7 = fopen (name, "w");
  foreach()
    if (boundaryVelocity[] > 0 && boundaryVelocity[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv7, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv7);
  fclose (fv7);
}


event snapshot (t = t_end) {
  scalar omega[];
  vorticity (u, omega);

  char name[80];
  sprintf (name, "vort_final-%d.png", Re);
  view (fov = 2, tx = -0.275, ty = -0.50,
        width = 3000, height = 1500); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  save (name);

  sprintf (name, "xvelo-%d.png", Re);
  clear();
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  squares ("u.x", min = 0, max = 1.5, map = cool_warm);
  save (name);

  sprintf (name, "yvelo-%d.png", Re);
  clear();
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  squares ("u.y", min = -0.5, max = 0.5, map = cool_warm);
  save (name);

  sprintf (name, "pressure-%d.png", Re);
  clear();
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  squares ("u.y", min = -0.5, max = 0.5, map = cool_warm);
  save (name);

  sprintf (name, "vectors-%d.png", Re);
  clear();
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  vectors ("u", scale =  0.01);
  squares ("u.x", min = -0.5, max = 1, map = cool_warm);
  save (name);

  sprintf (name, "pressureiso-%d.png", Re);
  clear();
  view (fov = 1, tx = -0.275, ty = -0.50,
	    width = 3000, height = 1500);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  isoline ("p", n = 20, min = -0.5, max = 0.5, lc = {1,0.647,0});
  save (name);

  sprintf (name, "veloiso-%d.png", Re);
  clear();
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  isoline ("u.x", n = 25, min = -.1, max = 1.15, lc = {1,0.647,0});
  save (name);
}


event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}

int count = 0;
event frequency (i++) {
  if (t >= tf_start) {
    char name[80];
    sprintf (name, "freq-7.5.dat");
    FILE * fp = fopen (name, "a");
    foreach_point(7.500, ci.y) {
      fprintf (fp, "%d %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], p[]);
    }
    fclose (fp);

    sprintf (name, "freq-10.dat");
    FILE * fp1 = fopen (name, "a");
    foreach_point(10.000, ci.y) {
      fprintf (fp1, "%d %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], p[]);
    }
    fclose (fp1);
    count++;
  }
}

event dump (t += 10) {
  char name[80];

  sprintf (name, "%d-dump-%g", Re, t);
  dump (file = name);
}

/*
event movie (t += 0.1; t <= t_end) {
  scalar omega[];
  vorticity (u, omega);
  char name[80];
  sprintf (name, "vort-%g.png", t);
  FILE * fp1 = fopen (name, "w");
  view (fov = 4, tx = -.375, ty = -.488, bg = {1,1,1,},
        width = 1024, height = 512);
  clear();
  isoline ("omega", n = 15, min = -10, max = 10);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  squares ("omega", map = cool_warm, min = -10, max = 10);
  save (fp = fp1);
  fflush (fp1);
  fclose (fp1);
}
*/

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}
