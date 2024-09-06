#include "navier-stokes/centered.h"
#include "../immersed-new.h" 
#include "view.h"

#define L0 15
#define D 0.5
#define LEVEL 10

int maxlevel = 10;
int Re;
double U0 =  1.; // inlet velocity
double t_end = 50;
double tf_start = 25;
coord ci = {L0/4, L0/2}; // initial coordinates of cylinder
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

u.n[bottom] = neumann (0);

int cdcount = 0;
double avgCD = 0, avgCL = 0;
int count = 0;

int main() {
  size(L0);
  init_grid (1 << (LEVEL-2));
  mu = muv;
  TOLERANCE = 1.e-6; 
  DT = 0.01;

  Re = 20;
  run();

  avgCD = 0, avgCL = 0;
  count = 0, cdcount = 0;

  Re = 40;
  run();

  avgCD = 0, avgCL = 0;
  count = 0, cdcount = 0;
  
  Re = 100;
  run();
  
  avgCD = 0, avgCL = 0;
  count = 0, cdcount = 0;

  Re = 185;
  run();

  avgCD = 0, avgCL = 0;
  count = 0, cdcount = 0;

  Re = 200;
  run();

}


event moving_cylinder (i++) {
  solid (vof, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}


event logfile (i++) {
  coord Fp, Fp2, Fp3, Fp4;
  coord Fmu, Fmu2, Fmu3, Fmu4;
  immersed_forcev3 (vof, p, u, mu, &Fp, &Fmu); 
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(D));
  avgCD += t > tf_start? CD : 0;
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(D));
  avgCL += t > tf_start? CL : 0;
  cdcount += t > tf_start? 1: 0;

  coord sumF = {0};
  foreach()
    foreach_dimension()
        sumF.x += forceTotal.x[]*dv();
  
  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g\n",
          i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, // 7 
          CD, avgCD/(cdcount + 1.e-6) , CL, avgCL/(cdcount + 1.e-6), // 11
          Fp.x, Fp.y, Fmu.x, Fmu.y, sumF.x/(0.5*sq(U0)*(D)), sumF.y/(0.5*sq(U0)*(D))); // 17
  
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

  sprintf (name, "left-boundary-%d.dat", Re);
  FILE * fpv = fopen (name, "w");
  foreach_boundary(right)
    fprintf (fpv, "%g %g %g %g %g %g\n", x, y, Delta, u.x[], u.y[], p[]);
  fflush (fpv);
  fclose (fpv);

  sprintf (name, "surface-%d", Re);
  FILE * fv5 = fopen (name, "w");
  foreach() {
    if (vof[] > 0 && vof[] < 1) {
      coord pc, ucb;
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);
      coord dudne = embed_gradient (point, u, pc, n, vc);
      double omegaE = embed_vorticity (point, u, pc, n);
      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      bilinear_interpolation (point, u, pc, &ucb);
      double pcb = scalar_bilinear_interpolation (point, p, pc);
      n = facet_normal (point, vof, sf);
      normalize (&n);
      double omega = ibm_vorticity (point, u, pc, n);
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      double mag = sqrt(sq(ucb.x) + sq(ucb.y));
      double cp = pcb / (0.5 * sq(U0));

      coord dudn = ibm_gradientv2 (point, u, pc, n);
      double magdudn = distance (dudn.x, dudn.y);
      
      double Fp1 = -(totalDesired.x[]*n.x + totalDesired.y[]*n.y) * area;
      double Fp2 = -(forceTotal.x[]*n.x + forceTotal.y[]*n.y) * area;

      coord avgu = {0};
      int counter = 0;
      foreach_neighbor()
        if (vof[] == 0) {
            foreach_dimension()
                avgu.x += u.x[];
            counter++;
        }

      fprintf (fv5, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", 
                    theta, mag, cp, desiredForce.x[], desiredForce.y[], // 5
                    pressureDrag.x[], pressureDrag.y[], frictionDrag.x[], frictionDrag.y[], // 9
                    omega,dudn.x, dudn.y, magdudn, avgu.x/counter, avgu.y/counter, // 15
                    omegaE, dudne.x, dudne.y, Fp1, Fp2); // 20
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
  view (fov = 2, tx = -0.235, ty = -0.465,
        width = 3000, height = 1500); 
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
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  isoline ("p", n = 20, min = -0.5, max = 0.5, lc = {1,0,0});
  save (name);

  sprintf (name, "veloiso-%d.png", Re);
  clear();
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  isoline ("u.x", n = 25, min = -.1, max = 1.15, lc = {1,0,0});
  save (name);

  sprintf (name, "vortiso-%d.png", Re);
  clear();
  isoline ("omega", n = 15, min = -3, max = 3, lc = {1,0,0});
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  save (name);

}


event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}


event frequency (i++) {
  if (t >= tf_start && Re >= 100) {
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

event dump (t += 20) {
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
