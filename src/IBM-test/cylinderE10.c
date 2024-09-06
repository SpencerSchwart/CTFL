#include "../my-embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define L0 15
#define D 0.5
#define LEVEL 11

int Re;
double U0 =  1.; // inlet velocity
double t_end = 50 [0,1];
double tf_start = 20 [0,1];
double xi = 4.8266;
double yi = 3.0835;
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
  DT = 0.01;

  Re = 40;
  run();

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


event init (t = 0) {
  // mask(y > 6 ? top: y < -6 ? bottom : none);
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

 /* 
  double E = 0;
  double E_p = 0;
  boundary ({u.x, u.y});
  scalar omega[];
  vorticity (u , omega);
  foreach(){
    double vort = omega[];
    double area = dv();
    if (cs[] < 1. && cs[] > 0){
      coord b, n;
      area *= embed_geometry (point, &b, &n);
      vort = embed_vorticity (point, u, b, n);
    }
    E += area*sq(vort);
    E_p += sq(vort);

  }
*/
  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g %g\n",
	   i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, Fp.x, Fp.y, Fmu.x, Fmu.y);
	   
}



event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}

event profile (t = t_end) {
  int k = 0;
  double delta = L0/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "vprofx1-%d", Re); // x = 4.8125
  FILE * fv = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (ci.x+(delta*20), i) {
      if (cs[] > 0 && cs[] < 1)
        k = 2.;
      else if (cs[] == 1)
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
      if (cs[] > 0 && cs[] < 1)
        k = 2.;
      else if (cs[] == 1)
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
      if (cs[] > 0 && cs[] < 1)
        k = 2.;
      else if (cs[] == 1)
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
      if (cs[] > 0 && cs[] < 1)
        k = 2.;
      else if (cs[] == 1)
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
      if (cs[] > 0 && cs[] < 1)
        k = 2.;
      else if (cs[] == 1)
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
  scalar omega[];

  foreach() {
    if (cs[] > 0. && cs[] < 1.) {
      coord b, n;
      embed_geometry (point, &b, &n);
      double xe = x + b.x*Delta, ye = y + b.y*Delta;

      omega[] = embed_vorticity (point, u, b, n);
      coord dudn = embed_gradient (point, u, b, n);
      double dudnMag = sqrt(sq(dudn.x) + sq(dudn.y));

      double theta = atan2 ((ye - ci.y), (xe - ci.x)) * (180/M_PI);
      double cp = embed_interpolate (point, p, b) / (0.5 * sq(U0));

      coord avgu = {0};
      int counter = 0;
      foreach_neighbor()
        if (cs[] == 1) {
            foreach_dimension()
                avgu.x += u.x[];
            counter++;
        }

      fprintf (fv5, "%g %g %g %g %g %g %g %g %g %g %g %g\n", theta, cp, // 2
                     pressureDrag.x[], pressureDrag.y[], frictionDrag.x[], frictionDrag.y[], // 6
                     omega[], dudn.x, dudn.y, dudnMag, avgu.x/counter, avgu.y/counter); // 12
    }                     
  }
  fflush (fv5);
  fclose (fv5);

  scalar boundaryVelocity[];
  fraction (boundaryVelocity, - sq(x - ci.x) - sq(y - ci.y) + sq(0.51758/2));
  sprintf (name, "surfaceVelocity1-%d", Re);
  FILE * fv6 = fopen (name, "w");
  foreach()
    if (boundaryVelocity[] > 0 && boundaryVelocity[] < 1) {
      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);
      fprintf (fv6, "%g %g %g %g %g %g\n", x, y, theta, u.x[], u.y[], p[]);
    }
  fflush (fv6);
  fclose (fv6);

  fraction (boundaryVelocity, - sq(x - ci.x) - sq(y - ci.y) + sq(0.53516/2));
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
  view (fov = 2, tx = -0.275, ty = -0.50,
	width = 3000, height = 1500); 
  sprintf (name, "xvelo-%d.png", Re);
  clear();
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.x", min = 0, max = 1.5, map = cool_warm);
  save (name);

  sprintf (name, "yvelo-%d.png", Re);
  clear();
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", min = -0.5, max = 0.5, map = cool_warm);
  save (name);

  sprintf (name, "pressure-%d.png", Re);
  clear();
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  squares ("u.y", min = -0.5, max = 0.5, map = cool_warm);
  save (name);

  sprintf (name, "vectors-%d.png", Re);
  clear();
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  vectors ("u", scale =  0.01);
  squares ("u.x", min = -0.5, max = 1, map = cool_warm);
  save (name);

  sprintf (name, "pressureiso-%d.png", Re);
  clear();
  view (fov = 1, tx = -0.275, ty = -0.50,
	    width = 3000, height = 1500);
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  isoline ("p", n = 20, min = -0.5, max = 0.5);
  save (name);

  sprintf (name, "veloiso-%d.png", Re);
  clear();
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  isoline ("u.x", n = 25, min = -.1, max = 1.15);
  save (name);

  sprintf (name, "vortiso-%d.png", Re);
  clear();
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("cs", "fs", filled = -1, lw = 5);
  save (name);
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


event movie (t += 0.01; t <= t_end) {
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}
