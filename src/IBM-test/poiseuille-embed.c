#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define L0 25
int maxlevel = 8;

const double U0 = 1;
double Re = 50;
double t_end = 50;
double H = 2.5;
coord vc = {0, 0};

scalar vof[];
face vector sf[];
face vector muv[];

u.n[left] = dirichlet (fabs(y) <= H? U0: 0);
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

u.t[embed] = dirichlet (0);
u.n[embed] = dirichlet (0);


int main() {
  size(L0);
  init_grid (2 << (8));
  origin (0, -(L0/2));
  mu = muv;
  TOLERANCE = 1.e-5;

  DT = 0.01;

  run();
}


event moving_wall (t = 0) {
  solid (cs, fs, fabs(y) < (H));
}


event propertires (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(H)/(Re);
  boundary ((scalar *) {muv});
}

event profile (t = t_end) {
  double delta = L0/pow(2,maxlevel);
  char name[80];

  sprintf (name, "vprof");
  FILE * fv = fopen (name, "w");
  for (double i = -L0/2; i <= L0/2; i += delta)
    foreach_point (L0-(3*delta), i)
      fprintf (fv, "%g %g %g\n", x, y, u.x[]);
  fclose (fv);
}

event logfile (i++, t <= t_end) {
  fprintf (stderr, "%d %g\n", i, t);
}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = maxlevel, minlevel = maxlevel - 4);
}


