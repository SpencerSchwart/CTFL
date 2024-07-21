#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"

#define L0 10
int maxlevel = 8;

const double U0 = 1;
double Re = 50;
double t_end = 40;
double H = 5;
coord vc = {U0, 0};

face vector muv[];

// u.n[left] = dirichlet (y > 5? U0: 0);
// u.n[left] = dirichlet (0);
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = dirichlet (0);
// u.t[top] = dirichlet (U0);
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
  init_grid (1 << (maxlevel));
  mu = muv;
  TOLERANCE = 1.e-5;

  run();
}


event moving_wall (t = 0) {
  solid (cs, fs, y < (H));
}

/*
event init (t = 0) {
  foreach()
    u.x[] = cs[]*1;
}
*/

event propertires (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(H)/(Re);
  boundary ((scalar *) {muv});
}

event logfile (i++, t <= t_end) {
  fprintf (stderr, "%d %g\n", i, t);
}

event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = maxlevel, minlevel = 2);
}


