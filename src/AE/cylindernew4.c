#include "navier-stokes/centered.h"
#include "../immersed-new.h" 
#include "view.h"

#define L0 25.
#define D 1.
#define LEVEL 9
int maxlevel = 9;

int Re;
double U0 =  1.; // inlet velocity
double t_end = 20;
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

int main() {
  size(L0);
  init_grid (1 << (LEVEL-2));
  mu = muv;
  TOLERANCE = 1.e-6; 
  CFL = 0.8;

  Re = 40;
  run();

}


event moving_cylinder (i++) {
  solid (vof, sf, - sq(x - ci.x) - sq(y - ci.y) + sq(D/2));
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   //boundary ((scalar *) {muv});
}


event logfile (i++) {
  coord F = ibm_force();
  double CD = (-F.x)/(0.5*sq(U0)*(D));
  double CL = (-F.y)/(0.5*sq(U0)*(D));
  
  fprintf (stderr, "%d %g %d %d %d %d %d %g %g\n",
          i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL);
}

event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}

event snapshot (t = t_end) {
  FILE * file = fopen("data.txt", "w");
  output_field (list = {u.x, u.y, p}, fp = file, n = 512);
  fclose(file);
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}
