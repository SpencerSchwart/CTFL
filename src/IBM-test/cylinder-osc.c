#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"
#include "../immersed.h" // IBM
#include "view.h"

int LEVEL = 9;
int maxlevel = 10;

double RE, KC, Amp, freq = 1 [*];
const double U0 =  1 [1,-1];
const double D = 0.5;
const double T_max = 30; // max # of periods

coord ci = {0}; // initial coordinates of cylinder
coord vc = {0}; // velocity of cylinder
coord xc = {0}; // displacement of cylinder

scalar vof[];
face vector sf[];
face vector muv[];

u.n[left] = dirichlet (0);
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = dirichlet (0);
u.t[right] = dirichlet (0);
p[right]   = neumann (0);
pf[right]  = neumann (0);

u.n[top] = neumann (0);
p[top] = dirichlet (0);
pf[top] = dirichlet (0);

u.n[bottom] = neumann (0);
p[bottom] = dirichlet (0);
pf[bottom] = dirichlet (0);

double SDF (double x, double y) {
   // return sqrt(sq(x - ci.x) + sq(y - ci.y)) - (D/2);
   return sq (x - ci.x - xc.x) + sq (y - ci.y - xc.y) - sq(D/2);
}

int j, T;
int main() {
  L0 = 12.8*D;
  ci.x = L0/2;
  ci.y = L0/2;
  size(L0);
  init_grid (2 << (LEVEL-4));
  mu = muv;
  TOLERANCE = 1.e-5 [*]; 
  DT = 5e-4;

  /*
  T = 0;
  j = 0;
  RE = 100;
  KC = 5;
  run();
*/
  T = 0;
  j = 1;
  RE = 200;
  KC = 10;
  L0 = 25.6*D;
  size(L0);
  LEVEL = 10;
  init_grid (2 << (LEVEL-4));
  ci.x = L0/2;
  ci.y = L0/2;
  run();
}

event init (t = 0) {
  Amp = KC*D/(2*M_PI);
  freq = U0/(D*KC);
  foreach()
    foreach_dimension()
      u.x[] = 0.;
}

scalar ref[];
face vector fref[];
event moving_cylinder (i++) {
  vc.y = -Amp*2*M_PI*freq*cos(2*M_PI*freq*t);
  xc.y = -Amp*sin(2*M_PI*freq*t);
  solid (vof, sf, - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2));
  // solid (ref, fref, - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2));
/*
  coord nv;
  foreach() {
    normal_vector (point, vof, &nv);
    double lambda = fabs(nv.x) + fabs(nv.y);
    if (lambda != 0.) {
      double eta = 0.065*(1 - sq(lambda)) + 0.39;
      double num = SDF (x,y); 
      double den = lambda * eta * sqrt(2) * Delta;
      vof[] = den != 0? 0.5*(1 - tanh (num/den)): vof[];
    }
  }
  */
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(RE);
   boundary ((scalar *) {muv});
}

event period (t += 1/freq) {
  T++;
}

vector Uc[];
event logfile (i++; T <= T_max){
  coord Fp = {0};
  immersed_forcev2 (vof, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*D);
  double CL = (Fp.y)/(0.5*sq(U0)*D);
  
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

      bilinear_interpolation (point, u, pc, &ucb);

      foreach_dimension() {
	Uc.x[] = ucb.x;
	usum.x += ucb.x;
      }
    }
    else {
      foreach_dimension()
	Uc.x[] = 0;
    }
  }
 double uavg_x = usum.x/counter;
 double uavg_y = usum.y/counter;

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g %g %g %d\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, xc.y, vc.y, uavg_x, uavg_y, dt, T);
}

event movie (t += 0.05) {
  int FOV = j > 0? 5:10;
  scalar omega[];
  vorticity (u, omega);
  char name[80];
  sprintf (name, "%d-snap-%g.png", j, t);
  FILE * fp1 = fopen(name, "w");
  view (fov = FOV, tx = -0.5, ty = -0.5, width = 600, height = 675);
  draw_vof ("vof", "sf", lw = 3);
  squares ("omega", map = cool_warm);
  cells();
  save (fp = fp1);
  fclose (fp1);
}

event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}

