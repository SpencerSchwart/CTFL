#include "navier-stokes/centered.h"
#include "../immersed.h"
#include "view.h"

#define L0 10
int maxlevel = 7;

const double U0 = 5;
double Re = 50;
double t_end = 100;
double H = 5;
const double delta = 10/pow(2,7);
double tt = delta * (1.75);
double y_min = -0.15;
double y_max = -0.10;

coord vc = {0, 0};

scalar vof[];
face vector sf[];


/*
u.n[left] = dirichlet (y < y_min? 0: y);
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);
*/

u.n[top] = dirichlet (0);
u.t[top] = dirichlet (U0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = dirichlet (0);
u.t[bottom] = dirichlet (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);


int main() {
  size(L0);
  init_grid (1 << (maxlevel));
  origin (0, -5);
  mu = fm;
  periodic (right);
  // stokes = true;
  TOLERANCE = 1.e-5;

  run();
}


event moving_wall (i++) {
  solid (vof, sf, y < y_max);
  foreach() 
    vof[] = y < y_max && y > y_min? 0.25: vof[];
}

vector uref[];
event init (t = 0) {
  foreach() {
    u.x[] = y > 0? y: 0;
    uref.x[] = U0*(y + tt)/(H + tt);
  }
}



event profile (t = t_end) {
  char name[80];

  sprintf (name, "vprof3");
  FILE * fv = fopen (name, "w");
  for (double i = -L0/2; i <= L0/2; i += delta)
    foreach_point (0.5*L0, i) {
      coord ucb = {0};
      if (vof[] > 0 && vof[] < 1) {
        coord pc, ucb;
        coord m = {x, y};
        coord n = interface_normal (point, vof);
        double alpha = plane_alpha (vof[], n);
        double area = plane_area_center(n, alpha, &pc);

        foreach_dimension()
          pc.x = m.x + pc.x*Delta; 

        bilinear_interpolation (point, u, pc, &ucb);
      } else if (vof[] == 1 && vof[0,1] == 0){
        double xc = x;
        double yc = y + delta/2;
        coord pc = {xc, yc}, ucb;
        bilinear_interpolation (point, u, pc, &ucb);
    }
        fprintf (fv, "%g %g %g %g\n", x, y, u.x[], ucb.x);
    }
  fclose (fv);
}

scalar e[];
event logfile (i++, t <= t_end) {
  trash ({e});
  int counter = 0;
  coord usum = {0};
  coord pe = {0};
  double max = -HUGE;
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

      foreach_dimension() 
        usum.x += ucb.x;

    } else if (vof[] == 1 && vof[0,1] == 0){
      counter++;
      double xc = x;
      double yc = y + delta/2;
      coord pc = {xc, yc}, ucb;
      bilinear_interpolation (point, u, pc, &ucb);

      foreach_dimension() 
        usum.x += ucb.x;
    }
    double ref = U0*(y + tt)/(H + tt);
    e[] = vof[] < 1? u.x[] - ref: 0;
    if (fabs(e[]) > max) {
      max = e[];
      pe.x = x;
      pe.y = y;
    }
  }
  norm n = normf (e);

  double uavg_x = counter > 0? usum.x/counter: 0;
  double uavg_y = counter > 0? usum.y/counter: 0;
  fprintf (stderr, "%d %g %g %g %g %g %g %g %g %g\n", i, t, uavg_x, uavg_y, n.avg, n.rms, n.max, max, pe.x, pe.y);
}

/*
event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = maxlevel, minlevel = maxlevel - 2);
}
*/

