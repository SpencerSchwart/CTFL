#include "navier-stokes/centered.h"
#include "../immersed.h"
#include "view.h"

#define L0 10
int maxlevel = 6;

const double U0 = 1;
double t_end = 100;
double H = 5;
const double delta = 10/pow(2,6);
double y_min = -0.23;
double y_max = 0;
const double int_frac = 0.25;
double tt = delta * (1 - int_frac);


coord vc = {0, 0};

scalar vof[];
face vector sf[];


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
  TOLERANCE = 1.e-5;

  run();
}


event moving_wall (i++) {
  solid (vof, sf, y < y_max);
  foreach() 
    vof[] = y < y_max && y > y_min? int_frac: vof[];
}

vector uref[];
event init (t = 0) {
  foreach() {
    u.x[] = y > y_max? y: 0;
    uref.x[] = U0*(y + tt)/(H + tt);
  }
}



event profile (t = t_end) {
  char name[80];

  sprintf (name, "vprof-%g", int_frac);
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
        fprintf (fv, "%g %g %g\n", pc.x, pc.y, ucb.x);
      } else if (vof[] == 1 && vof[0,1] == 0){
        double xc = x;
        double yc = y + delta/2;
        coord pc = {xc, yc}, ucb;
        bilinear_interpolation (point, u, pc, &ucb);
        fprintf (fv, "%g %g %g\n", pc.x, pc.y, ucb.x);
    }
        fprintf (fv, "%g %g %g\n", x, y, u.x[]);
    }
  fclose (fv);
}

scalar e[];
scalar perr[];
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
    e[] = vof[] <= int_frac? fabs(u.x[] - ref): 0;
    perr[] = vof[] <= int_frac? fabs(u.x[] - ref) / ref: 0;
    if (fabs(e[]) > max) {
      max = e[];
      pe.x = x;
      pe.y = y;
    }
  }
  norm n = normf (e);
  norm np = normf (perr);

  double uavg_x = counter > 0? usum.x/counter: 0;
  double uavg_y = counter > 0? usum.y/counter: 0;
  fprintf (stderr, "%d %g %g %g %g %g %g %g %g %g %g %g\n", i, t, uavg_x, uavg_y, n.avg, n.rms, n.max, max, pe.x, pe.y, np.avg, np.max);
}

/*
event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = maxlevel, minlevel = maxlevel - 2);
}
*/

