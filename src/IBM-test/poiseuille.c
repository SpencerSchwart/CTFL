#include "navier-stokes/centered.h"
#include "../immersed.h"
#include "view.h"

#define L0 20
int maxlevel = 7;

const double U0 = 1;
double t_end = 100;
double H = 5;
coord vc = {0, 0};
const double delta = 20/pow(2,7);
const double int_frac = 1;
double y_min = 2.5;
double y_max = 2.73;
double tt = delta * (1 - int_frac);

scalar vof[];
face vector sf[];
face vector muv[];

u.n[left] = neumann (0);
p[left]   = dirichlet (1);
pf[left]  = dirichlet (1);

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


int main() {
  size(L0);
  init_grid (1 << (maxlevel));
  origin (0, -(L0/2));
  mu = fm;
  TOLERANCE = 1.e-5;

  run();
}


event moving_wall (i++) {
  solid (vof, sf, fabs(y) > (H/2));
  foreach()
    vof[] = fabs(y) < y_max && fabs(y) > y_min? int_frac: vof[];
}

vector uref[];
event init (t = 0) {
  double U_max = (sq(H)*1)/(8*1*L0);
  foreach() {
    u.x[] = (1 - vof[])*U0;
    uref.x[] = fabs(y) <= y_min? U_max*(1 - sq((2*y/(H + 2*tt)))): 0;
  }
}

event profile (t = t_end) {
  char name[80];

  sprintf (name, "vprof-%g", int_frac);
  FILE * fv = fopen (name, "w");
  for (double i = -L0/2; i <= L0/2; i += delta)
    foreach_point (L0 - (delta), i) {
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
      } else if (vof[] == 1 && (vof[0,1] == 0 || vof[0,-1] == 0)){
        double xc = x;
        double yc = y + delta/2;
        coord pc = {xc, yc}, ucb;
        bilinear_interpolation (point, u, pc, &ucb);
        fprintf (fv, "%g %g %g\n", pc.x, pc.y, ucb.x);
    }
        fprintf (fv, "%g %g %g\n", x, y, u.x[]);
    }
  fclose (fv);

  sprintf (name, "vprofref-%g", int_frac);
  FILE * fref = fopen (name, "w");
  for (double i = -L0/2; i <= L0/2; i += delta)
    foreach_point (L0 - (delta), i)
      fprintf (fref, "%g %g %g\n", x, y, uref.x[]);
  fclose (fref);
}

scalar e[];
event logfile (i++, t <= t_end) {
  trash ({e});
  int counter = 0;
  coord usum = {0};
  coord pe = {0};
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
  }
  norm n = normf (e);

  double uavg_x = counter > 0? usum.x/counter: 0;
  double uavg_y = counter > 0? usum.y/counter: 0;

  double U_max = (sq(H)*1)/(8*1*L0);
  foreach() {
    double uref = U_max*(1 - sq((2*y/(H + 2*tt))));
    e[] = vof[] <= int_frac? fabs(u.x[] - uref): 0;
  }
  fprintf (stderr, "%d %g %g %g %g %g %g\n", i, t, uavg_x, uavg_y, n.avg, n.rms, n.max);
}

