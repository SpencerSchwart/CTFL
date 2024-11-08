#include "navier-stokes/centered.h"
#include "../immersed-new.h"
#include "view.h"

#define L0 29.257
#define D 1.
#define LEVEL 11

int maxlevel = 11;
int Re;
double U0 =  1.; // inlet velocity
double t_end = 200;
double tf_start = 100;
coord ci = {L0/5, L0/2}; // initial coordinates of cylinder
coord vc = {0, 0}; // velocity of cylinder
coord xc = {0, 0};
const double A = 0.2*D;
double freq = 0.156;

scalar vof[], * interfaces = {vof};

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

u.n[top]   = neumann (0);
u.n[bottom] = neumann (0);


int main() {
  size(L0);
  init_grid (1 << (LEVEL - 3));
  mu = muv;
  TOLERANCE = 1.e-7; 
  // DT = 0.004;
  CFL = 0.5;

  Re = 185;
  run();
}


event moving_cylinder (t = 0) {
  xc.y = -A*cos(2*M_PI*freq*t);
  vc.y = A*2*M_PI*freq*sin(2*M_PI*freq*(t)); // + dt or not?
  solid (vof, sf, - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2));
  refine (vof[] > 0 && vof[] < 1 && level < LEVEL);
  solid (vof, sf, - sq(x - ci.x - xc.x) - sq(y - ci.y - xc.y) + sq(D/2));
}

event init (t = 0) {
  refine (vof[] > 0 && vof[] < 1 && level < maxlevel);
}


event properties (i++) {
  xc.y = -A*cos(2*M_PI*freq*t);
  vc.y = A*2*M_PI*freq*sin(2*M_PI*freq*(t)); // + dt or not?

  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}

double avgCD = 0, avgCL = 0;
int count = 0;
scalar e[];
scalar p01[];
vector eu[];
vector u0[];
vector ef[];
vector f0[];

event logfile (i++, t <= t_end){

  coord Fu = {0};
  foreach() {
    e[] = p[] - p01[];
    foreach_dimension() {
      eu.x[] = u.x[] - u0.x[];
      Fu.x += u.x[] * vof[] * dv();
      ef.x[] = forceTotal.x[] - f0.x[];
    }
  }

  coord Fp;
  coord Fmu; 

  // immersed_force (vof, p, u, mu, &Fp, &Fmu);
  coord F = ibm_force();
  double CD = (F.x + Fu.x)/(0.5*sq(U0)*(D));
  avgCD += t > tf_start? CD: 0;
  double CL = (F.y + Fu.y)/(0.5*sq(U0)*(D));
  avgCL += t > tf_start? CL: 0;
  count += t > tf_start? 1:0;

  double volume = 0.;
  coord Fd = {0};
  foreach() {
    volume += vof[]*dv();
    p01[] = p[];
    foreach_dimension() {
        u0.x[] = u.x[];
        Fd.x += forceTotal.x[] * dv();
        f0.x[] = forceTotal.x[];
    }
  }
/*
  fprintf (stderr, "%d %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t, Re, CD, CD2, CD3, CD4, avgCD/(count + 1e-6), CL, CL2, CL3, CL4, avgCL/(count + 1e-6), // 13
       Fp.x, Fp.y, Fp2.x, Fp2.y, Fp3.x, Fp3.y, Fp4.x, Fp4.y,                                       // 21
       Fmu.x, Fmu.y, Fmu2.x, Fmu2.y, Fmu3.x, Fmu3.y, Fmu4.x, Fmu4.y,                               // 29
       xc.y, vc.y, Fd.x, Fd.y);                                                                    // 33
*/       

  fprintf (stderr, "%d %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
           i, t, Re, CD, avgCD/(count + 1e-6), CL, avgCL/(count + 1e-6), // 7
           Fp.x, Fp.y, Fmu.x, Fmu.y, xc.y, vc.y, // 13
           2.*Fd.x, 2.*Fd.y, Fu.x, Fu.y, volume); // 17
}



event profile (t += 195.51) {
  double delta = 20/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "surface-%d", Re);
  FILE * fv5 = fopen (name, "w");
  foreach() {
    if (vof[] > 0 && vof[] < 1) {
      coord pc, ucb, fcb;
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);
      double pcb1 = embed_interpolate (point, p, pc);
      double cp4 = extrapolate_scalar (point, vof, pc, n, p);
      double cp5 = extrapolate_scalarv2 (point, vof, pc, n, p);

      cp4 /= (0.5 * sq(U0));
      cp5 /= (0.5 * sq(U0));

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      bilinear_interpolation (point, u, pc, &ucb);
      // bilinear_interpolation (point, fc, pc, &fcb);
      double pcb = scalar_bilinear_interpolation (point, p, pc);

      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI) + 180;
      double umag = sqrt(sq(ucb.x) + sq(ucb.y));
      // double fmag = sqrt(sq(fcb.x) + sq(fcb.y));
      double cp = pcb / (0.5 * sq(U0));
      double cp1 = pcb1 / (0.5 * sq(U0));
      double omega = ibm_vorticity (point, u, pc, n);
      fprintf (fv5, "%g %g %g %g %g %g %g %g %g %g\n", xc.y, theta, umag, cp, cp1, omega, cp4, cp5);
    }
  }
  fflush (fv5);
  fclose (fv5);

  scalar omega[];
  vorticity (u, omega);

  sprintf (name, "vort_final-%d.png", Re);
  view (fov = 2, tx = -0.26, ty = -0.50,
        width = 3000, height = 1500); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  save (name);
}


event profile2 (t += 197.12) {
  double delta = 20/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "surface1-%d", Re);
  FILE * fv5 = fopen (name, "w");
  foreach() {
    if (vof[] > 0 && vof[] < 1) {
      coord pc, ucb, fcb;
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);
      double pcb1 = embed_interpolate (point, p, pc);
      double cp4 = extrapolate_scalar (point, vof, pc, n, p);
      double cp5 = extrapolate_scalarv2 (point, vof, pc, n, p);

      cp4 /= (0.5 * sq(U0));
      cp5 /= (0.5 * sq(U0));

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      bilinear_interpolation (point, u, pc, &ucb);
      // bilinear_interpolation (point, fc, pc, &fcb);
      double pcb = scalar_bilinear_interpolation (point, p, pc);

      double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI) + 180;
      double umag = sqrt(sq(ucb.x) + sq(ucb.y));
      // double fmag = sqrt(sq(fcb.x) + sq(fcb.y));
      double cp = pcb / (0.5 * sq(U0));
      double cp1 = pcb1 / (0.5 * sq(U0));
      double omega = ibm_vorticity (point, u, pc, n);
      fprintf (fv5, "%g %g %g %g %g %g %g %g\n", xc.y, theta, umag, cp, cp1, omega, cp4, cp5);
    }
  }
  fflush (fv5);
  fclose (fv5);

  scalar omega[];
  vorticity (u, omega);

  sprintf (name, "vort_final1-%d.png", Re);
  view (fov = 2, tx = -0.26, ty = -0.50,
        width = 3000, height = 1500); 
  isoline ("omega", n = 15, min = -3, max = 3);
  draw_vof ("vof", "sf", filled = 1, lw = 5);
  save (name);
}


/*
event movie (t += 0.1; t <= t_end) {
  scalar omega[];
  vorticity (u, omega);
  char name[80];
  sprintf (name, "vort-%g.png", t);
  FILE * fp1 = fopen (name, "w");
  view (fov = 2, tx = -.190, ty = -.5, bg = {1,1,1,},
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

event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = 2);
}


event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}
