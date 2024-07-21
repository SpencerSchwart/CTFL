@if _XOPEN_SOURCE < 700
  @undef _XOPEN_SOURCE
  @define _XOPEN_SOURCE 700
@endif
@if _GNU_SOURCE
@include <stdint.h>
@include <string.h>
@include <fenv.h>
@endif
#define _CATCH
#define dimension 2
#define BGHOSTS 2
#include "common.h"
#include "grid/quadtree.h"
#ifndef BASILISK_HEADER_0
#define BASILISK_HEADER_0
#line 1 "pitching.c"
#include "embed.h"
#undef EMBED
#include "navier-stokes/centered.h"
// #include "navier-stokes/double-projection.h"
#include "../immersed.h" // IBM
#include "../interfaceforce.h" // for CD and CL
#include "curvature.h"
#include "view.h"

#define MAX_THICKNESS 0.15
#define CHORD_LENGTH 1
#define L0 20.
#define Re (10000.)
#define LEVEL 12

double U0 =  1.0; // inlet velocity
// NACA XXXX airfoil properties
double rr = 1.1019*sq(MAX_THICKNESS); // Radius of leading edge
double mm = 0.02; // maximum camber
double pp = 0.4; // % pos of max camber
double p_ts = (0.2); // start time of pitching
double p_w0 = 0.6 [0,-1]; // angular velocity
double theta_i = 5; // starting angle

double t_end = 83 [0,1];
coord vc = {0.,0.}; // the velocity of the cylinder
coord ci = {5, 10}; // initial coordinates of airfoil
coord cr = {0.25*(CHORD_LENGTH), 0.}; // center of rotation

int j;

scalar airfoil[];
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

#define airfoil_thickness(x) (5 * MAX_THICKNESS * ((0.2969*(sqrt(x)))	\
						   -(0.1260*x)		\
						   -(0.3516*(sq(x)))	\
						   +(0.2843*(cube(x)))	\
						   -(0.1036*(pow(x, 4.)))))


void airfoil_shape (scalar c, face vector f, double theta)
{
  double yt, yc = 0;
  vertex scalar phi[];
  
  double chord = CHORD_LENGTH;
  
  foreach_vertex() {
    
    double XX = cr.x + (x - ci.x)*cos (theta) - (y - ci.y)*sin (theta);
    double YY = cr.y + (x - ci.x)*sin (theta) + (y - ci.y)*cos (theta);

    if (XX < 0) {
      // leading edge approximation
      phi[] =  - sq(XX-rr) - sq(YY) + sq(rr);
    }
    else if (XX >= 0. && XX <= chord) {
      // basic airfoil thickness
      yt = airfoil_thickness(XX);
      if (XX < pp*chord) {
	// airfoil camber (front)
	yc = (mm/sq(pp))*((2*pp*XX) - sq(XX));
      }
      else {
        // airfoil camber (back)
        yc = (mm/sq(1-pp))*((1-(2*pp))+(2*pp*XX)-sq(XX));
        }
      phi[] = YY > yc? - sq(YY) + sq(yc + yt): - sq(YY) + sq(yc - yt);
    }
    else
      phi[] = -1.;
  } 
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f);
  c.refine = c.prolongation = fraction_refine;
}

double SDF (double xx, double yy, double theta_p) {
  double XX = cr.x + (xx - ci.x)*cos (theta_p) - (yy - ci.y)*sin (theta_p);
  double YY = cr.y + (xx - ci.x)*sin (theta_p) + (yy - ci.y)*cos (theta_p);

  if (XX < 0) {
      // leading edge approximation
      return - sq(XX-rr) - sq(YY) + sq(rr);
    }
  else if (XX >= 0. && XX <= CHORD_LENGTH) {
      // basic airfoil thickness
      double yt = airfoil_thickness(XX);
      return YY > 0? - sq(YY) + sq(yt): - sq(YY) + sq(yt);
    }
  else
      return -1.;
}	


int main(){
  size(L0);
  init_grid (2 << (LEVEL-4));
  mu = muv;
  TOLERANCE = 1.e-5 [*];

  j = 0; // SPM
  run();

  j = 1; // No SPM
  run();
}


event init (t = 0) {
  // initial mesh refinement 
  astats ss;
  int ic = 0;
  do {
    ic++;
    airfoil_shape (airfoil, sf, theta_i*pi/180);
    ss = adapt_wavelet ({airfoil}, (double[]) {1.e-30},
			maxlevel = LEVEL, minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
  airfoil_shape(airfoil, sf, theta_i*pi/180);
}

scalar ref[];
face vector rf[];
vertex scalar phia[];
vector nv[];

double theta_p [0, 1];
event moving_cylinder (i++) {
  double time = t;
  theta_p = t >= p_ts? p_w0*(1- exp(-4.6*time)): 0;
  theta_p *= pi/180;
  airfoil_shape (airfoil, sf, theta_p);
  airfoil_shape (ref, rf, theta_p);
 
 /* 
  if (j == 0) { // SPM
  normal_vector(ref, nv);

  foreach() {
      double lambda = fabs(nv.x[]) + fabs(nv.y[]);
      double eta = 0.065*(1 - sq(lambda)) + 0.39;
      double num = - SDF (x,y);
      double den = lambda * eta * sqrt(2) * Delta;
      double scale = 12;
      double vof = airfoil[];
      airfoil[] = den != 0? 0.5*(1 - tanh ((num/den)*scale)):airfoil[];
      // fprintf (stderr, "i=%d t=%g n.x=%g n.y=%g lam=%g eta=%g vof_b=%g vof_a=%g x=%g y=%g delta=%g num=%g den=%g tanh(%g)\n",
// 		      i, t, nv.x[], nv.y[], lambda, eta, vof, airfoil[], x, y, Delta, num, den, (num/den)*scale);
    }
  }
  */
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(CHORD_LENGTH)/(Re);
  boundary ((scalar *) {muv});
}


event logfile (i++; t <= t_end) {

  coord Fp, Fmu;
  interface_force (airfoil, p, u, mu, &Fp, &Fmu);
  double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(CHORD_LENGTH));
  double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(CHORD_LENGTH));

  double E = 0;
  double E_p = 0;
  boundary ({u.x, u.y});
  scalar omega[];
  vorticity (u , omega);
  foreach(){
    double vort = omega[];
    double area = dv();
    E_p += sq(vort);
    if (airfoil[] < 1. && airfoil[] > 0){
      coord b, n;
      area *= embed_geometryo (point, &b, &n);
      vort = embed_vorticityo (point, u, b, n);
    }
    E += area*sq(vort);
  }

  
  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, E, E_p, theta_p*(180/pi));
  
}

/*
event movie (t += 1e-2; t <= 10)
{
  scalar omega[];
  vorticity (u, omega);
  view (fov = 17, camera = "front",
	tx = -0.5, ty = -0.5, bg = {1,1,1},
	width = 3200, height = 3200);
  clear();
  draw_vof ("airfoil", "sf",filled = 1, lw = 3);
  squares ("u.x", map = cool_warm);
  //squares ("omega", map = cool_warm);;
  save ("vorticity.mp4");
}
*/

event movie (t += 0.01; t <= t_end) {

}

event screenshot (t += 0.5; t <= t_end) {
  view (fov =2, camera = "front", 
        tx = -0.25, ty = -0.5, bg = {1,1,1},
	width = 3200, height = 3200);
  clear();
  draw_vof ("ref", "rf", filled = 1, lw = 3);
  isoline ("p", n = 120, min = -1, max = 1);
  if (j == 0)
    save("piso-0.png");
  else // j == 1
    save("piso-1.png");
}

event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = LEVEL - 6);
}


#endif
