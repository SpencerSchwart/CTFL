#include "navier-stokes/centered.h"
#include "../immersed.h" // IBM
#include "view.h"

#define MAX_THICKNESS 0.15
#define CHORD_LENGTH 1
#define L0 20.
#define Re (10000.)
#define LEVEL 12
int maxlevel = 12;

double U0 =  1.0; // inlet velocity
// NACA XXXX airfoil properties
double rr = 1.1019*sq(MAX_THICKNESS); // Radius of leading edge
double mm = 0.02; // maximum camber
double pp = 0.4; // % pos of max camber
double p_ts = (10); // start time of pitching
double p_w0 = 0.6 [0,-1]; // angular velocity
double theta_i = 45; // starting angle

double t_end = 83 [0,1];
coord vc = {0.,0.}; // the velocity of the cylinder
coord ci = {5, 10}; // initial coordinates of airfoil
coord cr = {0.25*(CHORD_LENGTH), 0.}; // center of rotation

int j;

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

u.n[top]   = neumann (0);
u.n[bottom] = neumann (0);

#define airfoil_thickness(x) (5 * MAX_THICKNESS * ((0.2969*(sqrt(x)))	\
						   -(0.1260*x)		\
						   -(0.3516*(sq(x)))	\
						   +(0.2843*(cube(x)))	\
						   -(0.1036*(pow(x, 4.)))))


/*
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
*/

void airfoil_shape (scalar c, face vector f, double theta)
{
  double yt, yc = 0;
  double chord = CHORD_LENGTH;
  vertex scalar phi[]; 
  foreach_vertex() {
    
    double XX = cr.x + (x - ci.x)*cos (theta) - (y - ci.y)*sin (theta);
    double YY = cr.y + (x - ci.x)*sin (theta) + (y - ci.y)*cos (theta);

    if (XX < 0) {
      // leading edge approximation
      phi[] = - sq(XX-rr) - sq(YY) + sq(rr);
    }
    else if (XX >= 0. && XX <= chord) {
      // basic airfoil thickness
      yt = airfoil_thickness(XX);
      
      phi[] = YY > yc? - sq(YY) + sq(yt): - sq(YY) + sq(yt);
    }
    else
      phi[] = -1.;
  } 
  boundary ({phi});
  fractions (phi, c, f);
  // fractions_cleanup (c, f);
  c.refine = c.prolongation = fraction_refine;
}


int main(){
  size(L0);
  init_grid (1 << (LEVEL-4));
  mu = muv;
  TOLERANCE = 1.e-5 [*];

  run();
}


event init (t = 0) {
  // initial mesh refinement 
  astats ss;
  int ic = 0;
  do {
    ic++;
    airfoil_shape (vof, sf, theta_i*pi/180);
    ss = adapt_wavelet ({vof}, (double[]) {1.e-30},
			maxlevel = LEVEL, minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
  airfoil_shape(vof, sf, theta_i*pi/180);
}


double theta_p;
event moving_cylinder (i++) {
  theta_p = t >= p_ts? p_w0*t + theta_i: theta_i;
  theta_p *= pi/180;
  airfoil_shape (vof, sf, theta_p);
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(CHORD_LENGTH)/(Re);
  boundary ((scalar *) {muv});
}


event logfile (i++; t <= t_end) {
  coord Fp = {0}; 
  immersed_forcev2 (vof, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*(CHORD_LENGTH));
  double CL = (Fp.y)/(0.5*sq(U0)*(CHORD_LENGTH));
  
  fprintf (stderr, "%d %g %d %d %d %d %d %g %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, theta_p*(180/pi));
  
}


event screenshot (t += 1; t <= t_end) {
  char name[80];

  scalar omega[];
  vorticity(u, omega);

  view (fov = 2, camera = "front", 
        tx = -0.25, ty = -0.5, bg = {1,1,1},
	width = 1000, height = 1000);
  clear();
  draw_vof ("vof", "sf", filled = 1, lw = 3);
  isoline ("p", n = 120, min = -1, max = 1);
  sprintf (name, "piso-%g.png", t);
  save (name);

  view (fov = 2, camera = "front", 
        tx = -0.275, ty = -0.5, bg = {1,1,1},
	width = 1000, height = 1000);
  clear();
  draw_vof ("vof", "sf", filled = 1, lw = 3);
  squares ("omega", map = cool_warm);
  sprintf (name, "foil-%g.png", t);
  save (name);
  
}

event adapt (i++) {
  adapt_wavelet ({vof,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = LEVEL - 6);
}

