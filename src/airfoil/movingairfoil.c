#include "navier-stokes/centered.h"
// #include "navier-stokes/double-projection.h"
#include "../immersed.h" // IBM
// #include "curvature.h"
#include "view.h"

#define MAX_THICKNESS 0.12
#define CHORD_LENGTH 1
#define L0 20.
#define Re (2000.)
#define LEVEL 13

double U0 =  1.0; // inlet velocity
double rr = 1.1019*sq(MAX_THICKNESS); // Radius of leading edge
// double theta_p = 5*pi/180; // aoa = 5 degrees
double theta_p = 20*pi/180; // aoa = 20 degrees

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


void airfoil_shape (scalar c, face vector f, double theta, vertex scalar phii = {0})
{
  double yt, yc = 0;
  vertex scalar phi = automatic (phii);
  
  double chord = CHORD_LENGTH;
  
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
  fractions_cleanup (c, f);
  c.refine = c.prolongation = fraction_refine;
}


double SDF (double xx, double yy) {
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
    airfoil_shape (airfoil, sf, theta_p);
    ss = adapt_wavelet ({airfoil}, (double[]) {1.e-30},
			maxlevel = LEVEL, minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
  airfoil_shape(airfoil, sf, theta_p);
}

scalar ref[];
face vector rf[];
vertex scalar phia[];

event moving_cylinder (i++) {
  airfoil_shape (airfoil, sf, theta_p, phia);
  airfoil_shape (ref, rf, theta_p);

  if (j == 0) { // SPM
  foreach() {
    coord nv;
    normal_vector (point, sf, &nv);
      double lambda = fabs(nv.x) + fabs(nv.y);
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

}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(CHORD_LENGTH)/(Re);
  boundary ((scalar *) {muv});
}


event logfile (i++; t <= 15) {

  coord Fp;
  immersed_force (airfoil, sf, &Fp);
  double CD = (Fp.x)/(0.5*sq(U0)*(CHORD_LENGTH));
  double CL = (Fp.y)/(0.5*sq(U0)*(CHORD_LENGTH));
/*
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
*/

  fprintf (stderr, "%d %g %d %d %d %d %d %g %g\n",
	   i, t, j, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL);
  
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

event movie (t += 0.01; t <= 10) {

}

event screenshot (t += 0.5; t <= 15) {
  char name[80];
  sprintf (name, "%d-piso-%g.png", j, t);
  FILE * fp1 = fopen (name, "w");
  view (fov = 2, camera = "front", 
        tx = -0.25, ty = -0.5, bg = {1,1,1},
	width = 3200, height = 3200);
  clear();
  draw_vof ("ref", "rf", filled = 1, lw = 3);
  isoline ("p", n = 120, min = -1, max = 1);
  save (fp = fp1);
  fclose (fp1);
}

event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-2,3e-3,3e-3},
		 maxlevel = LEVEL, minlevel = LEVEL - 6);
}
