// #include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
// #include "../contact-embed3D.h"
#include "tag.h"
#include "view.h"

#define LEVEL 9
#define RHOR 80

int j;
double uemax = 1.;
double t_end = 2.5 [0,1];

const double H = 1 [1];  // gas nozzle 
const double R = 3*H;  // liquid inlet
const double RHO1 = 1;
const double SIGMA = 1.;

double MUl[4] = {2.e-2, 5.e-3, 5.e-3, 5.e-3};
double MUg[4] = {1.e-1, 2.5e-2, 2.5e-2, 2.5e-2};
double Ul[4] = {2, 2, 1, 0.5};
double Ug[4] = {40, 160, 160, 160};
double theta0 = 70;
double thetai = 45;

#define cylinder(a, b, r) (sq(a) + sq(b) - sq(r)) // or cone if x = r
#define cone(a, b, r, angle) (((sq(a) + sq(b))/sq(tan(angle))) - sq(r))
#define distance(a, b) (sqrt(sq(a) + sq(b)))

u.n[left] = dirichlet	(distance(y, z) <= 1.1*R? Ul[j]: Ug[j]);
u.t[left] = dirichlet	(0);
u.r[left] = dirichlet	(0);
p[left]   = neumann	(0);
pf[left]  = neumann	(0);
f[left]	  = dirichlet	(distance(y, z) <= 1.1*R? 1: 0);

u.n[right] = neumann 	(0);
uf.n[right] = neumann 	(0);
p[right] = dirichlet 	(0);
pf[right] = dirichlet 	(0);
f[right] = neumann	(0);

u.n[top] = dirichlet	(0);
u.t[top] = dirichlet 	(0);
u.r[top] = dirichlet	(0);
p[top] = neumann 	(0);
pf[top] = neumann 	(0);
f[top] = neumann	(0);

u.n[bottom] = dirichlet	(0);
u.t[bottom] = dirichlet (0);
u.r[bottom] = dirichlet	(0);
p[bottom] = neumann 	(0);
pf[bottom] = neumann 	(0);
f[bottom] = neumann	(0);

u.n[front] = dirichlet	(0);
u.t[front] = dirichlet 	(0);
u.r[front] = dirichlet	(0);
p[front] = neumann 	(0);
pf[front] = neumann 	(0);
f[front] = neumann	(0);

u.n[back] = dirichlet	(0);
u.t[back] = dirichlet	(0);
u.r[back] = dirichlet	(0);
p[back] = neumann 	(0);
pf[back] = neumann 	(0);
f[back] = neumann	(0);


// double u0 = 0 [1, -1];
// u.n[embed] = dirichlet (u0); // gives error somehow?
uf.n[embed] = dirichlet (0);
u.t[embed] = dirichlet	(0.);
uf.t[embed] = dirichlet	(0.);
p[embed] = neumann	(0);
pf[embed] = neumann	(0);
f[embed] = neumann	(0); // = 90 degree C.A

void solid_domain(scalar s, face vector sf) {
  double chamfer = 1*H/tan(thetai*M_PI/180); // was 2*H
  double theta_cone = (180 - (thetai*2)) / 2; // or 90 - thetai
  double offset = (tan(theta_cone*M_PI/180) * (5*H)) - chamfer;

  vertex scalar phi[];
  foreach_vertex() {
    if (x < -chamfer)
      phi[] = -min(max(cylinder(y, z, 2*R), -cylinder(y, z, 5*H)),cylinder(y, z, R));
    else if (x >= -chamfer && x < 0) {
      phi[] = -min(max(cylinder(y, z, 6*H), -(cone(y, z, x - offset, thetai*M_PI/180))),cylinder(y, z, R));
    }
    else if (x >= 0 && x <= H)
      phi[] = -cylinder(y, z, 2*R);
    else if (x > H && x <= 2*H)
      phi[] = -cylinder(y, z, R);
    else
      phi[] = 1;
  }
  boundary({phi});
  fractions (phi, s, sf);
  fractions_cleanup (s, sf);
  s.refine = s.prolongation = fraction_refine;
}


int main (int argc, char * argv[]) {
  L0 = 35*H;
  size (L0);
  init_grid (2 << (LEVEL - 3));
  origin (-5*R, -(L0/2), -(L0/2));
  
  rho1 = RHO1;
  rho2 = rho1 / RHOR;
  f.sigma = SIGMA;

  j = 3;
  mu1 = MUl[j];
  mu2 = MUg[j];

  // const scalar c[] = theta0*pi/180;
  // contact_angle = c;

  TOLERANCE = 1e-4 [*];

  run();
}


scalar f1[];
scalar cs1[];
scalar cs2[];

event init (t = 0) {
 //  mask(y > 7*H? top: y < -7*H ? bottom : none);
 //  mask(z > 7*H? front: z < -7*H? back: none);

  astats ss;
  int ic = 0;
  do {
    ic++;
    solid_domain(cs, fs);
    fraction (f, x < -H && distance (y, z) < 3.1*H);
    ss = adapt_wavelet ({cs, f}, (double[]) {1.e-30, 1.e-30},
		        maxlevel = LEVEL , minlevel = LEVEL - 4); 
  } while ((ss.nf || ss.nc) && ic < 100);

  refine (((cs[] > 0 && cs[] < 1) || (f[] > 0 && f[] < 1)) && level < LEVEL);
  solid_domain(cs, fs);
  fraction (f, x < -H && distance(y, z) < 3.1*H);

  foreach()
    foreach_dimension()
      u.x[] = 0.;
}

scalar error[];
int e = 0;
event logfile (i++; t <= t_end) {
  foreach_face(reduction(max:e))
      if (!fm.x[] && uf.x[] != 0.) {
	error[] = 1.;
        fprintf (stderr, "%g %g %g %g %g %g %g %g\n", x, y, z, f[], cs[], fm.x[], fm.y[], cm[]);
	e = 1;
      }
      else
        error[] = 0.;
  if (e == 1) {
    char name[80];
    sprintf (name, "error-%d", i);
    FILE * fp1 = fopen(name, "w");
     
    view (fov = 15, camera = "iso", tx = -0.15);
    clear();
    draw_vof ("error");
    save (fp = fp1);
  }
  fprintf (stderr, "%d %g %d\n", i, t, j);
}


event movie (t += 0.0125; t <= t_end) {
  scalar umag[];
  foreach() {
    f1[] = x < 15.9*H? f[]: 0;
    cs1[] = distance(y, z) <= 2.1*R? cs[]: 1;
    cs2[] = x <= 0 && distance(y, z) <= 2.1*R? cs[]: 1;
    umag[] = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
  }

  char name[80];
  sprintf (name, "img/vof-%g.png", t);
  view (fov = 15, camera = "iso", tx = -0.15);
  clear();
  draw_vof ("f", fc = {1, 0.522, 0.439});
  draw_vof ("cs1");
  save (file = name);

  sprintf (name, "img/2vof-%g.png", t);
  view (fov = 15, camera = "iso", tx = -0.15);
  clear();
  draw_vof ("f", fc = {1, 0.522, 0.439});
  draw_vof ("cs2");
  save (file = name);

  sprintf (name, "img/velo-%g.png", t);
  view (fov = 20, camera = "front", tx=-0.075, bg = {1,1,1} );
  clear();
  squares ("umag");
  save (file = name);

  sprintf (name, "img/f-%g.png", t);
  clear();
  squares ("f");
  save (file = name);

  sprintf (name, "img/cells-%g.png", t);
  clear();
  squares ("umag");
  cells();
  save (file = name);

  sprintf(name, "dumps/dump-%g", t);
  dump (file = name);
}

/*
event vprof (t += 0.5) {
  char name[80];
  sprintf (name, "vprof1-%g", t);
  FILE * fpv1 = fopen (name, "w");
  double delta = L0/pow(2,LEVEL);
  for (int j = R + 2*H; j < R + 3*H; j += delta)
    foreach_point(-H, j, 0)
      fprintf (fpv1, "%g %g %g %g %g %g\n", x, y, z, u.x[], u.y[], u.z[]);
  fflush(fpv1);
  fclose(fpv1);

  sprintf (name, "vprof2-%g", t);
  FILE * fpv2 = fopen (name, "w");
  for (int i = 0; i < H; i += delta)
    foreach_point(i, R + (H/2), 0)
      fprintf (fpv2, "%g %g %g %g %g %g\n", x, y, z, u.x[], u.y[], u.z[]);
  fflush(fpv2);
  fclose(fpv2);
}
*/


event droplets (i+=5) {
  remove_droplets (f);
  // fractions_cleanup (cs, fs);
}


event adapt (i++) {	
  adapt_wavelet ({cs,f,u}, (double[]){1.e-4, 1.e-4, uemax,uemax,uemax}, 
		  maxlevel = LEVEL, minlevel = LEVEL - 5);
  unrefine (x > 16*H && level > 3);
}

