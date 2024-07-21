#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "../immersed.h"
#include "distance.h"
#include "lambda2.h"
#include "view.h"

#define LEVEL 9

double t_end = 1 [0,1];
double U0 = 1; // inlet velocity
double Re;
double uemax = 0.5;

coord vc = {0,0};

face vector muv[];
scalar airfoil[];
face vector sf[];

// Boudary Conditions

// ceiling
u.n[front] = dirichlet(U0);
u.t[front] = dirichlet(0);
u.r[front] = dirichlet(0);
p[front] = neumann(0);
pf[front] = neumann(0);

// outlet
u.n[back] = neumann(0);
p[back] = dirichlet(0);
pf[back] = dirichlet(0);


void fraction_from_stl (scalar cs, face vector fs, FILE * fp) {
  coord * p = input_stl (fp); // import CAD
  coord min, max;
  bounding_box (p, &min, &max);
  double maxl = -HUGE;
  foreach_dimension()
    if (max.x - min.x > maxl)
      maxl = max.x - min.x;

  // STL points to domain points
  scalar d[];
  distance (d, p);
  while (adapt_wavelet ({d}, (double[]){5e-4}, LEVEL, 5).nf);

  vertex scalar phi[];
  foreach_vertex(){;
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
	     d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.;
  }
  fractions (phi, cs, fs);
  fractions_cleanup (cs, fs);
}

void rotate_vof (scalar s, face vector fs) {
  vertex scalar phi[];
  foreach_vertex() {
    
  }
}

int main() {
  double L0 = 1 [1];
  size(L0);
  init_grid(2 << (LEVEL-5));
  origin (-L0/2, -L0/2, -L0/2);
  mu = muv;
 
  run();
}


event init (t = 0) {
    if (!restore (file = "restart")) {
    FILE * fp = fopen ("/home/spencer/basilisk/CTFL/files/prop_test1.stl", "r");
    if (fp == NULL)
      fprintf(stderr, "STL file is NULL\n");
    else {
      fprintf(stderr, "STL is NOT null\n");
      fraction_from_stl (airfoil, sf, fp);
      fclose (fp);
    }
  }
}

event logfile (i++; t <= t_end) {
  fprintf (stderr, "%d %g", i, t);
}

event adapt (i++) {
  adapt_wavelet ({airfoil,u}, (double[]){1.e-2, (uemax), (uemax), (uemax)},	
		  maxlevel = LEVEL, minlevel = 2);
}

