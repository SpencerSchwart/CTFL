#ifndef BASILISK_HEADER_19
#define BASILISK_HEADER_19
#line 1 "./../immersed.h"
extern coord vc; // solid velocity
extern scalar airfoil;
extern face vector sf;
vector aF[]; // body force acceleration term

#define vof_part_x(a,i) (a[i] > 0 && a[i] < 1? 1:0)
#define vof_part_y(a,i) (a[0,i] > 0 && a[0,i] < 1? 1:0)
#define vof_full_x(a) (a[] == 1 || a[1] == 1? 1:0)
#define vof_full_y(a) (a[0,0] == 1 || a[0,1] == 1? 1:0)


event end_timestep (i++) {
  foreach()
    foreach_dimension() {
      aF.x[] = airfoil[]*(vc.x-u.x[])/dt;
    }
  correction(-dt);
  foreach()
    foreach_dimension()
      g.x[] += aF.x[];
  correction(dt);
}

/*
event end_timestep (i++) {
  foreach_face()
    aF.x[] = face_value(airfoil,0)*(vc.x-face_value(u.x,0))/(dt);
  a = aF;
}
*/


void normal_vector (Point point, face vector fc, coord * nv) {
  coord n;
  foreach_dimension() {
    n.x = fc.x[] - fc.x[1];
  }
  double nn = sqrt(sq(n.x) + sq(n.y));
  if (nn > 0)
    foreach_dimension()
      n.x /= nn;
  else
    foreach_dimension()
      n.x /= 1./dimension;
  *nv = n;
}

/*
double int_area (Point point, coord * b, coord * n) {
  *n = facet_normal (point, airfoil, sf);
  double alpha = plane_alpha (airfoil[], *n);
  double area = plane_area_center (*n, alpha, b);
  normalize (n);
  return area;
}

void interface_area_field (scalar c, scalar a) {
  foreach() {
    if (c[] > 0 && c[] < 1) {
      coord n = interface_normal (point, c), p;
      double alpha = plane_alpha (c[], n);
      a[] = pow(Delta, dimension - 1)*plane_area_center (n, alpha, &p);
    }
  }
}
*/
void immersed_force (scalar c, face vector fc, coord * F) {
  coord Fi = {0, 0};
   foreach() {
    double area  = c[]*dv(); // dv() = sq(Delta) in 2D
    foreach_dimension()
      Fi.x += - fm.x[]*(aF.x[])*area;
   }
  *F = Fi; 
}

trace
int fractions_cleanup (scalar c, face vector s,
		       double smin = 0., bool opposite = false) {

  int changed = 1, schanged = 0, i;
  for (i = 0; i < 100 && changed; i++) {

  foreach_face()
      if (s.x[] && ((!c[] || !c[-1]) || s.x[] < smin))
	s.x[] = 0.;

    changed = 0;
    foreach(reduction(+:changed))
      if (c[] > 0. && c[] < 1.) {
	int n = 0;
	foreach_dimension() {
	  for (int i = 0; i <= 1; i++)
	    if (s.x[i] > 0.)
	      n++;

	  if (opposite && s.x[] == 0. && s.x[1] == 0.)
	    c[] = 0., changed++;
	}  
	if (n < dimension)
	  c[] = 0., changed++;
      }

    schanged += changed;
  }
  if (changed)
    fprintf (stderr, "WARNING: fractions_cleanup() did not converge after "
	     "%d iterations\n", i);
  return schanged;
}

#endif
