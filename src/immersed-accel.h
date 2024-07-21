#include "fractions.h"

extern coord vc; // solid velocity
extern scalar vof;
extern face vector sf;
extern int maxlevel;

face vector afc[];
vector fc[]; // body force term (acceleration units)
vector Pc[];  // coordinates of interface center (i.c)
vector Fc[];  // acutal force at i.c
vector Fd[];  // desired force at i.c

double x2 = 4.75342;
double y2 = 2.99561;


void bilinear_interpolation (Point point, vector uv, coord pc, coord * uc) {
  coord uci;
  double xc = pc.x;
  double yc = pc.y;

  double xx = x - xc;
  double yy = y - yc;

  double x2 = x + sign(xx)*Delta;
  double y2 = y + sign(yy)*Delta;

  foreach_dimension() {
  uci.x = (x2 - xc)*(y2 - yc) * uv.x[] + \
		    (xc - x) * (y2 - yc) * uv.x[sign(xx)] + \
		    (x2 - xc) * (yc - y) * uv.x[0,sign(yy)] + \
		    (xc - x) * (yc - y) * uv.x[sign(xx),sign(yy)];
  uci.x /=  ((x2 - x)*(y2 - y));
  }

  *uc = uci;
}

double phi_func (double x, double h) {
  double r = x / h;
  double phi;
  if (fabs(r) <= 1)
    phi = (1./8.) * (3 - (2 * fabs(r)) + sqrt (1 + (4 * fabs(r)) - (4 * sq(r))));
  else if (fabs(r) > 1 && fabs(r) <= 2)
    phi = (1./8.) * (5 - (2 * fabs(r)) - sqrt (-7 + (12 * fabs(r)) - ( 4 * sq(r))));
  else
    phi = 0;
  return phi;
}

double delta_func (double x, double y, double xc, double yc, double Delta) {  // modify to include 3D
  double phi_x = phi_func (x - xc, Delta);
  double phi_y = phi_func (y - yc, Delta);
  return (phi_x * phi_y) / sq(Delta);
}

void printx (Point point, int j) {

  if (x > x2*0.999 && x < x2*1.001 && y > y2*0.999 && y < y2*1.001 && vof[] < 1 && vof[] > 0) {
      coord pc, uc;
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      coord sum = {0};
      foreach_neighbor() {
        double delta_u = delta_func (x, y, pc.x, pc.y, Delta);
        foreach_dimension()
	  sum.x += u.x[]*delta_u*dv();
      }
      foreach_dimension()
        uc.x = sum.x;
      
    fprintf(stderr, "|| %d x=%g y=%g pc.x=%g pc.y =%g uc.x=%g uc.y=%g u.x=%g u.y=%g fc.x=%g fc.y=%g Fd.x=%g Fd.y=%g Fc.x=%g Fc.y=%g\n", j, x, y, pc.x, pc.y, uc.x, uc.y, u.x[], u.y[], fc.x[], fc.y[], Fd.x[], Fd.y[], Fc.x[], Fc.y[]);
  }

}


event acceleration (i++) {
  coord uc, pc;
  trash({fc, Pc, Fc, Fd});

  foreach() {
    if (vof[] < 1 && vof[] > 0) {
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      // find u_c
      coord sum = {0};
      foreach_neighbor() {
        double delta_u = delta_func (x, y, pc.x, pc.y, Delta);
        foreach_dimension()
	  sum.x += u.x[]*delta_u*dv();
      }
      
      foreach_dimension() {
	uc.x = sum.x;
        Fc.x[] = Fd.x[] = (vc.x - uc.x)/dt; // find f_c
	Pc.x[] = pc.x;
      }
    }
    else  // cell is empty or full
      foreach_dimension() {
	Fc.x[] = Fd.x[] = 0;
	Pc.x[] = 0;
      }
    printx (point, 1);
  }
  
  boundary({Fc, Pc, vof});

  // spread F_c to f_c
  foreach() {
    coord sum = {0};
    if (level == maxlevel) {
      double x1 = x;
      double y1 = y;
      double delta_h;
      foreach_neighbor() {
        if (vof[] > 0 && vof[] < 1 && level == maxlevel) {
          delta_h = delta_func(x1, y1, Pc.x[], Pc.y[], Delta);
          foreach_dimension()
            sum.x += Fd.x[] * delta_h * dv();
	}
      }
    }
    foreach_dimension()
      fc.x[] = sum.x;
    printx (point, 2);
  }

  // correct F_c
  coord num = {0};
  coord den = {0};
  coord correct;
  foreach() {
    coord pc, sum = {0};
    if (vof[] > 0 && vof[] < 1 && Pc.x[]) {
      double delta_sum = 0;
      pc.x = Pc.x[];
      pc.y = Pc.y[];
      foreach_neighbor() { // ensure its max level to avoid ghost cells?
        double delta_f = delta_func (x, y, pc.x, pc.y, Delta);
        foreach_dimension()
	  sum.x += fc.x[] * delta_f * dv();
	// delta_sum += delta_f * dv();
      }
      foreach_dimension() {
        Fc.x[] = sum.x;
        num.x += (Fc.x[] * Fd.x[]);
        den.x += (sq(Fc.x[]));
      // fprintf(stderr, "|| point: %g %g %g %g %g %g %g\n", x, y, Fc.x[], Fd.x[], num.x, den.x, delta_sum);
    }
    printx (point, 3);
    }
  }
  foreach_dimension()
    correct.x = den.x != 0? num.x/den.x: 1; 
  fprintf (stderr, "|| 3.5 correct.x=%g correct.y=%g num.x=%g den.x=%g\n", correct.x, correct.y, num.x, den.x);

  // spread F_c to f_c with correction
  foreach() {
    coord sum = {0};
    if (level == maxlevel) {
      double x1 = x;
      double y1 = y;
      double delta_h;
      foreach_neighbor() {
        if (vof[] > 0 && vof[] < 1 && level == maxlevel) {
          delta_h = delta_func(x1, y1, Pc.x[], Pc.y[], Delta);
          foreach_dimension()
            sum.x += Fd.x[] * delta_h * dv();
	}
      }
    }
    foreach_dimension() {
      // fc.x[] = vof[] == 1? (vc.x - u.x[])/dt: correct.x * sum.x;
      fc.x[] = correct.x * sum.x;
    }
    printx (point, 4);
  }

  /*
  correction(-dt);
  foreach()
    foreach_dimension()
      g.x[] += (fc.x[]);
  correction(dt);
*/

  foreach_face()
    afc.x[] = face_value(fc.x,0);
  a = afc;

  foreach()
    printx (point, 5);
}


void normal_vector_face (Point point, face vector fc, coord * nv) {
  coord n;
  foreach_dimension()
    n.x = fc.x[] - fc.x[1];
  double nn = sqrt(sq(n.x) + sq(n.y));
  if (nn > 0)
    foreach_dimension()
      n.x /= nn;
  else
    foreach_dimension()
      n.x /= 1./dimension;
  *nv = n;
}

void normal_vector (Point point, scalar s, coord * nv) {
  coord n;
  foreach_dimension()
      n.x = (s[1] - s[-1])/(2.*Delta);
  double mag = 0;
  foreach_dimension()
    mag += sq(n.x);
  mag = sqrt(mag);
  if (mag > 0)
    foreach_dimension()
      n.x /= mag;
  *nv = n;
}

/*
void immersed_force (scalar c, face vector fc, coord * F) {
  coord Fi = {0, 0};
   foreach() {
    double area  = c[]*dv(); // dv() = sq(Delta) in 2D
    foreach_dimension()
      Fi.x += - fm.x[]*(fc.x[])*area;
   }
  *F = Fi; 
}
*/


void immersed_forcev2 (scalar c, coord * F) {
   coord Fi = {0, 0};
   foreach() {
    double area  = dv(); // dv() = sq(Delta) in 2D
    foreach_dimension()
      Fi.x += - fm.x[]*(fc.x[])*area;
   }
  *F = Fi; 
}


/* trace
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
*/
