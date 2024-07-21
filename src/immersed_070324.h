#include "fractions.h"

extern coord vc; // solid velocity
extern scalar vof;
extern face vector sf;
extern int maxlevel;

vector fc[]; // body force term (acceleration units)

double x2 = 5.19287;
double y2 = 3.17139;


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

void printx (Point point, int j, vector Fc) {
  if (x > x2*0.999 && x < x2*1.001 && y > y2*0.999 && y < y2*1.001) {
      coord pc, uc;
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      bilinear_interpolation (point, u, pc, &uc); // find u_c

      double temp = (vc.x - uc.x);

    //  fprintf(stderr, "|| %d x=%g y=%g pc.x=%g pc.y =%g uc.x=%g uc.y=%g, u.x=%g, u.y=%g fc.x=%g fc.y=%g Fc.x=%g Fc.y=%g\n", j, x, y, pc.x, pc.y, uc.x, uc.y, u.x[], u.y[], fc.x[], fc.y[], Fc.x[], Fc.y[]);

  }
}

event solid_force (i++) {
  vector Pc[];  // coordinates of interface center (i.c)
  vector Fc[];  // force at i.c

  coord uc, pc;
  trash({fc});

  for (int i = 0; i < 3; i++) {
  foreach() {
    if (vof[] > 0 && vof[] < 1) {
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      bilinear_interpolation (point, u, pc, &uc); // find u_c
      
      foreach_dimension() {
        Fc.x[] = (vc.x - uc.x)/dt; // find f_c
	Pc.x[] = pc.x;
      }

      printx(point, 0, Fc);
    }
    else if (vof[] == 1) {
      foreach_dimension() {
        // Fc.x[] = (vc.x - u.x[])/dt;
	Fc.x[] = 0.;
        Pc.x[] = 0.;
      }
    }
    else  // cell is empty
      foreach_dimension() {
	Fc.x[] = 0;
	Pc.x[] = 0;
      }
  }
  
  boundary({Fc, Pc, vof});

  foreach() {
    coord sum = {0};
   //if (vof[] > 0 && vof[] < 1)
    //  foreach_dimension()
      //  sum.x = Fc.x[];
   printx (point, 1, Fc);
   if (level == maxlevel) {
      double x1 = x;
      double y1 = y;
      int trig =  (x > x2*0.999 && x < x2*1.001 && y > y2*0.999 && y < y2*1.001)? 1:0;
      double delta_h;
      foreach_neighbor() {
        if (vof[] > 0 && vof[] < 1 && level == maxlevel) {
          delta_h = delta_func(x1, y1, Pc.x[], Pc.y[], Delta);
          foreach_dimension()
            sum.x += Fc.x[] * delta_h * dv();
        
	if(trig)
         fprintf(stderr,"|| 2 x=%g y=%g fc.x=%g fc.y=%g dh=%g sum.x=%g sum.y=%g r.x=%g r.y=%g\n", x, y, fc.x[], fc.y[], delta_h, sum.x, sum.y, (x1 - Pc.x[])/Delta, (y1 - Pc.y[])/Delta);
	}
      }
    }
    foreach_dimension()
      fc.x[] = sum.x;
  printx (point, 3, Fc);
  }
  

  correction(-dt);
  foreach()
    foreach_dimension()
      g.x[] += (fc.x[]);
  correction(dt);

  }
  foreach()
    if (vof[] > 0 && vof[] < 1) {
      printx (point, 4, Fc);
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

      foreach_dimension()
        pc.x = m.x + pc.x*Delta;

      bilinear_interpolation (point, u, pc, &uc); // find u_c
  if (x > x2*0.999 && x < x2*1.001 && y > y2*0.999 && y < y2*1.001) {      
       fprintf (stderr, "|| 5 x=%g y=%g pc.x=%g pc.y=%g uc.x=%g uc.y=%g Fc.x=%g Fc.y=%g fc.x=%g fc.y=%g u.x=%g u.y=%g\n", x, y, pc.x, pc.y, uc.x, uc.y, Fc.x[], Fc.y[], fc.x[], fc.y[], u.x[], u.y[]);
    }
    }
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
