#include "fractions.h"

extern coord vc; // solid velocity
extern scalar vof;
extern int maxlevel;

#define distance(a, b) sqrt(sq(a) + sq(b))

vector fc[]; // body force term (acceleration units)
vector Pc[];  // coordinates of interface center (i.c)
vector Fc[];  // acutal force at i.c
vector Fd[];  // desired force at i.c
vector us[];  // velocity w/o solid force

double x2 = 4.7750;
double y2 = 10.1150;

void bilinear_interpolation (Point point, vector uv, coord pc, coord * uc) {
  coord uci;
  double xc = pc.x;
  double yc = pc.y;

  double xx = x - xc;
  double yy = y - yc;

  double x2 = x - sign(xx)*Delta;
  double y2 = y - sign(yy)*Delta;

  foreach_dimension() {
    uci.x = (x2 - xc)*(y2 - yc) * uv.x[] +
            (xc - x) * (y2 - yc) * uv.x[-sign(xx)] +
	    (x2 - xc) * (yc - y) * uv.x[0,-sign(yy)] +
	    (xc - x) * (yc - y) * uv.x[-sign(xx),-sign(yy)];
    uci.x /=  ((x2 - x)*(y2 - y));
  }

  *uc = uci;
}

double scalar_bilinear_interpolation (Point point, scalar p, coord pc) {
  double pci;
  double xc = pc.x;
  double yc = pc.y;

  double xx = x - xc;
  double yy = y - yc;

  double x2 = x - sign(xx)*Delta;
  double y2 = y - sign(yy)*Delta;

  pci = (x2 - xc)*(y2 - yc) * p[] +
            (xc - x) * (y2 - yc) * p[-sign(xx)] +
	    (x2 - xc) * (yc - y) * p[0,-sign(yy)] +
	    (xc - x) * (yc - y) * p[-sign(xx),-sign(yy)];
  pci /=  ((x2 - x)*(y2 - y));

  return pci;
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

bool empty_neighbor (Point point, coord * pc) {
  coord pc_temp;
  double temp_vof = vof[];
  double xc = x;
  double yc = y;
  double max_d = 1e6;
  int neighbor = 0;
  foreach_neighbor(1)
    if (vof[] == 0 && temp_vof == 1 && distance(x - xc, y - yc) < max_d) {
      pc_temp.x = (xc + x) / 2;
      pc_temp.y = (yc + y) / 2;
      max_d = distance (x - xc, y - yc);
      neighbor = 1;
      *pc = pc_temp;
    }
   return neighbor;
}

void printx (Point point, int j, int i) {

  if (x > x2*0.9999 && x < x2*1.0001 && y > y2*0.9999 && y < y2*1.0001) {
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
      
    fprintf(stderr, "||  %d %d x=%g y=%g vof= %g pc.x=%g pc.y =%g uc.x=%g uc.y=%g u.x=%g u.y=%g fc.x=%g fc.y=%g Fd.x=%g Fd.y=%g Fc.x=%g Fc.y=%g\n", i, j,
		    x, y, vof[], Pc.x[], Pc.y[], uc.x, uc.y, u.x[], u.y[], fc.x[], fc.y[], Fd.x[], Fd.y[], Fc.x[], Fc.y[]);
  }
}


event end_timestep (i++) {
  coord uc, pc;
  trash({fc, Pc, Fc, Fd});
  vector nv[];

  foreach() {
    foreach_dimension() {
      us.x[] = u.x[];
    }
    if (vof[] < 1 && vof[] > 0) {
      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);

      foreach_dimension() {
        pc.x = m.x + pc.x*Delta;
	nv.x[] = n.x;
      }

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
    else if (empty_neighbor(point, &pc)) {
      coord n = interface_normal (point, vof);
      coord sum = {0};
      foreach_neighbor() {
        double delta_u = delta_func (x, y, pc.x, pc.y, Delta);
        foreach_dimension()
          sum.x += u.x[]*delta_u*dv();
      }
      foreach_dimension() {
        uc.x = sum.x;
	Fc.x[] = Fd.x[] = (vc.x - uc.x)/dt;
        Pc.x[] = pc.x;
	nv.x[] = n.x;
      }
    }
    else
      foreach_dimension() {
	Fc.x[] = Fd.x[] = 0; 
	Pc.x[] = nv.x[] = 0;
      }

    printx (point, 1, i);
  }
  
  boundary({Fc, Pc, vof});

  // spread F_c to f_c
  foreach() {
    coord sum = {0};
    if (level == maxlevel) {
      double x1 = x;
      double y1 = y;
      double delta_h;
      foreach_neighbor()
        if (Pc.x[] && level == maxlevel) { 
          delta_h = delta_func(x1, y1, Pc.x[], Pc.y[], Delta);
          foreach_dimension()
            sum.x += Fd.x[] * delta_h * dv();
	    }
    }
     foreach_dimension()
       fc.x[] = sum.x;

    printx (point, 2, i);
  }

  // correct F_c
  coord num = {0};
  coord den = {0};
  coord correct;

  foreach() {
    coord pc, sum = {0};
    if (Fc.x[]) {
      pc.x = Pc.x[];
      pc.y = Pc.y[];
      foreach_neighbor() {
        double delta_f = level == maxlevel? delta_func (x, y, pc.x, pc.y, Delta): 0;
        foreach_dimension()
	  sum.x += fc.x[] * delta_f * dv();
      }
      foreach_dimension() {
        Fc.x[] = sum.x;
        num.x += (Fc.x[] * Fd.x[]);
        den.x += (sq(Fc.x[]));
      }
    }
    printx (point, 3, i);
  }
  foreach_dimension()
    correct.x = den.x != 0? num.x/den.x: 1;
  // fprintf (stderr, "|| correct.x=%g correct.y=%g\n", correct.x, correct.y);

  // spread F_c to f_c with correction
  foreach() {
    double delta_sum = 0.;
    coord sum = {0};
    if (level == maxlevel) {
      double x1 = x;
      double y1 = y;
      double delta_h;
      foreach_neighbor() {
        if (Pc.x[] && level == maxlevel) {
          delta_h = delta_func(x1, y1, Pc.x[], Pc.y[], Delta);
          foreach_dimension()
            sum.x += correct.x * Fd.x[] * delta_h * dv();
	      delta_sum += delta_h;
	    }
      }
    }
    foreach_dimension() {
      // fc.x[] = vof[] >= 0.5? (vc.x - u.x[])/dt: correct.x * sum.x;
      //  fc.x[] = vof[] == 1 && fc.x[] == 0? (vc.x - u.x[])/dt: correct.x * sum.x;
      Fc.x[] *= correct.x;
      fc.x[] = sum.x;
    }
    printx (point, 4, i);
  }


  foreach()
    printx (point, 5, i);

  correction(-dt);
  foreach()
    foreach_dimension()
      g.x[] += (fc.x[]);
  correction(dt);

  foreach() {
    if (vof[] >= 0.5)
      foreach_dimension() {
        // fc.x[] = (vc.x - u.x[])/dt;
	u.x[] = vc.x;
      }
  }

  foreach()
    printx (point, 6, i);

}


void immersed_forcev2 (scalar c, coord * desiredForce, coord * markerForce, coord * interpolateForce, coord * gridForce) {
  coord desiredForcei = {0};
  coord markerForcei = {0};
  coord interpolateForcei = {0};
  coord gridForcei = {0};

  foreach() {
    if (fc.x[]) {
      coord n = interface_normal (point, c), p, pc;
      coord sum = {0};

      double alpha = plane_alpha (c[], n);
      double area = pow(Delta, dimension - 1) * plane_area_center (n, alpha, &p);
      pc.x = Pc.x[];
      pc.y = Pc.y[];
      foreach_neighbor() {
        double delta_f = level == maxlevel? delta_func (x, y, pc.x, pc.y, Delta): 0;
        foreach_dimension() {
	      sum.x += fc.x[] * delta_f * dv();
        }
      }
      // bilinear_interpolation (point, fc, pc, &sum);
      foreach_dimension() {
        desiredForcei.x += - (Fd.x[])*area;
        markerForcei.x += - (Fc.x[])*area;
        interpolateForcei.x += - (sum.x)*area;
        gridForcei.x += - (fc.x[])*dv();
      }
    }
  }
  *desiredForce = desiredForcei;
  *markerForce = markerForcei;
  *interpolateForce = interpolateForcei;
  *gridForce = gridForcei;
}
/*
void immersed_forcev2 (scalar c, coord * F) {
   coord Fi = {0, 0};
   foreach() {
    double area  = dv(); // dv() = sq(Delta) in 2D
    foreach_dimension()
      Fi.x += - fm.x[]*(fc.x[])*area;
   }
  *F = Fi;
}
*/

