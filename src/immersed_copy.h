extern coord vc; // solid velocity
extern scalar vof;
extern face vector sf;
extern int maxlevel;

vector fc[]; // body force term (acceleration units)


double x2 = 4.75342;
double y2 = 2.99561;
void pointPrint (Point point, vector Fc, vector Pc, coord uc, int j, double deltah) {
  if (x >= x2*0.999 && x <= x2*1.001 && y >= y2*0.999 && y <= y2*1.001)
    fprintf (stderr, "%d x=%g y=%g u.x=%g u.y=%g fc.x=%g fc.y=%g Fc.x=%g Fc.y=%g Pc.x=%g Pc.y=%g uc.x=%g uc.y=%g dh=%g vof=%g\n", j, x, y, u.x[], u.y[], fc.x[], fc.y[], Fc.x[], Fc.y[], Pc.x[], Pc.y[], uc.x, uc.y, deltah, vof[]);
}

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
  // fprintf (stderr, "2 x=%g h=%g r=%g absr=%g\n", x, h, r, fabs(r));
  if (fabs(r) <= 1) {
  //  fprintf (stderr, "2.1 yes\n");
    phi = (1./8.) * (3 - (2 * fabs(r)) + sqrt (1 + (4 * fabs(r)) - (4 * sq(r))));

  }
  else if (fabs(r) > 1 && fabs(r) <= 2) {
    // fprintf (stderr, "2.2 yes\n");
    phi = (1./8.) * (5 - (2 * fabs(r)) - sqrt (-7 + (12 * fabs(r)) - ( 4 * sq(r))));
  }
  else {
    // fprintf (stderr, "2.3 yes\n");
    phi = 0;
  }
  return phi;
}

double delta_func (double x, double y, double xc, double yc, double Delta) {  // modify to include 3D
 // fprintf (stderr, "1 x=%g y=%g xc=%g yc=%g delta=%g\n", x, y, xc, yc, Delta); 
  double phi_x = phi_func (x - xc, Delta);
  double phi_y = phi_func (y - yc, Delta);
// fprintf (stderr, "3 phix=%g phiy=%g\n", phi_x, phi_y);
  return (phi_x * phi_y) / sq(Delta);
}

event solid_force (i++) {
  vector Pc[];  // coordinates of interface center (i.c)
  vector Fc[];  // force at i.c
  coord uc, pc;
  int j = 0;
  // trash({fc});

  foreach() {
    if (vof[] > 0 && vof[] < 1) {
      j = 0;
      pointPrint (point, Fc, Pc, uc, j, 0);

      coord m = {x, y};
      coord n = interface_normal (point, vof);
      double alpha = plane_alpha (vof[], n);
      double area = plane_area_center(n, alpha, &pc);
      foreach_dimension()
        pc.x = m.x + pc.x*Delta;


      j = 1;
      pointPrint (point, Fc, Pc, uc, j, 0);

      bilinear_interpolation (point, u, pc, &uc); // find u_c

      j = 2;
      pointPrint (point, Fc, Pc, uc, j, 0);

      foreach_dimension() {
        Fc.x[] = (vc.x - uc.x)/dt; // find f_c
	Pc.x[] = pc.x;
      }
      // fprintf (stderr, "1 x=%g y=%g xc=%g yc=%g delta=%g vof=%g\n", x, y, Pc.x[], Pc.y[], Delta, vof[]); 


      j = 3;
      pointPrint (point, Fc, Pc, uc, j, 0);    
    }

    else if (vof[] == 1) {
      foreach_dimension() {
        Fc.x[] += (vc.x - u.x[])/dt;
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

  foreach_dimension() {
    Pc.x.refine = Fc.x.refine = refine_linear;
  }

  foreach() {
    coord sum = {0};
    j = 40;
    pointPrint (point, Fc, Pc, uc, j, 0);
    if (vof[] == 1)
      foreach_dimension()
        sum.x = Fc.x[];
    else if (level == maxlevel) {
      double x1 = x;
      double y1 = y;   
      // fprintf (stderr, "2 i=%d x=%g y=%g xc=%g yc=%g delta=%g vof=%g\n", i, x, y, Pc.x[], Pc.y[], Delta, vof[]); 
      foreach_neighbor() {
      // fprintf (stderr, "2.5 x=%g y=%g xc=%g yc=%g delta=%g vof=%g\n", x, y, Pc.x[], Pc.y[], Delta, vof[]); 
        if (vof[] > 0 && vof[] < 1 && level == maxlevel) {

        j = 4;
        pointPrint (point, Fc, Pc, uc, j, 0);

        double delta_h = delta_func(x1, y1, Pc.x[], Pc.y[], Delta);

        // fprintf (stderr, "3 i=%d x=%g y=%g xc=%g yc=%g delta=%g vof=%g dh=%g\n", i, x, y, Pc.x[], Pc.y[], Delta, vof[], delta_h);
        j = 42;
        pointPrint (point, Fc, Pc, uc, j, delta_h);

        foreach_dimension()
            sum.x += Fc.x[] * delta_h * dv();
        }
      }
    }
    foreach_dimension()
      fc.x[] = sum.x;

    j = 5;
    pointPrint (point, Fc, Pc, uc, j, 0);
  }

  correction(-dt);
  foreach()
    foreach_dimension()
      g.x[] += (fc.x[]/dv());
  correction(dt);
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

void immersed_force (scalar c, face vector fc, coord * F) {
  coord Fi = {0, 0};
   foreach() {
    double area  = c[]*dv(); // dv() = sq(Delta) in 2D
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
