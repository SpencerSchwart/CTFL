#ifndef BASILISK_HEADER_20
#define BASILISK_HEADER_20
#line 1 "./../ibm-utils.h"
extern scalar vof;
extern face vector sf;
#define distance(a,b) sqrt(sq(a) + sq(b))

void bilinear_interpolation (Point point, vector uv, coord pc, coord * uc)
{
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

double scalar_bilinear_interpolation (Point point, scalar p, coord pc) 
{
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

double phi_func (double x, double h) 
{
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

double delta_func (double x, double y, double xc, double yc, double Delta)
{
    double phi_x = phi_func (x - xc, Delta);
    double phi_y = phi_func (y - yc, Delta);
    return (phi_x * phi_y) / sq(Delta);
}

bool empty_neighbor (Point point, coord * pc, scalar vof)
{
    coord pc_temp;
    double temp_vof = vof[];
    double xc = x;
    double yc = y;
    double max_d = 1e6;
    int neighbor = 0;

    foreach_neighbor(1)
        if (vof[] == 0 && temp_vof == 1 && (distance(x - xc, y - yc) < max_d)) {
            pc_temp.x = (xc + x) / 2.;
            pc_temp.y = (yc + y) / 2.;
            max_d = distance(x - xc, y - yc);
            neighbor = 1;
            *pc = pc_temp;
        }
   return neighbor;
}

double marker_point (Point point, scalar vof, coord * markerPoint)
{
    coord cellCenter = {x, y};
    coord n = interface_normal (point, vof);
    double alpha = plane_alpha (vof[], n);
    double area = plane_area_center (n, alpha, markerPoint);

    foreach_dimension()
        markerPoint->x = cellCenter.x + markerPoint->x*Delta;
    return area;
}

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar cs,
					   coord n, coord p, coord bc,
					   double * coef)
{
    double d[2] = {0,0}, v[2] = {nodata,nodata};
    bool defined = true;
    for (int l = 0; l <= 1; l++) {
        int i = (l + 1)*sign(n.x);
        d[l] = (i - p.x)/n.x;
        double y1 = p.y + d[l]*n.y;
        int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
        y1 -= j;
        if (cs[i,j-1] < 0.5 && cs[i,j] < 0.5 && cs[i,j+1] < 0.5)
            v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
    }
    if (v[0] == nodata) {
        d[0] = max(1e-3, fabs(p.x/n.x));
        *coef = - 1./(d[0]*Delta);
        return bc.x/(d[0]*Delta);
    }
    *coef = 0.;
    double gradient = 0;
    if (v[1] != nodata) // third-order gradient
        gradient = (d[1]*(bc.x - v[0])/d[0] - d[0]*(bc.x - v[1])/d[1])/((d[1] - d[0])*Delta);
    else
        gradient = (bc.x - v[0])/(d[0]*Delta); // second-order gradient
    return gradient;
}


double dirichlet_gradient (Point point, scalar s, scalar cs,
			   coord n, coord p, coord bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradient_x (point, s, cs, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradient_y (point, s, cs, n, p, bc, coef);
  return dirichlet_gradient_z (point, s, cs, n, p, bc, coef);
#endif // dimension == 3
  return nodata;
}

double embed_geometry (Point point, coord * b, coord * n)
{
  *n = facet_normal (point, vof, sf);
  double alpha = plane_alpha (vof[], *n);
  double area = plane_area_center (*n, alpha, b);
  normalize (n);
  return area;
}

double embed_interpolate (Point point, scalar s, coord p)
{
  assert (dimension == 2);
  int i = sign(p.x), j = sign(p.y);
  if (vof[i] < 0.5 && vof[0,j] < 0.5 && vof[i,j] < 0.5)
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (vof[i])
	val += fabs(p.x)*(s[i] - s[]);
      else if (vof[-i])
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}

static inline
coord embed_gradient (Point point, vector u, coord p, coord n, coord bc)
{
  coord dudn;
  foreach_dimension() {
    bool dirichlet = true;
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, vof, n, p, bc, &val);
    }
    else // Neumann
      dudn.x = bc.x;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}


/*
void ibm_dirichlet_gradient (Point point, vector u, scalar vof, coord n, coord pm, coord vb, coord * val)
{
    double d[2] = {nodata, nodata};
    coord v[2];
    double x1, y1;
    for (int l = 0; l <= 1; l++) {
        if (fabs(n.x) >= fabs(n.y)) {
            int i = sign(n.x);
            x1 = (Delta)*((l + 1) + pm.x);
            y1 = sign(n.y)*fabs(atan(n.y/n.x)) * x1;
            d[l] = distance(x1, y1);
            int j = fabs(y1 + (Delta*pm.y)) < Delta/2? 0: y1 + (Delta*pm.y) > Delta/2? 1: -1;
            if (vof[i*(l + 1), j] < 0.5 && vof[i*(l + 1), j+1] < 0.5 && vof[i*(l + 1), j-1] < 0.5)
                foreach_dimension()
                    v[l].x = quadratic(y1, (u.x[i*(l+1),j]), (u.x[i*(l+1),j+1]), (u.x[i*(l+1), j-1]));
            else
                foreach_dimension()
                    v[l].x = 0.;
            printy (point, vof, u, n, d, v, x1, y1, i, j, l, pm);
        }
        else if (fabs(n.x) < fabs(n.y)) {
            int j = sign(n.y);
            y1 = (Delta)*((l + 1) - pm.y);
            x1 = y1 / (sign(n.x)*fabs(atan(n.y/n.x)));
            d[l] = distance(x1, y1);
            int i = fabs(x1 + (Delta*pm.x)) < Delta/2? 0: x1 + (Delta*pm.x) > Delta/2? 1: -1;
            if (vof[i, j*(l + 1)] < 0.5 && vof[i+1,j*(l + 1)] < 0.5 && vof[i-1,j*(l + 1)] < 0.5)
                foreach_dimension()
                    v[l].x = quadratic (x1, (u.x[i,j*(l+1)]), (u.x[i+1,j*(l+1)]), (u.x[i-1,j*(l+1)]));
            else
                foreach_dimension()
                    v[l].x = 0.;
            printy (point, vof, u, n, d, v, x1, y1, i, j, l, pm);
        }
    }
    coord gradient;
    foreach_dimension()
        gradient.x = (d[1]/d[0]*(vb.x - v[0].x) - (d[0]/d[1])*(vb.x - v[1].x))/(d[1] - d[0]);
    *val = gradient;
    printx (point, vof, u, n, d, v, gradient);
}


coord ibm_gradient (Point point, vector u, coord pm, coord n, coord vb, scalar vof)
{
    coord dudn;
    ibm_dirichlet_gradient (point, u, vof, n, pm, vb, &dudn);
    if (dudn.x == nodata)
        dudn.x = -0.;
    return dudn;
}
*/
















/*
foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar vof, face vector sf,
                                           coord n, coord p, double bc, double * coef)
{
    foreach_dimension()
        n.x = - n.x;
    double d[2], v[2] = {nodata,nodata};
    bool defined = true;
    foreach_dimension()
        if (defined && !sf.x[(n.x > 0.)])
            defined = false;
    if (defined)
        for (int l = 0; l <= 1; l++) {
            int i = (l + 1)*sign(n.x);
            d[l] = (i - p.x)/n.x;
            double y1 = p.y + d[l]*n.y;
            int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
            y1 -= j;
#if dimension == 2
            if (sf.x[i + (i < 0),j] && sf.y[i,j] && sf.y[i,j+1] &&
                vof[i,j-1] && vof[i,j] && vof[i,j+1])
                v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
#endif         
        }
        if (v[0] == nodata) {
            d[0] = max(1e-3, fabs(p.x/n.x));
            *coef = - 1./(d[0]*Delta);
            return bc/(d[0]*Delta);
        }
    *coef = 0.;
    if (v[1] != nodata) // third-order gradient
        return (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
    return (bc - v[0])/(d[0]*Delta); // second-order gradient
}

double dirichlet_gradient (Point point, scalar s, scalar vof, face vector sf,
			               coord n, coord p, double bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradient_x (point, s, vof, sf, n, p, bc, coef);
#endif
  return nodata;
}
*/

#endif
