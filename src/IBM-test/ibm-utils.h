#ifndef BASILISK_HEADER_21
#define BASILISK_HEADER_21
#line 1 "./../ibm-utils.h"

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

double quadratic_interpolation (Point point, double xc, scalar u, int i, int j, coord n)
{
    if (fabs(n.x) >= fabs(n.y)) {
        double x1 = 1; 
    }
 
}


void printx (Point point, scalar vof, vector u, coord n, double d[], coord v[], coord gradient)
{
    int k = fabs(n.x) > fabs(n.y)? 1: 2;
    fprintf (stderr, "|| %d %g %g %g u.x=%g u.y=%g n.x=%g n.y=%g d[0]=%g d[1]=%g v0.x=%g v0.y=%g v1.x=%g v1.y=%g grad.x=%g grad.y=%g\n",
             k, x, y, vof[], u.x[], u.y[], n.x, n.y, d[0], d[1], v[0].x, v[0].y, v[1].x, v[1].y, gradient.x, gradient.y); 
}

void printy (Point point, scalar vof, vector u, coord n, double d[], coord v[], double xc, double yc, int i, int j, int l, coord pm)
{
    int k = fabs(n.x) > fabs(n.y)? 1: 2;
    fprintf (stderr, "|| %d %d %g %g %g u.x=%g u.y=%g n.x=%g n.y=%g d[0]=%g d[1]=%g v0.x=%g v0.y=%g v1.x=%g v1.y=%g x%d=%g y%d=%g\n|| i = %d j = %d pm.x=%g pm.y=%g atan2=%g\n",
             k, l, x, y, vof[], u.x[], u.y[], n.x, n.y,
             d[0], d[1], v[0].x, v[0].y, v[1].x, v[1].y, l,xc,l,yc,
             i, j, pm.x, pm.y, atan2(n.y, n.x)); 
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
