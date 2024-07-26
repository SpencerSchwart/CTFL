
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

bool empty_neighbor (Point point, coord * pc) 
{
    coord pc_temp;
    double temp_vof = vof[];
    double xc = x;
    double yc = y;
    double max_d = 1e6;
    int neighbor = 0;

    foreach_neighbor(1) {
        if (vof[] == 0 && temp_vof == 1 && distance(x - xc, y - yc) < max_d) {
            pc_temp.x = (xc + x) / 2;
            pc_temp.y = (yc + y) / 2;
            max_d = distance (x - xc, y - yc);
            neighbor = 1;
            *pc = pc_temp;
        }
    }
   return neighbor;
}

extern double x2, y2;

void printx (Point point, int j, int i)
{
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
      
        fprintf(stderr, "||  %d %d x=%g y=%g vof= %g pc.x=%g pc.y =%g uc.x=%g uc.y=%g u.x=%g u.y=%g fc.x=%g fc.y=%g Fd.x=%g Fd.y=%g Fc.x=%g Fc.y=%g\n", 
        i, j, x, y, vof[], Pc.x[], Pc.y[], uc.x, uc.y, u.x[], u.y[],
        fc.x[], fc.y[], Fd.x[], Fd.y[], Fc.x[], Fc.y[]);
  }
}
