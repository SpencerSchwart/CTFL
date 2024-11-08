#define distance(a,b) sqrt(sq(a) + sq(b))

extern scalar fused;
extern scalar vof;
extern face vector sf;

typedef struct fragment {
    coord n;
    double alpha;
    double c;  // vof
} fragment;

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


double centroid_point (Point point, scalar vof, coord * boundInter)
{
#if dimension == 3
    coord cellCenter = {x, y, z};
#else
    coord cellCenter = {x, y};
#endif
    coord n = interface_normal (point, vof);
    double alpha = plane_alpha (vof[], n);
    double area = plane_area_center (n, alpha, boundInter);

    foreach_dimension()
        boundInter->x = cellCenter.x + boundInter->x*Delta;
    return area;
}

bool fluid_neighbor (Point point, scalar vof)
{
    // check left and right neighbors
    for(int i = -1; i <= 1; i++)
        if (vof[i] < 0.5)
            return true;

    // check top and bottom neighbors
    for(int j = -1; j <= 1; j++)
        if (vof[0, j] < 0.5)
            return true;

    return false;
}

void fill_fragment (double c, coord n, fragment * frag)
{
    frag->c = c;
    frag->n = n;
    frag->alpha = plane_alpha (c, n);
}

/*
The function below fills frag with the normal vector n, alpha, and the volume fraction of the
cell that is closest to the ghost cell. It also returns the coordinates of the fragment's midpoint 
and cell center.
*/

coord closest_interface (Point point, vector midPoints, scalar vof, vector normals, fragment * frag, coord * fluidCell)
{
    fragment temp_frag;
    coord temp_midPoint, temp_fluidCell = {0,0};
    coord n;
    double min_distance = 1e6;

     for(int i = -1; i <= 1; i++) {
        if (midPoints.x[i] && distance(midPoints.x[i] - x, midPoints.y[i] - y) < min_distance) {
            foreach_dimension()
                temp_midPoint.x = midPoints.x[i];
            // coord n = {normals.x[i], normals.y[i]};
            n.x = normals.x[i];
            n.y = normals.y[i];
            fill_fragment (vof[i], n, &temp_frag);
            temp_fluidCell.x = i*Delta + x;
            temp_fluidCell.y = y;
            min_distance = distance(midPoints.x[i] - x, midPoints.y[i] - y);
        }
     }

     for(int j = -1; j <= 1; j++) {
        double a = midPoints.x[0,j] - x;
        double b = midPoints.y[0,j] - y;
        if (midPoints.x[0,j] && distance(a, b) < min_distance) {
            foreach_dimension()
                temp_midPoint.x = midPoints.x[0,j];
            n.x = normals.x[0,j];
            n.y = normals.y[0,j];   
            // coord n = {normals.x[0,j], normals.y[0,j]};
            fill_fragment (vof[0,j], n, &temp_frag);
            temp_fluidCell.x = x;
            temp_fluidCell.y = j*Delta + y;
            min_distance = distance(a, b);
        }
     }
    // fprintf(stderr, "|| n.x=%g n.y=%g frag.n.x=%g, frag.n.y=%g theta=%g\n", n.x, n.y, temp_frag.n.x, temp_frag.n.y, atan2(temp_frag.n.y, temp_frag.n.x) * 180/M_PI);
    *fluidCell = temp_fluidCell;
    *frag = temp_frag;

    return temp_midPoint;
}

coord boundary_int (Point point, fragment frag, coord fluidCell, scalar vof)
{
    coord boundaryInt;
    coord ghostCell = {x, y};

    double slope = -(frag.n.x / frag.n.y);
    
    double coeff[3];
    if (fabs(frag.n.x) > fabs(frag.n.y) && vof[] == 1) {
        coeff[0] = 1.;
        coeff[1] = slope;
        coeff[2] = sign(frag.n.x)*1.;
    }
    else if (fabs(frag.n.x) < fabs(frag.n.y) && vof[] == 1) {
        coeff[0] = 1./slope;
        coeff[1] = 1;
        coeff[2] = sign(frag.n.y)*1.;
    }
    else {
        coeff[0] = 1./slope;
        coeff[1] = 1;
        coeff[2] = 0.;
    }

    boundaryInt.x = (frag.n.y*coeff[2] - coeff[1]*-frag.alpha);
    boundaryInt.x /= (frag.n.x*coeff[1] - coeff[0]*frag.n.y);

    boundaryInt.y = (-frag.alpha*coeff[0] - coeff[2]*frag.n.x);
    boundaryInt.y /= (frag.n.x*coeff[1] - coeff[0]*frag.n.y);
   
    // Bring point in global coordinate system
    foreach_dimension()
        boundaryInt.x = fluidCell.x + boundaryInt.x*Delta;

    double dist = distance(boundaryInt.x - ghostCell.x, boundaryInt.y - ghostCell.y);

    double max_x, max_y, min_x, min_y;
    foreach_dimension() {
        max_x = fluidCell.x + 0.5*Delta;
        min_x = fluidCell.x - 0.5*Delta;
    }

/*
    if (boundaryInt.x > max_x || boundaryInt.x < min_x || boundaryInt.y > max_y || boundaryInt.y < min_y) {
        fprintf (stderr,"|| ERROR: boundary intercept (%g, %g) is outside of cell domain [%g, %g]x[%g, %g] w/center = (%g,%g)\n",
                    boundaryInt.x, boundaryInt.y, min_x, min_y, max_x, max_y, fluidCell.x, fluidCell.y);
        fprintf (stderr, "|| 1 x=%g y=%g fx=%g fy=%g m=%g d=%g Delta=%g\n", x, y, fluidCell.x, fluidCell.y, slope, dist, Delta);
        fprintf (stderr, "|| 2 n.x=%g n.y=%g alpha=%g c[0]=%g c[1]=%g c[2]=%g\n",
                    frag.n.x, frag.n.y, frag.alpha, coeff[0], coeff[1], coeff[2]);
        fused[] = 1;
    }
*/
    return boundaryInt;
}


coord image_point (coord boundaryInt, coord ghostCell)
{
     double dist = distance(boundaryInt.x - ghostCell.x, boundaryInt.y - ghostCell.y);

     double dx = boundaryInt.x - ghostCell.x;
     double dy = boundaryInt.y - ghostCell.y;

     coord imagePoint = {ghostCell.x + 2*dx, ghostCell.y + 2*dy};

     return imagePoint;
}


coord image_velocity (Point point, vector u, coord imagePoint)
{
    // 2D bilinear interpolation
    coord temp_velo;
    bilinear_interpolation (point, u, imagePoint, &temp_velo);

    return temp_velo;
    /*
    double coeff[6];
    double V[6][6]; // 6x6 Vandermonde matrix
    double values[6];

    for (int i = 0; i < 6; i++) {
        values[i];
        for (int j = 0; j < 6; j++) {

        }
    }
    */
}

#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

foreach_dimension()
static inline double dirichlet_gradient_x (Point point, scalar s, scalar cs,
					   coord n, coord p, coord bc, double * coef)
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
  coord markerCoord, boundaryCondition, dudn;

  centroid_point(point, vof, &markerCoord);
  bilinear_interpolation(point, u, markerCoord, &boundaryCondition);

  foreach_dimension() {
    bool dirichlet = true;
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, vof, n, p, boundaryCondition, &val);
      dudn.x += u.x[]*val;
    }
    else // Neumann
      dudn.x = bc.x;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }

    // fprintf (stderr, "|| x=%g y=%g mC.x=%g mC.y=%g bC.x=%g bC.y=%g dudn.x=%g dudn.y=%g\n",
    //         x, y, markerCoord.x, markerCoord.y, boundaryCondition.x, boundaryCondition.y, dudn.x, dudn.y);

  return dudn;
}


double quadratic_interpolation (scalar uc, coord normalPoint, int type)
{
    double xc[3] = {0}; // x representing x or y depending on the type
    double v[3] = {0};
    double xp = type == 0? normalPoint.y: normalPoint.x;

    foreach_point (normalPoint.x, normalPoint.y) {
        for (int i = -1; i <= 1; i++) {
            if (type == 0)
                xc[i+1] = y + i*Delta;
            else // type = 1
                xc[i+1] = x + i*Delta;
            v[i+1] = type == 0? uc[0,i]: uc[i];
        }
    }

    double interpolate = v[0]*((xp - xc[1])*(xp - xc[2]))/((xc[0] - xc[1])*(xc[0]-xc[2]));
    interpolate += v[1]*((xp - xc[0])*(xp - xc[2]))/((xc[1] - xc[0])*(xc[1] - xc[2]));
    interpolate += v[2]*((xp - xc[0])*(xp - xc[1]))/((xc[2] - xc[0])*(xc[2] - xc[1]));

    return interpolate;
}

double ibm_dirichlet_gradientv2 (Point point, scalar uc, coord n, coord markerCoord, double boundaryCondition)
{
    coord cellCenter = {x, y}, distance, normalPoint;
    double d[2], v[2];
    int type;
    for (int num = 0; num <= 1; num++) {
        foreach_dimension()
            distance.x = ((1 + num) * Delta * n.x) - (markerCoord.x - cellCenter.x);

        if (fabs(n.x) >= fabs(n.y)) {
            type = 0;
            normalPoint.x = cellCenter.x + sign(n.x)*(Delta * (1 + num));
            normalPoint.y = markerCoord.y + distance.y;
        }
        else {  // fabs(n.x) < fabs(n.y)
            type = 1;
            normalPoint.x = markerCoord.x + distance.x;
            normalPoint.y = cellCenter.y + sign(n.y)*(Delta * (1 + num));
        }

        // v[num] = scalar_bilinear_interpolation(point, uc, normalPoint);
        v[num] = quadratic_interpolation (uc, normalPoint, type);
        d[num] = distance(markerCoord.x - normalPoint.x, markerCoord.y - normalPoint.y);
    }
    
    double gradient = ((boundaryCondition - v[0])*(d[1]/d[0]) - (boundaryCondition - v[1])*(d[0]/d[1]));
    gradient /= d[1] - d[0]; 

    return gradient;
}


coord ibm_gradientv2 (Point point, vector u, coord markerCoord, coord n)
{
    coord dudn, boundaryCondition;

    centroid_point (point, vof, &markerCoord);
    bilinear_interpolation (point, u, markerCoord, &boundaryCondition);

    foreach_dimension()
        dudn.x = ibm_dirichlet_gradientv2 (point, u.x, n, markerCoord, boundaryCondition.x);
         
    return dudn;
}


/*
coord boundary_int (Point point, fragment frag, coord fluidCell)
{
    coord boundaryInt;
    coord ghostCell = {x, y};
    coord ghostP[2];
    coord fluidP[2];

    facets (frag.n, frag.alpha, fluidP);
    
    double ghostAlpha = plane_alpha (0.5, frag.n);
    facets (frag.n, ghostAlpha, ghostP);

    
    // Put points in global coordinate system
    for (int i = 0; i < 2; i++)
        foreach_dimension() {
            ghostP[i].x = ghostCell.x + ghostP[i].x*Delta;
            fluidP[i].x = fluidCell.x + fluidP[i].x*Delta;
        }

    double m = (ghostP[1].y - ghostP[0].y) / (ghostP[1].x - ghostP[0].x);  // slope of interface fragment
    double d = 0;  // minimum distance from the ghost cell to the interface
    double theta = atan2(frag.n.y, frag.n.x);

    organize_points(ghostP, fluidP);

    if (fabs(frag.n.y) > fabs(frag.n.x)) {
        d = fabs(ghostP[0].y - fluidP[0].y) / sqrt(1 + sq(m));
    }
    else
        d = fabs(ghostP[0].x - fluidP[0].x) / sqrt(1 + sq(1/m));

    double phi = (M_PI) - fabs(theta);
    if (fabs(phi) > M_PI/2)
        phi -= M_PI/2;
    phi *= sign(theta);
    if (phi < 0)
        phi = fabs(phi) + M_PI;

    boundaryInt.x = d * cos(phi) + ghostCell.x;
    boundaryInt.y = d * sin(phi) + ghostCell.y;

    double maxX = fluidCell.x + 0.5*Delta;
    double minX = fluidCell.x - 0.5*Delta;
    double maxY = fluidCell.y + 0.5*Delta;
    double minY = fluidCell.y - 0.5*Delta;

    if (boundaryInt.x > maxX || boundaryInt.x < minX || boundaryInt.y > maxY || boundaryInt.y < minY) {
        fprintf (stderr,"|| ERROR: boundary intercept (%g, %g) is outside of cell domain [%g, %g]x[%g, %g] w/center = (%g,%g)\n",
                    boundaryInt.x, boundaryInt.y, minX, minY, maxX, maxY, fluidCell.x, fluidCell.y);

    fprintf (stderr, "|| 1 x=%g y=%g fx=%g fy=%g m=%g d=%g Delta=%g\n", x, y, fluidCell.x, fluidCell.y, m, d, Delta);
    fprintf (stderr, "|| 2 g1.x=%g g1.y=%g g2.x=%g g2.y=%g f1.x=%g f1.y=%g f2.x=%g f2.y=%g\n", ghostP[0].x, ghostP[0].y, ghostP[1].x, ghostP[1].y, fluidP[0].x, fluidP[0].y, fluidP[1].x, fluidP[1].y);
    fprintf (stderr, "|| 3 theta=%g phi=%g n.x=%g n.y=%g galpha=%g falpha=%g\n", theta*180/M_PI, phi*180/M_PI, frag.n.x, frag.n.y, ghostAlpha, frag.alpha);

                    fused[] = 1;
    }
    return boundaryInt;
}
*/
