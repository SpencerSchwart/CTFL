extern scalar vof;
extern face vector sf;

void boundary_cells (scalar c, scalar ink) {

  scalar phi[];

  foreach() {
    if(c[] > 1e-6 && c[] < 1. - 1e-6)
      phi[] = 1.;
    else
      phi[] = -1.;
  }
  boundary ({phi});
  fractions (phi, ink);
}


#define quadratic(x,a1,a2,a3) \
  (((a1)*((x) - 1.) + (a3)*((x) + 1.))*(x)/2. - (a2)*((x) - 1.)*((x) + 1.))

void printx (int k, int l, Point point, coord n, coord p, int i, int j, double d[], double y1, double v[], scalar vof)
{
 //   fprintf (stderr, "|| %d %d %g %g %g n.x=%g n.y=%g p.x=%g p.y=%g i=%d j=%d d[%d]=%g y1=%g v[%d]=%g\n",
 //            k, l , x, y, vof[], n.x, n.y, p.x, p.y, i, j, l, d[l], y1, l, v[l]);
}

void printy (int k, Point point, coord n, coord p, double gradient, scalar vof, double d[], double v[], double bc)
{
 //   fprintf (stderr, "|| %d %g %g %g n.x=%g n.y=%g p.x=%g p.y=%g grad=%g\n",
 //            k, x, y, vof[], n.x, n.y, p.x, p.y, gradient);
}



foreach_dimension()
static inline double dirichlet_gradiento_x (Point point, scalar s, scalar cs,
					   coord n, coord p, double bc,
					   double * coef)
{
  double d[2] = {0,0}, v[2] = {nodata,nodata};
  bool defined = true;
  // foreach_dimension()
    // if (defined && sf.x[(n.x > 0.)]) // what to do here?
    //  defined = false;
  if (defined)
    for (int l = 0; l <= 1; l++) {
      int i = (l + 1)*sign(n.x);
      d[l] = (i - p.x)/n.x;
      double y1 = p.y + d[l]*n.y;
      int j = y1 > 0.5 ? 1 : y1 < -0.5 ? -1 : 0;
      y1 -= j;
#if dimension == 2
   if (cs[i,j-1] < 0.5 && cs[i,j] < 0.5 && cs[i,j+1] < 0.5) {
	    v[l] = quadratic (y1, (s[i,j-1]), (s[i,j]), (s[i,j+1]));
        printx (1, l, point, n, p, i, j, d, y1, v, vof);
    }
#else // dimension == 3
      double z = p.z + d[l]*n.z;
      int k = z > 0.5 ? 1 : z < -0.5 ? -1 : 0;
      z -= k;
      bool defined = sf.x[i + (i < 0),j,k];
      for (int m = -1; m <= 1 && defined; m++)
	if (!sf.y[i,j,k+m] || !sf.y[i,j+1,k+m] ||
	    !sf.z[i,j+m,k] || !sf.z[i,j+m,k+1] ||
	    !cs[i,j+m,k-1] || !cs[i,j+m,k] || !cs[i,j+m,k+1])
	  defined = false;
      if (defined)
	// bi-quadratic interpolation
	v[l] =
	  quadratic (z,
		     quadratic (y1,
				(s[i,j-1,k-1]), (s[i,j,k-1]), (s[i,j+1,k-1])),
		     quadratic (y1,
				(s[i,j-1,k]),   (s[i,j,k]),   (s[i,j+1,k])),
		     quadratic (y1,
				(s[i,j-1,k+1]), (s[i,j,k+1]), (s[i,j+1,k+1])));
#endif // dimension == 3
      else {
//        fprintf(stderr, "bye\n");
	    break;
        }
    }
  if (v[0] == nodata) {
    d[0] = max(1e-3, fabs(p.x/n.x));
    *coef = - 1./(d[0]*Delta);
    // fprintf (stderr, "yes0\n");
    return bc/(d[0]*Delta);
  }
  *coef = 0.;
  double gradient = 0;
  if (v[1] != nodata) { // third-order gradient 
    // fprintf (stderr,"|| v[0]=%g v[1]=%g d[0]=%g d[1]=%g delta=%g bc=%g\n", v[0], v[1], d[0], d[1], Delta, bc);
    gradient = (d[1]*(bc - v[0])/d[0] - d[0]*(bc - v[1])/d[1])/((d[1] - d[0])*Delta);
    }
  else {
    // fprintf (stderr, "yes2\n");
    gradient = (bc - v[0])/(d[0]*Delta); // second-order gradient
    }
  printy (2, point, n, p, gradient, vof, d, v, bc);
  return gradient;
}


double dirichlet_gradiento (Point point, scalar s, scalar cs,
			   coord n, coord p, double bc, double * coef)
{
#if dimension == 2
  foreach_dimension()
    if (fabs(n.x) >= fabs(n.y))
      return dirichlet_gradiento_x (point, s, cs, n, p, bc, coef);
#else // dimension == 3
  if (fabs(n.x) >= fabs(n.y)) {
    if (fabs(n.x) >= fabs(n.z))
      return dirichlet_gradiento_x (point, s, cs, n, p, bc, coef);
  }
  else if (fabs(n.y) >= fabs(n.z))
    return dirichlet_gradiento_y (point, s, cs, n, p, bc, coef);
  return dirichlet_gradiento_z (point, s, cs, n, p, bc, coef);
#endif // dimension == 3
  return nodata;
}


double embed_geometryo (Point point, coord * b, coord * n)
{
  *n = facet_normal (point, vof, sf);
  double alpha = plane_alpha (vof[], *n);
  double area = plane_area_center (*n, alpha, b);
  normalize (n);
  return area;
}

double embed_interpolateo (Point point, scalar s, coord p)
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
coord embed_gradiento (Point point, vector u, coord p, coord n)
{
  coord dudn;
  foreach_dimension() {
    bool dirichlet = true;
    double vb = 0.;
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradiento (point, u.x, vof, n, p, vb, &val);
    }
    else // Neumann
      dudn.x = vb;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}


void interface_force (scalar c, scalar p, vector u,
			face vector mu, coord * Fp, coord * Fmu, double t)
{
//  FILE * fp = fopen ("gradient.dat", "w");
  double Fn = 0.;
  coord Fps = {0}, Fmus = {0};
  int counter = 0;
  foreach () {
    if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      
      coord n, b;
      double area = embed_geometryo (point, &b, &n);
      area *= pow (Delta, dimension - 1);

      Fn = area*embed_interpolateo (point, p, b);
      foreach_dimension()
	    Fps.x -= Fn*n.x;

      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fm.x[] + fm.x[1];
	}
	mua /= fa;
	assert (dimension == 2);
    counter++;
	coord dudn = embed_gradiento (point, u, b, n);
    // if (t >= 39.9) 
    //    fprintf (stderr, "%d %g %g %g %g %g\n", counter, x, y, vof[], dudn.x, dudn.y);
	foreach_dimension()
	  Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
      }
    }
  }
  *Fp = Fps; *Fmu = Fmus;

}

double embed_vorticityo (Point point, vector u, coord p, coord n)
{
  coord dudn = embed_gradiento (point, u, p, n);

  return dudn.y*n.x - dudn.x*n.y;
}

