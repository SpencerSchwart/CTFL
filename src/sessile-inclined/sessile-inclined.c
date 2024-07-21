#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "../contact-embed.h"

double theta0;

int main()
{
  size (2.);
  origin (-1, 0);

  mu1 = mu2 = 0.1;

  f.sigma = 1.;

  for (theta0 = 15; theta0 <= 165; theta0 += 15) {
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
}

event init (t = 0)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (y - x - 0.97);
  boundary ({phi});
  fractions (phi, cs, fs);
  fraction (f, - (sq(x - 0) + sq(y - 1.) - sq(0.25)));
}

event logfile (i++; t <= 20)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (y - x - 0.97);
  boundary ({phi});
  fractions (phi, cs, fs);
  fraction (f, - (sq(x - 0) + sq(y - 1.) - sq(0.25)));
}

double equivalent_contact_angle (double R, double V)
{
  double x0 = 0., x1 = pi;
  while (x1 - x0 > 1e-4) {
    double x = (x1 + x0)/2.;
    double f = V - sq(R)*(x - sin(x)*cos(x));
    if (f > 0.)
      x0 = x;
    else
      x1 = x;
  }
  return (x0 + x1)/2.;
}
  
event end (t = end)
{
  output_facets (f, stdout);
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
    
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g %.4g\n", N, theta0, R/sqrt(V/pi), s.stddev,
	   equivalent_contact_angle (R, V)*180./pi);
}
