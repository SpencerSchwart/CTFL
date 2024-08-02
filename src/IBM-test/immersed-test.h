#ifndef BASILISK_HEADER_20
#define BASILISK_HEADER_20
#line 1 "./../immersed-test.h"
#include "fractions.h"
#include "ibm-utils.h"

extern coord vc;
extern scalar vof;
extern face vector sf;
extern int maxlevel;

vector cellForce[];
vector markerForce[];
vector desiredForce[];
vector markerCoord[];

event end_timestep (i++)
{
    trash({cellForce, markerForce, desiredForce, markerCoord});

    // 1. calculate the force at the marker point
    foreach() {
        coord markerVelocity = {0}, markerPoint;
        if (vof[] > 0 && vof[] < 1) {
            marker_point (point, vof, &markerPoint);

            // interpolate to find velocity at marker point
            foreach_neighbor() {
                double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta);
                foreach_dimension()
                    markerVelocity.x += u.x[] * delta_u * dv();
            }

            // calculate the desired force at the marker point
            foreach_dimension() {
                desiredForce.x[] = (vc.x - markerVelocity.x) / dt;
                markerCoord.x[] = markerPoint.x;
            }
        }
        else if (empty_neighbor(point, &markerPoint, vof)) {
            foreach_neighbor() {
                double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta);
                foreach_dimension()
                    markerVelocity.x += u.x[] * delta_u * dv();
            }
            foreach_dimension() {
                desiredForce.x[] = (vc.x - markerVelocity.x) / dt;
                markerCoord.x[] = markerPoint.x;
            }
        }
        else
            foreach_dimension()
                desiredForce.x[] = markerCoord.x[] = 0.;
    }

    // 2. spread the force at the marker point to the nearby cell centers
    foreach() {
        coord forceSum = {0};
        if (level == maxlevel) {
            double x1 = x, y1 = y;
            foreach_neighbor()
                if (markerCoord.x[] && level == maxlevel) {
                    double delta_h = delta_func (x1, y1, markerCoord.x[], markerCoord.y[], Delta);
                    foreach_dimension()
                        forceSum.x += desiredForce.x[] * delta_h * dv();
                }
        }
        foreach_dimension()
            cellForce.x[] = forceSum.x;
    }

    // 3. correct interfacial velocity
    correction(-dt);
    foreach()
        foreach_dimension()
            g.x[] += (cellForce.x[]);
    correction(dt);

    // 4. mask inside velocity
    foreach()
        if (vof[] >= 0.5)
            foreach_dimension()
                u.x[] = vc.x;
}

double immersed_interpolate (Point point, scalar s, coord p)
{
  assert (dimension == 2);
  int i = sign(p.x), j = sign(p.y);
  if (vof[i] && vof[0,j] && vof[i,j])
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

/*
void ibm_force (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
    coord Fps = {0}, Fmus = {0};
    foreach (reduction(+:Fps) reduction (+:Fmus), nowarning)
        if (vof[] > 0. && vof[] < 1.) {
            coord cellCenter = {x,y}, pm;
            coord n = interface_normal (point, vof), pc;
            double alpha = plane_alpha (vof[], n);
            double area = pow(Delta, dimension - 1) * plane_area_center (n, alpha, &pc);
            foreach_dimension()
                pm.x = cellCenter.x + pc.x*Delta;
            // double Fn = area*scalar_bilinear_interpolation (point, p, pm);
            double Fn = area*embed_interpolate (point, p, pc);
            foreach_dimension()
                Fps.x += Fn*n.x;
            if (constant(mu.x) != 0.) {
                double mua = 0., fa = 0.;
                foreach_dimension() {
                mua += mu.x[] + mu.x[1];
                fa += fm.x[] + fm.x[1];
                }
                mua /= fa;
                assert (dimension == 2);
                coord dudn = ibm_gradient (point, u, pc, n, vc, vof);
                foreach_dimension()
                    Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
            }

        }
    *Fp = Fps; *Fmu = Fmus;
}
*/

#endif
