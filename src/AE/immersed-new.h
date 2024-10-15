#ifndef BASILISK_HEADER_19
#define BASILISK_HEADER_19
#line 1 "./../immersed-new.h"
#include "fractions.h"
#include "ibm-utils.h"

extern coord vc;          // object's imposed velocity
extern scalar vof;
extern face vector sf;
extern int maxlevel;
// #define ROTATION 1

#ifdef ROTATION
extern coord centerRot;
extern double w0;
#endif

vector desiredForce[];    // force calculated at the marker point
vector cellForce[];       // force located at the cell center (after spreading)
vector markerCoord[];     // field to store the coordinates of all marker points

vector velocityGrad[];

face vector faceForce[];  // for averaging the cell force to get the face values

vector utemp[];
vector forceTotal[];

event init (t = 0)
{
    foreach()
        foreach_dimension() { 
            velocityGrad.x[] = 0.;
                if (vof[] == 1)
                    u.x[] = vc.x;
        }
}


event acceleration (i++)
{
    trash({cellForce, desiredForce, markerCoord, faceForce, utemp, forceTotal});

    foreach() {
        foreach_dimension() {
            utemp.x[] = u.x[] - dt*(p[] - p[-1])/(Delta); // second order central difference?
            forceTotal.x[] = 0.;
        }
    }

    for (int counter = 0; counter < 10; counter++) { 

        // 1. calculate the force at the marker point
        foreach() {

            coord markerVelocity = {0}, desiredVelocity, markerPoint;

            if (vof[] > 0 && vof[] < 1) {

                marker_point (point, vof, &markerPoint);

                // interpolate to find velocity at marker point
                foreach_neighbor() {
#if dimension == 3                
                    double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta, z, markerPoint.z);
#else
                    double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta);
#endif                    
                    foreach_dimension()
                        markerVelocity.x += utemp.x[] * delta_u * dv();
                }

                // calculate the desired force at the marker point
#ifdef ROTATION
                desiredVelocity.x = w0 * (markerPoint.y - centerRot.y);
                desiredVelocity.y = -w0 * (markerPoint.x - centerRot.x);
#else
                desiredVelocity = vc;
#endif

                foreach_dimension() {
                    desiredForce.x[] = (desiredVelocity.x - markerVelocity.x) / dt;
                    markerCoord.x[] = markerPoint.x;
                }
            }
            else if (empty_neighbor(point, &markerPoint, vof)) {
                foreach_neighbor() {
#if dimension == 3                
                    double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta, z, markerPoint.z);
#else
                    double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta);
#endif                    
                    foreach_dimension()
                        markerVelocity.x += utemp.x[] * delta_u * dv();
                }
                coord desiredVelocity = vc;
#ifdef ROTATION
                desiredVelocity.x = w0 * (markerPoint.y - centerRot.y);
                desiredVelocity.y = -w0 * (markerPoint.x - centerRot.x);
#endif
                foreach_dimension() {
                    desiredForce.x[] = (desiredVelocity.x - markerVelocity.x) / dt;
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
#if dimension == 3
                double z1 = z;
#endif
                  foreach_neighbor()
                    if (markerCoord.x[] && level == maxlevel) {
#if dimension == 3                
                        double delta_h = delta_func (x1, y1, markerCoord.x[], markerCoord.y[], Delta, z1, markerCoord.z[]);
#else
                        double delta_h = delta_func (x1, y1, markerCoord.x[], markerCoord.y[], Delta);
#endif                    
                        foreach_dimension() {
                            forceSum.x += (desiredForce.x[] * delta_h * dv());
                        }
                        // fprintf (stderr, "|| %g %g %g %g delta=%g F.x=%g F.y=%g F.z=%g px=%g py=%g pz=%g\n", vof[], x1, y1, z1, delta_h, forceSum.x, forceSum.y, forceSum.z, markerCoord.x[], markerCoord.y[], markerCoord.z[]);
                    }
            // fprintf (stderr, "|| DONE: %g %g %g sum.x=%g sumy=%g sumz=%g\n", x, y, z, forceSum.x, forceSum.y, forceSum.z);
            }
            foreach_dimension() 
                cellForce.x[] = forceSum.x;

        }

        foreach()
            foreach_dimension() {
                forceTotal.x[] += cellForce.x[];
                utemp.x[] += dt*cellForce.x[];
            }
    }
    
    // 3. correct interfacial velocity
    foreach_face()
        faceForce.x[] = fm.x[]*(face_value (forceTotal.x, 0));
    a = faceForce;

}


//  g is used to find uf t+dt/2 at the next time step, so the contributions
//  from f should be subtracted

event end_timestep (i++)
{
    trash({a});
    centered_gradient (p, g);

    trash ({velocityGrad});
    foreach()
        foreach_dimension()
            velocityGrad.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
}


vector frictionDrag[];
vector pressureDrag[];

void immersed_force (scalar c, scalar p, vector u,
			face vector mu, coord * Fp, coord * Fmu)
{
    coord Fps = {0}, Fmus = {0};

     foreach (reduction(+:Fps) reduction(+:Fmus), nowarning) {

        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      
            coord n, b;
            double area = embed_geometry (point, &b, &n);
            area *= pow (Delta, dimension - 1);

            double Fn = area*embed_interpolate (point, p, b);
            foreach_dimension() {
	            Fps.x -= Fn*n.x;
                pressureDrag.x[] = -Fn*n.x;
            }

            if (constant(mu.x) != 0.) {
	            double mua = 0., fa = 0.;
                foreach_dimension() {
                    mua += mu.x[] + mu.x[1];
                    fa  += fm.x[] + fm.x[1];
                }
                mua /= fa;
            	assert (dimension == 2);

                coord dudn = ibm_gradientv2 (point, u, b, n);
            	foreach_dimension() {
                    Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
                    frictionDrag.x[] = -(area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y));
                }
            }
        }
    }
    *Fp = Fps; *Fmu = Fmus;
}

coord ibm_force ()
{
    coord ibmForce = {0};
    foreach(reduction(+:ibmForce))
        foreach_dimension()
            ibmForce.x += forceTotal.x[]*dv();
    return ibmForce;
}

double ibm_vorticity (Point point, vector u, coord p, coord n)
{
    coord dudn = ibm_gradientv2 (point, u, p, n);

    return dudn.y*n.x - dudn.x*n.y;
}

double embed_vorticity (Point point, vector u, coord p, coord n)
{
    coord dudn = embed_gradient (point, u, p, n, vc);

    return dudn.y*n.x - dudn.x*n.y;
}

double ibm_pressure (Point point, scalar p, coord markerPoint)
{
    double distance_min = 1e10;
    double pressure = 0.;
    foreach_neighbor() {
        if (vof[] == 0) {
            double d = distance (x - markerPoint.x, y - markerPoint.y);
            if (d <= distance_min) {
                distance_min = d;
                pressure = p[];
            }
        }
    }
    return pressure;
}

double ibm_pressurev2 (Point point, scalar p, coord markerPoint)
{
    double distance_min = 1e10;
    double pressure = 0.;
    foreach_neighbor() {
        if (vof[] < 0.5) {
            double d = distance (x - markerPoint.x, y - markerPoint.y);
            if (d <= distance_min) {
                distance_min = d;
                pressure = p[];
            }
        }
    }
    return pressure;
}

double ibm_pressurev3 (Point point, scalar p, coord markerPoint)
{
    double pressure = 0;
    foreach_neighbor() {
        if (vof[] < 0.5) {
            double delta_p = delta_func (x, y, markerPoint.x, markerPoint.y, Delta);
            foreach_dimension()
               pressure += p[] * delta_p * dv();
        }
    }

    return pressure;
}

#endif
