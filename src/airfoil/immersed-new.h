#ifndef BASILISK_HEADER_19
#define BASILISK_HEADER_19
#line 1 "./../immersed-new.h"
#include "fractions.h"
#include "ibm-utils.h"

extern coord vc;
extern scalar vof;
extern face vector sf;
extern int maxlevel;


// #if ROTATION
extern coord centerRot;   // to store the center of rotation cooridinates, if applicable
extern double w0;         // angular velocity (rad/s)
// #endif // no rotation

vector desiredForce[];    // force calculated at the marker point
vector cellForce[];       // force located at the cell center (after spreading)
vector markerCoord[];     // field to store the coordinates of all marker points

face vector faceForce[];  // for averaging the cell force to get the face values

event acceleration (i++)
{
    trash({cellForce, desiredForce, markerCoord, faceForce});
    
    // temporary velocity field with advection (n), viscous (n), and pressure (n-1)
    vector utemp[];
    foreach()
        foreach_dimension()
            utemp.x[] = u.x[] - dt*(alpha.x[]*(p[] - p[-1])/Delta);
    
    // 1. calculate the force at the marker point
    foreach() {
        coord markerVelocity = {0}, desiredVelocity, markerPoint;
        if (vof[] > 0 && vof[] < 1) {
            marker_point (point, vof, &markerPoint);

            // interpolate to find velocity at marker point
            foreach_neighbor() {
                double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta);
                foreach_dimension()
                    markerVelocity.x += utemp.x[] * delta_u * dv();
            }

            // calculate the desired force at the marker point
            // #if ROTATION
            double xr = centerRot.x - markerPoint.x, yr = centerRot.y - markerPoint.y;
            double theta = abs(atan(yr/xr));
            double r = distance(xr, yr);

            desiredVelocity.x = (r * w0)*sin(theta)*-sign(yr);
            desiredVelocity.y = (r * w0)*cos(theta)*sign(xr);

            // #else // no rotation

            // desiredVelocity = vc;
            // #endif

            foreach_dimension() {
                desiredForce.x[] = (desiredVelocity.x - markerVelocity.x) / dt;
                markerCoord.x[] = markerPoint.x;
            }
        }
        /*
        else if (empty_neighbor(point, &markerPoint, vof)) {
            foreach_neighbor() {
                double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta);
                foreach_dimension()
                    markerVelocity.x += utemp.x[] * delta_u * dv();
            }
            foreach_dimension() {
                desiredForce.x[] = (vc.x - markerVelocity.x) / dt;
                markerCoord.x[] = markerPoint.x;
            }
        }
        */
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
    foreach_face()
        faceForce.x[] = fm.x[]*(face_value (cellForce.x, 0));
    a = faceForce;
   
    // 4. mask inside velocity
    /*
    foreach()
        if (vof[] >= 0.5)
            foreach_dimension()
                u.x[] = vc.x;
                */
  
}


//  g is used to find uf t+dt/2 at the next time step, so the contributions
//  from f should be subtracted

event end_timestep (i++)
{
    coord desiredVelocity;
    trash(a);
    centered_gradient (p, g);
    foreach()
        if (vof[] >= 0.5) {
            // #if ROTATION 
            double xr = centerRot.x - x, yr = centerRot.y - y;
            double theta = fabs(atan(yr/xr));
            double r = distance(xr, yr);
            desiredVelocity.x = (r * w0)*sin(theta)*-sign(yr);
            desiredVelocity.y = (r * w0)*cos(theta)*sign(xr);


            // #else // no rotation

            // desiredVelocity = vc;
            // #endif
            foreach_dimension()
                u.x[] = desiredVelocity.x;
        }
}


void immersed_force (scalar c, scalar p, vector u,
			face vector mu, coord * Fp, coord * Fmu)
{
    double Fn = 0.;
    coord Fps = {0}, Fmus = {0};

    foreach () {

        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      
            coord n, b;
            double area = embed_geometry (point, &b, &n);
            area *= pow (Delta, dimension - 1);

            Fn = area*embed_interpolate (point, p, b);
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
            coord dudn = embed_gradient (point, u, b, n, vc);
        	foreach_dimension()
                Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
            }
        }
    }
    *Fp = Fps; *Fmu = Fmus;
}


#endif
