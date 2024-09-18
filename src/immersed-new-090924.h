#include "fractions.h"
#include "ibm-utils.h"

extern coord vc;          // object's imposed velocity
extern scalar vof;
extern face vector sf;
extern int maxlevel;

vector desiredForce[];    // force calculated at the marker point
vector cellForce[];       // force located at the cell center (after spreading)
vector markerCoord[];     // field to store the coordinates of all marker points

vector velocityGrad[];

face vector faceForce[];  // for averaging the cell force to get the face values

vector utemp[];
vector forceTotal[];
vector totalDesired[];

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
            utemp.x[] = u.x[] - dt*(p[] - p[-1])/(Delta);
            // utemp.x[] = u.x[];
            forceTotal.x[] = 0.;
            totalDesired.x[] = 0.;
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
                    double delta_u = delta_func (x, y, markerPoint.x, markerPoint.y, Delta);
                    foreach_dimension()
                        markerVelocity.x += utemp.x[] * delta_u * dv();
                }

                // calculate the desired force at the marker point
                desiredVelocity = vc;

                foreach_dimension() {
                    desiredForce.x[] = (desiredVelocity.x - markerVelocity.x) / dt;
                    markerCoord.x[] = markerPoint.x;
                }
            }
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
            else
                foreach_dimension()
                    desiredForce.x[] = markerCoord.x[] = 0.;
        }

        scalar hd[];
        foreach() {
            if (markerCoord.x[]) {
                double hd_sum = 0.;
                coord markerPoint = {markerCoord.x[], markerCoord.y[]};
                foreach_neighbor()
                    if (vof[] >= 0)
                        hd_sum += delta_func (x, y, markerPoint.x, markerPoint.y, Delta) * dv();
                hd[] = hd_sum;
            }
            else
                hd[] = 0.;
        }

        // 2. spread the force at the marker point to the nearby cell centers
        foreach() {
            coord forceSum = {0};
            if (level == maxlevel && vof[] >= 0) {
                double x1 = x, y1 = y;

                foreach_neighbor()
                    if (markerCoord.x[] && level == maxlevel) {
                        double delta_h = delta_func (x1, y1, markerCoord.x[], markerCoord.y[], Delta);
                        foreach_dimension() {
                            forceSum.x += (desiredForce.x[] * delta_h * dv()) / hd[];
                        }
                    }
            }
            foreach_dimension()
                cellForce.x[] = forceSum.x;
        }

        foreach()
            foreach_dimension() {
                totalDesired.x[] += desiredForce.x[];
                forceTotal.x[] += cellForce.x[];
                utemp.x[] += dt*cellForce.x[];
            }
    }
    
    // 3. correct interfacial velocity
    
    foreach_face()
        faceForce.x[] = fm.x[]*(face_value (forceTotal.x, 0));
    a = faceForce;
    
    /*
    correction(-dt);
    foreach()
        foreach_dimension()
            g.x[] += forceTotal.x[];
    correction(dt);
   */
}


//  g is used to find uf t+dt/2 at the next time step, so the contributions
//  from f should be subtracted

event end_timestep (i++)
{
    // trash({a});
    // centered_gradient (p, g);

    trash ({velocityGrad});
    foreach()
        foreach_dimension()
            velocityGrad.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
}


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
            coord marker = {markerCoord.x[], markerCoord.y[]};
            // Fn = area*scalar_bilinear_interpolation (point, p, marker);
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
            // coord dudn = ibm_gradient (point, u, b, n);
            // coord dudn = ibm_gradientv2 (point, u, b, n);
            // coord dudn = extrapolate_gradient (point, vof, marker, n, velocityGrad);
        	foreach_dimension()
                Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
            }
        }
    }
    *Fp = Fps; *Fmu = Fmus;
}

void immersed_forcev2 (scalar c, scalar p, vector u,
			face vector mu, coord * Fp, coord * Fmu)
{
    coord Fps = {0}, Fmus = {0};

    foreach (reduction(+:Fps) reduction(+:Fmus), nowarning) {

        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      
            coord n, b;
            double area = embed_geometry (point, &b, &n);
            area *= pow (Delta, dimension - 1);

            double Fn = area*embed_interpolate (point, p, b);
            coord marker = {markerCoord.x[], markerCoord.y[]};
            // Fn = area*scalar_bilinear_interpolation (point, p, marker);
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
            // coord dudn = embed_gradient (point, u, b, n, vc);
            coord dudn = ibm_gradient (point, u, b, n);
            // coord dudn = ibm_gradientv2 (point, u, b, n);
            // coord dudn = extrapolate_gradient (point, vof, marker, n, velocityGrad);
        	foreach_dimension()
                Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
            }
        }
    }
    *Fp = Fps; *Fmu = Fmus;
}


vector frictionDrag[];
vector pressureDrag[];

void immersed_forcev3 (scalar c, scalar p, vector u,
			face vector mu, coord * Fp, coord * Fmu)
{
    coord Fps = {0}, Fmus = {0};

     foreach (reduction(+:Fps) reduction(+:Fmus), nowarning) {

        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      
            coord n, b;
            double area = embed_geometry (point, &b, &n);
            area *= pow (Delta, dimension - 1);

            double Fn = area*embed_interpolate (point, p, b);
            coord marker = {markerCoord.x[], markerCoord.y[]};
            // Fn = area*scalar_bilinear_interpolation (point, p, marker);
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

void immersed_forcev4 (scalar c, scalar p, vector u,
			face vector mu, coord * Fp, coord * Fmu)
{
    coord Fps = {0}, Fmus = {0};

    foreach (reduction(+:Fps) reduction(+:Fmus), nowarning) {

        if (c[] > 1e-6 && c[] < 1. - 1e-6) {
      
            coord n, b;
            double area = embed_geometry (point, &b, &n);
            area *= pow (Delta, dimension - 1);

            double Fn = area*embed_interpolate (point, p, b);
            coord marker = {markerCoord.x[], markerCoord.y[]};
            // Fn = area*scalar_bilinear_interpolation (point, p, marker);
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
            // coord dudn = embed_gradient (point, u, b, n, vc);
            // coord dudn = ibm_gradient (point, u, b, n);
            // coord dudn = ibm_gradientv2 (point, u, b, n);
            coord dudn = extrapolate_gradient (point, vof, marker, n, velocityGrad);
        	foreach_dimension()
                Fmus.x -= area*mua*(dudn.x*(sq(n.x) + 1.) + dudn.y*n.x*n.y);
            }
        }
    }
    *Fp = Fps; *Fmu = Fmus;
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

