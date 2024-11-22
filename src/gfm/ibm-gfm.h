#ifndef BASILISK_HEADER_19
#define BASILISK_HEADER_19
#line 1 "./../ibm-gfm.h"
#include "fractions.h"
#include "ibm-gfm-utils.h"

extern scalar vof;
extern face vector sf;
extern coord vc;

scalar fused[];
scalar id[];
scalar q[];         // mass source/sink term

vector normals[];
vector midPoints[];

vector utemp[];

/*
event acceleration (i++)
{
    // 1. Initalize fields to hold interface normals and fragment midpoints

    foreach() {
        id[] = 0;  // -1 instead?
        if (vof[] < 1 && vof[] > 0) {
            coord midPoint;
            centroid_point (point, vof, &midPoint);
            coord n = interface_normal (point, vof);
            foreach_dimension() {
                midPoints.x[] = midPoint.x;
                normals.x[] = n.x;
            }
        }
        else
            foreach_dimension() {
                midPoints.x[] = 0;
                normals.x[] = 0;
            }
    }
    
    // 2. Find the total # of image points/ghost cells (to set list size)

    int numCells = 0;
    foreach() {
        if (fluid_neighbor(point, vof) && vof[] >= 0.5)
            numCells += 1;
    }
    
    coord imagePoints[numCells];
    coord intercepts[numCells];
    
    // 3. Find ghost cell, boundary intercepts, and image points

    int count = 0;
    foreach() {
         if (fluid_neighbor(point, vof) && vof[] >= 0.5) {
            fragment interFrag;
            coord fluidCell, ghostCell = {x,y};

            coord closestCell = closest_interface (point, midPoints, vof, normals, &interFrag, &fluidCell);
            coord boundaryIntercept = boundary_int (point, interFrag, fluidCell, vof);
            coord imagePoint = image_point (boundaryIntercept, ghostCell);

            count += 1;
            id[] = count;

            imagePoints[count-1] = imagePoint;
            intercepts[count-1] = boundaryIntercept;
            
            fprintf(stderr, "|| ip.x=%g ip.y=%g bi.x=%g bi.y=%g\n",
                    imagePoints[count-1].x, imagePoints[count-1].y,
                    intercepts[count-1].x, intercepts[count-1].y);
            
        }
    }

    // 4. Calculate image point values

    coord imageVelos[numCells];

    for (int i = 0; i < numCells; i++) {
        foreach_point(imagePoints[i].x, imagePoints[i].y) {
            imageVelos[i] = image_velocity (point, u, imagePoints[i]); 

            fprintf(stderr,"|| ip.x=%g ip.y=%g iv.x=%g iv.y=%g\n",
                            imagePoints[i].x, imagePoints[i].y,
                            imageVelos[i].x, imageVelos[i].y);
        }
    }
    
    // 5. Calculate and assign velocity at ghost cell

    foreach() {
        if (id[] > 0) {
            int index = id[];
            foreach_dimension() {
                u.x[] = -1*imageVelos[index-1].x;
                // u.x[] = -1*imageVelos[index-1].x;

            }
        }
        else if (vof[] == 1) {
            foreach_dimension()
                u.x[] = 0.;
        }
        foreach_dimension()
            utemp.x[] = u.x[];
    }
    boundary({u});
}
*/

event advection_term (i++)
{
    for (scalar s in {p, u, g, pf})
        foreach()
            if (vof[] >= 1)
	            s[] = 0.;
}

event end_timestep (i++)
{
    // 1. Initalize fields to hold interface normals and fragment midpoints

    foreach() {
        id[] = 0;  // -1 instead?
        if (vof[] < 1 && vof[] > 0) {
            coord midPoint;
            centroid_point (point, vof, &midPoint);
            coord n = interface_normal (point, vof);
            foreach_dimension() {
                midPoints.x[] = midPoint.x;
                normals.x[] = n.x;
            }
        }
        else
            foreach_dimension() {
                midPoints.x[] = 0;
                normals.x[] = 0;
            }
    }
    
    // 2. Find the total # of image points/ghost cells (to set list size)

    int numCells = 0;
    foreach() {
        if (fluid_neighbor(point, vof) && vof[] >= 0.5)
            numCells += 1;
    }
    
    coord imagePoints[numCells];
    coord intercepts[numCells];
    
    // 3. Find ghost cell, boundary intercepts, and image points

    int count = 0;
    foreach() {
         if (fluid_neighbor(point, vof) && vof[] >= 0.5) {
            fragment interFrag;
            coord fluidCell, ghostCell = {x,y};

            coord closestCell = closest_interface (point, midPoints, vof, normals, &interFrag, &fluidCell);
            coord boundaryIntercept = boundary_int (point, interFrag, fluidCell, vof);
            coord imagePoint = image_point (boundaryIntercept, ghostCell);

            count += 1;
            id[] = count;

            imagePoints[count-1] = imagePoint;
            intercepts[count-1] = boundaryIntercept;
/*            
            fprintf(stderr, "|| ip.x=%g ip.y=%g bi.x=%g bi.y=%g\n",
                    imagePoints[count-1].x, imagePoints[count-1].y,
                    intercepts[count-1].x, intercepts[count-1].y);
*/            
        }
    }

    // 4. Calculate image point values

    coord imageVelos[numCells];

    for (int i = 0; i < numCells; i++) {
        foreach_point(imagePoints[i].x, imagePoints[i].y) {
            imageVelos[i] = image_velocity (point, u, imagePoints[i]); 
/*
            fprintf(stderr,"|| ip.x=%g ip.y=%g iv.x=%g iv.y=%g\n",
                            imagePoints[i].x, imagePoints[i].y,
                            imageVelos[i].x, imageVelos[i].y);
*/
        }
    }
    
    // 5. Calculate and assign velocity at ghost cell

    foreach() {
        if (id[] > 0) {
            int index = id[];
            foreach_dimension() {
                u.x[] = -1*imageVelos[index-1].x;
                // u.x[] = -1*imageVelos[index-1].x;

            }
        }
        else if (vof[] == 1) {
            foreach_dimension()
                u.x[] = 0.;
        }
        foreach_dimension()
            utemp.x[] = u.x[];
    }
    boundary({u});

}


/*
event end_timestep (i++)
{
    // event ("advection_term1");
    
}
*/

face vector solidFlux[];
/*
event projection (i++)
{
    foreach() {
        if (vof[] > 0 && vof[] < 1) {
            coord n = interface_normal (point, vof);
            double alpha = plane_alpha (vof[], n);
            coord p[2];
            facets (n, alpha, p);
            
            coord cell = {x,y};

            if (p[0].x > p[1].x) {
                coord p_temp = p[1];
                p[1] = p[0];
                p[0] = p_temp;
            }

            

            solidFlux[].x = ;

            foreach_dimension()
                for (int i = 0; i < 2; i++)
                    p[i] = cell.x + p[i]*Delta;
        
            
        }
    }

    foreach() {
        if (id[] > 0) {
            
        }
    }
}
*/

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
            // coord marker = {markerCoord.x[], markerCoord.y[]};
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

#endif
