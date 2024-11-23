#include "fractions.h"
#include "ibm-gfm-utils.h"

scalar vof[];
face vector sf[];

#include "ibm-tree.h"

coord vc = {0,0};

extern face vector uf;

scalar fused[];
scalar id[];
scalar q[];         // mass source/sink term

vector normals[];
vector midPoints[];

vector utemp[];

#undef SEPS
#define SEPS 1e-30

#define face_condition(sf, vof)                 \
    (sf.x[i,j] > 0.5 &&                         \
     sf.y[i,j + (j < 0)] &&                     \
     sf.y[i-1,j + (j < 0)] &&                   \
     vof[i,j] && vof[i-1,j])

#define vof_avg(a,i,j,k)							\
  ((a[i,j,k]*(1.5 + vof[i,j,k]) + a[i-1,j,k]*(1.5 + vof[i-1,j,k]))/	\
   (vof[i,j,k] + vof[i-1,j,k] + 3.))

foreach_dimension()
static inline double ibm_face_gradient_x (Point point, scalar a, int i)
{
    int j = sign(sf.x[i,1] - sf.x[i,-1]);
    assert (vof[i] && vof[i-1]);

    if (face_condition (sf, vof))
        return ((1. + sf.x[i]) * (a[i] - a[i-1]) +
                (1. - sf.x[i]) * (a[i,j] - a[i-1,j])) / (2. * Delta);

    return (a[i] - a[i-1]) / Delta;
}

foreach_dimension()
static inline double ibm_face_value_x (Point point, scalar a, int i)
{
  int j = sign(sf.x[i,1] - sf.x[i,-1]);
  return face_condition (sf, vof) ?
    ((1. + sf.x[i])*vof_avg(a,i,0,0) + (1. - sf.x[i])*vof_avg(a,i,j,0))/2. :
    vof_avg(a,i,0,0);
}

#undef face_gradient_x
#define face_gradient_x(a,i)                    \
    (sf.x[i] < 1. && sf.x[i] > 0. ?             \
     ibm_face_gradient_x (point, a, i) :        \
     (a[i] - a[i-1]) / Delta)

#undef face_gradient_y
#define face_gradient_y(a,i)                    \
    (sf.y[0,i] < 1. && sf.y[0,i] > 0. ?         \
     ibm_face_gradient_x (point, a, i) :        \
     (a[0,i] - a[0,i-1]) / Delta)


#undef face_value
#define face_value(a,i)							\
  (sf.x[i] < 1. && sf.x[i] > 0. ?				\
   ibm_face_value_x (point, a, i) :					\
   vof_avg(a,i,0,0))


#undef center_gradient
#define center_gradient(a) (sf.x[] && sf.x[1] ? (a[1] - a[-1])/(2.*Delta) : \
			    sf.x[1] ? (a[1] - a[])/Delta :		    \
			    sf.x[]  ? (a[] - a[-1])/Delta : 0.)


static inline double bilinear_ibm (Point point, scalar s)
{
  if (!coarse(vof) || !coarse(vof,child.x))
    return coarse(s);
  #if dimension >= 2
  if (!coarse(vof,0,child.y) || !coarse(vof,child.x,child.y))
    return coarse(s);
  #endif
  return bilinear (point, s);
}

#define bilinear(point, s) bilinear_ibm(point, s)


event metric (i = 0)
{
  if (is_constant (fm.x)) {
    foreach_dimension()
      assert (constant (fm.x) == 1.);
    fm = sf;
  }
  foreach_face() {
    sf.x[] = 1.;
  }
  if (is_constant (cm)) {
    assert (constant (cm) == 1.);
    cm = vof;
  }
  foreach() {
    vof[] = 1.;
  }

    vof.refine = ibm_fraction_refine;
    vof.prolongation = fraction_refine;
    //vof.refine = vof.prolongation = fraction_refine;
    //vof.dirty = true;
    foreach_dimension() {
        sf.x.prolongation = ibm_face_fraction_refine_x;
        // sf.x.refine = sf.x.prolongation = fraction_refine;
        // sf.x.dirty = true;
    }
  restriction ({vof, sf});
}

/*
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

/*
event advection_term (i++)
{
    for (scalar s in {p, u, g, pf})
        foreach()
            if (vof[] >= 1)
	            s[] = 0.;
}
*/


/*
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
        }
    }

    // 4. Calculate image point values

    coord imageVelos[numCells];

    for (int i = 0; i < numCells; i++) {
        foreach_point(imagePoints[i].x, imagePoints[i].y) {
            imageVelos[i] = image_velocity (point, u, imagePoints[i]); 
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

event adapt (i++,last) {
  // fractions_cleanup (cs, fs);
  foreach_face()
    if (uf.x[] && !sf.x[])
      uf.x[] = 0.;
  event ("properties");
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
