#ifndef BASILISK_HEADER_19
#define BASILISK_HEADER_19
#line 1 "./../immersed-old.h"
#include "fractions.h"

extern coord vc; // solid velocity
extern scalar vof;
extern face vector sf;
extern int maxlevel;

vector fc[]; // body force term (acceleration units)

event end_timestep (i++) {
  correction (-dt);
  foreach()
    foreach_dimension() {
      fc.x[] = vof[]*(vc.x - u.x[])/dt;
      g.x[] += fc.x[];
    }
  correction (dt);
}

#endif
