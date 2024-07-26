# -*- coding: utf-8 -*-
# #include "embed.h"
# #include "navier-stokes/centered.h"
# #include "view.h"
# 
# #define L0 20
# #define D 0.5
# #define LEVEL 11
# 
# int Re;
# double U0 =  1.0; // inlet velocity
# double t_end = 20 [0,1];
# double tf_start = 20 [0,1];
# double xi = 4.8266;
# double yi = 3.0835;
# coord ci = {5, 10}; // initial coordinates of cylinder
# 
# face vector muv[];
# 
# u.n[left] = dirichlet ((U0));
# u.t[left] = dirichlet (0);
# p[left]   = neumann (0);
# pf[left]  = neumann (0);
# 
# u.n[right] = neumann (0);
# u.t[right] = neumann (0);
# p[right]   = dirichlet (0);
# pf[right]  = dirichlet (0);
# 
# u.n[top] = neumann (0);
# u.t[top] = dirichlet (U0);
# p[top] = neumann (0);
# pf[top] = neumann (0);
# 
# u.n[bottom] = neumann (0);
# u.t[bottom] = dirichlet (U0);
# p[bottom] = neumann (0);
# pf[bottom] = neumann (0);
# 
# double SDF (double x, double y) {
#    return - sq(x - ci.x) - sq(y - ci.y) + sq(D/2);
# }
# 
# int main() {
#   size(L0);
#   init_grid (1 << (LEVEL - 3));
#   mu = muv;
#   TOLERANCE = 1.e-6 [*]; 
# 
#   Re = 1;
#   run();
# 
#   Re = 2;
#   run();
# 
#   Re = 5;
#   run();
#   
#   Re = 10;
#   run();
#   
#   Re = 20;
#   run();
# 
#   Re = 40;
#   run();
# 
#   Re = 50;
#   run();
# }
# 
# 
# event init (t = 0) {
#   // mask(y > 6 ? top: y < -6 ? bottom : none);
#   refine (level <= LEVEL*(1. - sqrt(fabs(sq(x-ci.x) + sq(y-ci.y) - sq(D/2.)))/2.));
#   solid (cs, fs, sq(x - ci.x) + sq(y - ci.y) - sq(D/2));
#   foreach()
#     u.x[] = cs[] ? U0 : 0.;
#   
#   u.n[embed] = dirichlet(0);
#   u.t[embed] = dirichlet(0);  
# 
# }
# 
# 
# event properties (i++) {
#   foreach_face()
#     muv.x[] = fm.x[]*(U0)*(D)/(Re);
#  // boundary ((scalar *) {muv});
# }
# 
# 
# event logfile (i++){
#   coord Fp, Fmu;
#   embed_force (p, u, mu, &Fp, &Fmu);
#   double CD = (Fp.x + Fmu.x)/(0.5*sq(U0)*(D));
#   double CL = (Fp.y + Fmu.y)/(0.5*sq(U0)*(D));
# 
#  /* 
#   double E = 0;
#   double E_p = 0;
#   boundary ({u.x, u.y});
#   scalar omega[];
#   vorticity (u , omega);
#   foreach(){
#     double vort = omega[];
#     double area = dv();
#     if (cs[] < 1. && cs[] > 0){
#       coord b, n;
#       area *= embed_geometry (point, &b, &n);
#       vort = embed_vorticity (point, u, b, n);
#     }
#     E += area*sq(vort);
#     E_p += sq(vort);
# 
#   }
# */
#   fprintf (stderr, "%d %g %d %d %d %d %d %g %g\n",
# 	   i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL);
# 	   
# }
# 
# 
# event snapshot (t = t_end) {
#   scalar omega[];
#   vorticity (u, omega);
# 
#   char name[80];
#   sprintf (name, "vort-%d", Re);
#   FILE * fp1 = fopen (name, "w");
#   view (fov = 2, tx = -0.375, ty = -0.20,
# 	width = 800, height = 400); 
#   isoline ("omega", n = 15, min = -3, max = 3);
#   draw_vof ("cs", "fs", filled = 1, lw = 5);
#   save (fp = fp1);
# }
# 
# event adapt (i++) {
#   adapt_wavelet ({cs,u}, (double[]){1.e-2,3e-3,3e-3},
# 		 maxlevel = LEVEL, minlevel = 2);
# }
# 
# event profile (t = t_end) {
#   int k = 0;
#   double delta = L0/(pow(2,LEVEL));
#   char name[80];
# 
#   sprintf (name, "vprofx1-%d", Re); // x = 4.8125
#   FILE * fv = fopen(name, "w");
#   for(double i = 0; i <= L0; i += delta) {
#     foreach_point (4.8125, i) {
#       if (cs[] > 0 && cs[] < 1)
#         k = 2.;
#       else if (cs[] == 1)
# 	k = 1.;
#       else
# 	k = 0.;
#       fprintf (fv, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
#     }
#   }
#   fflush (fv);
#   fclose (fv);
# 
#   sprintf (name, "vprofx2-%d", Re); // x = 5
#   FILE * fv1 = fopen(name, "w");
#   for(double i = 0; i <= L0; i += delta) {
#     foreach_point (5, i) {
#       if (cs[] > 0 && cs[] < 1)
#         k = 2.;
#       else if (cs[] == 1)
# 	k = 1.;
#       else
# 	k = 0.;
#       fprintf (fv1, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
#     }
#   }
#   fflush (fv1);
#   fclose (fv1);
# 
#   sprintf (name, "vprofx3-%d", Re); // x = 10
#   FILE * fv2 = fopen(name, "w");
#   for(double i = 0; i <= L0; i += delta) {
#     foreach_point (10, i) {
#       if (cs[] > 0 && cs[] < 1)
#         k = 2.;
#       else if (cs[] == 1)
# 	k = 1.;
#       else
# 	k = 0.;
#       fprintf (fv2, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
#     }
#   }
#   fflush (fv2);
#   fclose (fv2);
# 
#   sprintf (name, "vprofy1-%d", Re); // y = 3
#   FILE * fv3 = fopen(name, "w");
#   for(double i = 0; i <= L0; i += delta) {
#     foreach_point (i, 3) {
#       if (cs[] > 0 && cs[] < 1)
#         k = 2.;
#       else if (cs[] == 1)
# 	k = 1.;
#       else
# 	k = 0.;
#       fprintf (fv3, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
#     }
#   }
#   fflush (fv3);
#   fclose (fv3);
# 
#   sprintf (name, "vprofy2-%d", Re); // y = 3.1875
#   FILE * fv4 = fopen(name, "w");
#   for(double i = 0; i <= L0; i += delta) {
#     foreach_point (i, 3.1875) {
#       if (cs[] > 0 && cs[] < 1)
#         k = 2.;
#       else if (cs[] == 1)
# 	k = 1.;
#       else
# 	k = 0.;
#       fprintf (fv4, "%d %g %g %g %g\n", k, x, y, u.x[], u.y[]);
#     }
#   }
#   fflush (fv4);
#   fclose (fv4);
#  
# }
# 
# int count = 0;
# event frequency (i++) {
#   if (t >= tf_start) {
#     char name[80];
#     sprintf (name, "freq-7.5.dat");
#     FILE * fp = fopen (name, "a");
#     foreach_point(7.500, ci.y) {
#       fprintf (fp, "%d %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], p[]);
#     }
#     fclose (fp);
# 
#     sprintf (name, "freq-10.dat");
#     FILE * fp1 = fopen (name, "a");
#     foreach_point(10.000, ci.y) {
#       fprintf (fp1, "%d %g %g %g %g %g %g\n", count, x, y, t, u.x[], u.y[], p[]);
#     }
#     fclose (fp1);
#     count++;
#   }
# }
# 
# event movie (t += 0.01; t <= t_end) {
# }
# 
# event stop (t = t_end) {
#   static FILE * fp = fopen("perf", "w");
#   timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
#   fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
#   fflush (fp);
#   return 1;
# }
