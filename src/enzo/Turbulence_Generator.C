/***********************************************************************
/
/  GENERATE TURBULENT VELOCITY FIELD
/
/  written by: Peng Wang
/  date:       April, 2007
/  modified1: Tom Abel (September 2009) moved normalization outside to 
/                            allow multiple grids and parallel setups
/
/
************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "macros_and_parameters.h"
extern int NumberOfGhostZones;


double Gaussian(double cs); 

 void Turbulence_Generator(float **vel, int dim0, int dim1, int dim2, int ind,  
			   float kmin, float kmax, float dk,
			   FLOAT **LeftEdge, FLOAT **CellWidth, int seed)
   /* 
      vel[3][ActiveSize]
      size: the grid dimension
      ind: index of the velocity power spectrum v_k^2 ~ k^{-n}, 
           n = 11/3 for Kolmogorov turbulence,
	   n = 4 for Larson relation
      kmin: the lower wave number cutoff
      kmax: the upper wave number cutoff 
   */
 {
   int igrid, i, j, k;
   igrid = 0;
   for (k = 0; k < dim2; k++) {
     for (j = 0; j < dim1; j++) {
       for (i = 0; i < dim0; i++, igrid++) {
	 vel[0][igrid] = 0.0;
	 vel[1][igrid] = 0.0;
	 vel[2][igrid] = 0.0;
       }
     }
   }

   printf("Turbulence_Generator: seed=%"ISYM", kmin=%"GSYM", kmax=%"GSYM", dk = %"GSYM"\n", 
	  seed, kmin, kmax, dk);
   srand(seed);

   double phix, phiy, phiz, Ax, Ay, Az, AA, Ak0;
   double pi = 4.0*atan(1.0);
   double k_wave;
   float kx, ky, kz, k2;
   for (kz = 0; kz <= kmax; kz+=dk) {
     printf("kz=%"GSYM"\n", kz);
     for (ky = 0; ky <= kmax; ky+=dk) {
       for (kx = 0; kx <= kmax; kx+=dk) {

	 k2 = kx*kx + ky*ky + kz*kz;

	 if (k2 < kmin*kmin || k2 > kmax*kmax) 
	   continue;

	 /* Generate complex vector A_k = Ak0*(Ax*exp(i*phix), Ay*exp(i*phiy), Az*exp(iphiz))  */

	 Ax = rand();
	 Ax /= RAND_MAX;
	 Ay = rand();
	 Ay /= RAND_MAX;
	 Az = rand();
	 Az /= RAND_MAX;
	 AA = sqrt(Ax*Ax+Ay*Ay+Az*Az);
	 Ax /= AA;
	 Ay /= AA;
	 Az /= AA;

	 k_wave = sqrt(k2)*2.0*pi;
	 Ak0 = Gaussian(pow(k_wave, -ind*0.5-1));

	 phix = rand();
	 phix = phix / RAND_MAX * 2.0 * pi;
	 phiy = rand();
	 phiy = phiy / RAND_MAX * 2.0 * pi;
	 phiz = rand();
	 phiz = phiz / RAND_MAX * 2.0 * pi;

	 /* no helicity ( in the case of magnetic fields) */ 
	 phix = phiy = phiz = 0.;

	 /* Fourier transform A_k to v(x) 
	  v(x) = Sum_k{Re(ik x A_k)}*/

	 igrid = 0;
	 for (k = 0; k < dim2; k++) {
	   for (j = 0; j < dim1; j++) {
	     for (i = 0; i < dim0; i++, igrid++) {
	       int ii = i + NumberOfGhostZones;
	       int jj = j + NumberOfGhostZones;
	       int kk = k + NumberOfGhostZones;
	       double x = LeftEdge[0][ii]+(0.5)*CellWidth[0][ii];
	       double y = LeftEdge[1][jj]+(0.5)*CellWidth[1][jj];
	       double z = LeftEdge[2][kk]+(0.5)*CellWidth[2][kk];
	       double kdotx = 2.0*pi*(kx*x + ky*y + kz*z);
	       vel[0][igrid] += Ak0*(-ky*Az*sin(kdotx+phiz) + kz*Ay*sin(kdotx+phiy));
	       vel[1][igrid] += Ak0*(-kz*Ax*sin(kdotx+phix) + kx*Az*sin(kdotx+phiz));
	       vel[2][igrid] += Ak0*(-kx*Ay*sin(kdotx+phiy) + ky*Ax*sin(kdotx+phix));
	     }
	   }
	 }
	   //	   printf("%"GSYM" \n", CellWidth[0][0]);

       }
     }
   }

   return ;

}

