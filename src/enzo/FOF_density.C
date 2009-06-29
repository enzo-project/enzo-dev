#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_nrutil.h"
#include "FOF_proto.h"


void density(FOFData &A)
{
  float  *r2list;
  int    *ngblist;
  int    i,k,ii,ind,signal;
  double h,h2,hinv3,hv_inv3,wk,u,r;
  
  set_sph_kernel(A);
  
//  if (debug)
//    printf("Computing densities...\n");

  for (i = 1, signal = 0; i <= A.NumInGroup; i++) {
//    if (debug) {
//      if (i > (signal/100.0)*A.NumInGroup) {
//	if ((signal%10) == 0)
//	  printf("%"ISYM,signal);
//	else
//	  printf(".",signal);
//	fflush(stdout);
//	signal++;
//      }
//    } // ENDIF debug

    h2 = ngb_treefind(A.P, A.P[i].Pos, A.DesDensityNgb, 0, &ngblist, &r2list); 
    
    h = sqrt(h2);

    hinv3 = 1.0/(h*h2);
	  
    A.Density[i] = 0;
      
    for (k = 0; k < A.DesDensityNgb; k++) {
      r = sqrt(r2list[k]);
      ind = ngblist[k];
      if (r < h) {
	u = r/h;
	ii = (int)(u*KERNEL_TABLE);
	if(ii >= KERNEL_TABLE)  /* may occur with pgcc for certain optimization levels */
	  wk= 0;
	else
	  wk =hinv3*( A.Kernel[ii] + (A.Kernel[ii+1]-A.Kernel[ii])*
		      (u-A.KernelRad[ii])*KERNEL_TABLE);
	A.Density[i] += A.P[ind].Mass*wk; 
	    }      
	}
    }
  
//  if (debug)
//    printf("\ndone.\n");
}




void set_sph_kernel(FOFData &A)
{
  int i;
  FILE *fd;

  for (i = 0; i <= KERNEL_TABLE; i++)
    A.KernelRad[i] = ((float) i) / KERNEL_TABLE;
  
  for (i = 0;i <= KERNEL_TABLE; i++) {
    if (A.KernelRad[i] <= 0.5) {
      A.Kernel[i] = 8/PI * (1-6*A.KernelRad[i]*A.KernelRad[i]* 
			    (1-A.KernelRad[i]));
      A.KernelDer[i] = 8/PI *( -12*A.KernelRad[i] + 
			       18*A.KernelRad[i]*A.KernelRad[i]);
    }
    else {
      A.Kernel[i] = 8/PI * 2*(1-A.KernelRad[i])*(1-A.KernelRad[i])*
	(1-A.KernelRad[i]);
      A.KernelDer[i] = 8/PI *( -6*(1-A.KernelRad[i])*(1-A.KernelRad[i]));
    }
  }
}

