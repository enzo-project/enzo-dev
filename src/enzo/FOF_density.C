#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "allvars.h"
#include "nrsrc/nrutil.h"
#include "proto.h"



void density(void)
{
  float  *r2list;
  int    *ngblist;
  int    i,k,ii,ind,signal;
  double h,h2,hinv3,hv_inv3,wk,u,r;
  
  set_sph_kernel();
  
  if(ThisTask==0)
    printf("Computing densities...\n");

  for(i=1,signal=0; i<=NumInGroup; i++)
    {
      if(ThisTask==0)
	{
	  if(i> (signal/100.0)*NumInGroup)
	    {
	      if((signal%10)==0)
		printf("%d",signal);
	      else
		printf(".",signal);
	      fflush(stdout);
	      signal++;
	    }
	}
      
      h2= ngb_treefind(P[i].Pos, DesDensityNgb, 0, &ngblist, &r2list); 
      
      h=sqrt(h2);

      hinv3 = 1.0/(h*h2);
	  
      Density[i]=0;
      
      for(k=0; k<DesDensityNgb; k++)
	{
	  r=sqrt(r2list[k]);
	  
	  ind=ngblist[k];
	  
	  if(r<h)
	    {
	      u = r/h;
	      ii = (int)(u*KERNEL_TABLE);
	      if(ii>=KERNEL_TABLE)  /* may occur with pgcc for certain optimization levels */
		wk= 0;
	      else
		wk =hinv3*( Kernel[ii]    + (Kernel[ii+1]-Kernel[ii])*(u-KernelRad[ii])*KERNEL_TABLE);
	      Density[i] += P[ind].Mass*wk; 
	    }      
	}
    }
  
  if(ThisTask==0)
    printf("\ndone.\n");
}




void set_sph_kernel(void)
{
  int i;
  FILE *fd;

  for(i=0;i<=KERNEL_TABLE;i++)
    KernelRad[i] = ((float)i)/KERNEL_TABLE;
      
  for(i=0;i<=KERNEL_TABLE;i++)
    {
      if(KernelRad[i]<=0.5)
	{
	  Kernel[i] = 8/PI *(1-6*KernelRad[i]*KernelRad[i]*(1-KernelRad[i]));
	  KernelDer[i] = 8/PI *( -12*KernelRad[i] + 18*KernelRad[i]*KernelRad[i]);
	}
      else
	{
	  Kernel[i] = 8/PI * 2*(1-KernelRad[i])*(1-KernelRad[i])*(1-KernelRad[i]);
	  KernelDer[i] = 8/PI *( -6*(1-KernelRad[i])*(1-KernelRad[i]));
	}
    }
}

