#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>


#include "allvars.h"
#include "proto.h"




void get_properties(struct particle_data *p, int len, float *pcm, float *pmtot, float *pmgas, float *pmstars, float *pmsfr, float *pmcold)
{
  int i,k;
  double s[3];
  double mtot, mgas, mstars, sfr, mcold;
  

  
  s[0]= s[1]= s[2]= 0;
  mtot= mgas= mstars= mcold= 0;
  sfr= 0;

  for(i=0; i<len; i++)
    {
      for(k=0; k<3; k++)
	s[k]+= p[i].Mass * periodic(p[i].Pos[k] - p[0].Pos[k]);

      mtot+= p[i].Mass;
      switch(p[i].Type)
	{
	case 0:
	  mgas+=    (p[i].Mass-p[i].Mfs-p[i].Mclouds);
	  mcold+=   p[i].Mclouds;
	  mstars+=  p[i].Mfs;
	  sfr+=  p[i].Sfr;
	  break;
	case 4:
	  mstars+= p[i].Mass; break;
	}
    }

  for(k=0; k<3; k++)
    {
      s[k]= s[k]/mtot;
      s[k]= periodic_wrap(s[k] + p[0].Pos[k]); 
    }

  for(k=0; k<3; k++)
    pcm[k]= s[k];
  
  *pmtot=   mtot;
  *pmgas=   mgas;
  *pmstars= mstars;
  *pmsfr=   sfr;
  *pmcold=  mcold;
}
