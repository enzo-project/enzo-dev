#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_proto.h"

/************************************************************************/

void get_properties(FOFData D, FOF_particle_data *p, int len, float *pcm, 
		    float *pmtot, float *pmgas, float *pmstars, float *pmsfr, 
		    float *pmcold)
{
  int i,k;
  double s[3];
  double mtot, mgas, mstars, sfr, mcold;
  
  
  s[0] = s[1] = s[2] = 0;
  mtot = mgas = mstars = mcold = 0;
  sfr = 0;

  for (i = 0; i < len; i++) {
    for (k = 0; k < 3; k++)
      s[k] += p[i].Mass * FOF_periodic(p[i].Pos[k] - p[0].Pos[k], D);

    mtot += p[i].Mass;
    switch (p[i].Type) {
    case 0:
      mgas +=    (p[i].Mass-p[i].Mfs-p[i].Mclouds);
      mcold +=   p[i].Mclouds;
      mstars +=  p[i].Mfs;
      sfr +=  p[i].Sfr;
      break;
    case 4:
      mstars += p[i].Mass; break;
    } // ENDSWITCH
  } // ENDFOR

  for(k=0; k<3; k++) {
    s[k] = s[k]/mtot;
    s[k] = FOF_periodic_wrap(s[k] + p[0].Pos[k], D);
  }

  for(k=0; k<3; k++)
    pcm[k] = s[k];
  
  *pmtot   = mtot;
  *pmgas   = mgas;
  *pmstars = mstars;
  *pmsfr   = sfr;
  *pmcold  = mcold;

}
