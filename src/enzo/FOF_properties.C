#include <stdio.h>
#include <stdlib.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_proto.h"

/************************************************************************/

void get_properties(FOFData D, FOF_particle_data *p, int len, float *pcm, 
		    float *pcmv, float *pmtot, float *pmstars)
{
  int i,k;
  double s[3], sv[3];
  double mtot, mstars;
  
  
  for (i = 0; i < 3; i++) {
    s[i] = 0;
    sv[i] = 0;
  }

  mtot = 0;
  mstars = 0;

  for (i = 0; i < len; i++) {
    for (k = 0; k < 3; k++) {
      s[k] += p[i].Mass * FOF_periodic(p[i].Pos[k] - p[0].Pos[k], D);
      sv[k] += p[i].Mass * p[i].Vel[k];
    }

    mtot += p[i].Mass;
    switch (p[i].Type) {
    case PARTICLE_TYPE_STAR:
    case PARTICLE_TYPE_SINGLE_STAR:
    case PARTICLE_TYPE_CLUSTER:
      mstars += p[i].Mass; 
      break;
    } // ENDSWITCH
  } // ENDFOR

  for(k=0; k<3; k++) {
    s[k] = s[k] / mtot;
    s[k] = FOF_periodic_wrap(s[k] + p[0].Pos[k], D);
    sv[k] = sv[k] / mtot;
  }

  for(k=0; k<3; k++) {
    pcm[k] = s[k];
    pcmv[k] = sv[k];
  }

  // Convert to solar masses, 0->1 domain, and km/s
  mtot *= 1e10;
  mstars *= 1e10;
  for (dim = 0; dim < 3; dim++) {
    pcm[dim] /= D.Boxsize;
    pcmv[dim] /= D.UnitVelocity_in_cm_per_s;
  }
  
  *pmtot   = mtot;
  *pmstars = mstars;

}
