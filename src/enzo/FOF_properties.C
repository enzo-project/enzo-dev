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
		    float *pmtot, float *pmstars)
{
  int i,k;
  double s[3];
  double mtot, mstars;
  
  
  s[0] = s[1] = s[2] = 0;
  mtot = 0;

  for (i = 0; i < len; i++) {
    for (k = 0; k < 3; k++)
      s[k] += p[i].Mass * FOF_periodic(p[i].Pos[k] - p[0].Pos[k], D);

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
    s[k] = s[k]/mtot;
    s[k] = FOF_periodic_wrap(s[k] + p[0].Pos[k], D);
  }

  for(k=0; k<3; k++)
    pcm[k] = s[k];
  
  *pmtot   = mtot;
  *pmstars = mstars;

}
