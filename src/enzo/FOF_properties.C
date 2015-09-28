#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_proto.h"

/************************************************************************/

void get_properties(FOFData D, FOF_particle_data *p, int len, bool subgroup,
		    float *pcm, 
		    float *pcmv, float *pmtot, float *pmstars, float *pmvir,
		    float *prvir, float *pL, float *pvrms, float *pspin)
{
  int i,k,dim, irvir, len4;
  double s[3], sv[3], L[3], delx[3], delv[3], vrms, spin, mvir, rvir, del;
  double mtot, mtot0, mstars, menc, rho, factor, rho178, r3, e2;
  float *radius;
  int *pindex;
  
  vrms = 0;
  for (i = 0; i < 3; i++) {
    s[i] = 0;
    sv[i] = 0;
    L[i] = 0;
  }

  mtot = 0;
  mstars = 0;

  for (i = 0; i < len; i++) {
    for (k = 0; k < 3; k++) {
      s[k] += p[i].Mass * FOF_periodic(p[i].Pos[k] - p[0].Pos[k], D.BoxSize);
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
    s[k] = FOF_periodic_wrap(s[k] + p[0].Pos[k], D.BoxSize);
    sv[k] = sv[k] / mtot;
  }

  for(k=0; k<3; k++) {
    pcm[k] = s[k];
    pcmv[k] = sv[k];
  }

  mtot0 = mtot;
  pindex = new int[len];
  radius = new float[len];

  /* Search for iterative center of mass, decreasing sphere radius by
     5% in each step */

  int ninside = len;
  float rmax;
  bool first = true;
  while (ninside > 2) {
    for (i = 0; i < len; i++) {
      radius[i] = 0;
      for (dim = 0; dim < 3; dim++) {
	del = FOF_periodic(p[i].Pos[dim] - pcm[dim], D.BoxSize);
	radius[i] += del*del;
      }
      radius[i] = sqrt(radius[i]);
    }
    indexx(len, radius-1, pindex-1);
    for (i = 0; i < len; i++)
      pindex[i]--;
    // Search for particle with r ~ 0.95*(previous search radius)
    ninside = len;
    if (!first) {
      rmax *= 0.95;
      while (radius[pindex[ninside-1]] > rmax)
	ninside--;
    } else {
      rmax = radius[pindex[ninside-1]];
    }
    // Re-calculate center of mass
    mtot = 0.0;
    for (dim = 0; dim < 3; dim++)
      pcm[dim] = 0;
    for (i = 0; i < ninside; i++) {
      for (dim = 0; dim < 3; dim++) {
	pcm[dim] += p[pindex[i]].Mass * 
	  FOF_periodic(p[pindex[i]].Pos[dim] - p[0].Pos[dim], D.BoxSize);
      }
      mtot += p[pindex[i]].Mass;
    }
    for (dim = 0; dim < 3; dim++) {
      pcm[dim] /= mtot;
      pcm[dim] = FOF_periodic_wrap(pcm[dim] + p[0].Pos[dim], D.BoxSize);
    }
    first = false;
  } // ENDWHILE (ninside > 2)

  /* For groups, find the virial radius (r200) and calculate virial mass */

  // Sort by radius and search for an enclosed density of 18*pi^2 times
  // the critical density.  Search outside-in.

  if (!subgroup) {
    for (i = 0; i < len; i++) {
      radius[i] = 0;
      for (dim = 0; dim < 3; dim++) {
	del = FOF_periodic(p[i].Pos[dim] - pcm[dim], D.BoxSize);
	radius[i] += del*del;
      }
      radius[i] = sqrt(radius[i]);
    }
    indexx(len, radius-1, pindex-1);

    // Convert one-based (from indexx) to zero-based.
    for (i = 0; i < len; i++)
      pindex[i]--;

    // Convert to Msun from 1e10 Msun and pre-compute the (4PI/3)
    // factor and comoving factors.  Rho will be in units of Msun /
    // kpc^3, as is rho_crit.
    e2 = D.OmegaLambda + D.Omega * pow(D.Time, -3.0);
    // 1/a^3 factor converts comoving radius to proper.
    factor = 1e10 / (4*M_PI/3.0) * pow(D.Time, -3.0);
    rho178 = 178.0 * D.RhoCritical0 * e2;
    len4 = len/4;

    mtot = mtot0;
    menc = mtot;
    for (i = len-1; i >= 0; i--) {
      menc -= p[pindex[i]].Mass;
      r3 = radius[pindex[i]] * radius[pindex[i]] * radius[pindex[i]];
      rho = factor * menc / max(r3, tiny_number);
      if (rho > rho178) 
	break;
    }

    irvir = max(i, 0);
    rvir = radius[pindex[irvir]] * D.Time;  // comoving -> proper
    mvir = menc + p[pindex[irvir]].Mass;  // Add the last particle removed

  } // ENDIF !subgroup

  else {

    // Skip finding the overdense sphere for subgroups

    irvir = len-1;
    rvir = 0;
    mvir = mtot;
    for (i = 0; i < len; i++)
      pindex[i] = i;

  } // ENDELSE !subgroup

  /* Calculate other quantities :: RMS velocity, angular momentum (Mpc * km/s),
     spin parameter */

  for (i = 0; i <= irvir; i++) {
    k = pindex[i];
    for (dim = 0; dim < 3; dim++) {
      delx[dim] = p[k].Pos[dim] - pcm[dim];
      delv[dim] = p[k].Vel[dim] - pcmv[dim];
      vrms += p[k].Mass * delv[dim] * delv[dim];
    }

    L[0] += p[k].Mass * ( delv[1]*delx[2] - delv[2]*delx[1]);
    L[1] += p[k].Mass * (-delv[0]*delx[2] + delv[2]*delx[0]);
    L[2] += p[k].Mass * ( delv[0]*delx[1] - delv[1]*delx[0]);

  } // ENDFOR particles

  /* Divide out weight (m_vir) */
  // 1e3 to convert from kpc to Mpc, comoving -> proper
  for (dim = 0; dim < 3; dim++)
    L[dim] /= mvir * 1e3 / D.Time;
  vrms /= mvir;
  vrms = sqrt(vrms);

  // Convert to solar masses, 0->1 domain
  mvir *= 1e10;
  mtot *= 1e10;
  mstars *= 1e10;
  for (dim = 0; dim < 3; dim++)
    pcm[dim] /= D.BoxSize;

  // Spin parameter
  float ang_mom = 0;
  float SpinUnits = (CM_PER_MPC * 1.0e5 * 1.0e5) / (GRAVITY * SOLAR_MASS);
  for (dim = 0; dim < 3; dim++)
    ang_mom += L[dim] * L[dim];
  ang_mom = sqrt(ang_mom);
  spin = SpinUnits * ang_mom * vrms / mvir;

  delete [] radius;
  delete [] pindex;

  *pmtot   = mtot;
  *pmstars = mstars;
  for (dim = 0; dim < 3; dim++)
    pL[dim] = L[dim];
  *pvrms = vrms;
  *pspin = spin;
  *pmvir = mvir;
  *prvir = rvir;

}
