#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>


#include "allvars.h"
#include "proto.h"




void get_properties(struct particle_data *p, int len, float *pcm, float *pmtot, 
		    float *pmgas, float *pmstars, float *pmsfr, float *pmcold,
		    int subgroup, float *pcmv, float *pmvir, float *prvir,
		    float *pL, float *pvrms, float *pspin)
{
  int i,k,dim, irvir;
  double s[3], sv[3], L[3], delx[3], delv[3], del, vrms, spin, mvir, rvir;
  double mtot, mgas, mstars, sfr, mcold, menc, rho, factor, rho200, r3;
  double dwidth_inv;
  float *radius;
  int *pindex;

  for (i = 0; i < 3; i++) {
    s[i] = 0.0;
    sv[i] = 0.0;
    L[i] = 0.0;
  }
    
  mtot= mgas= mstars= mcold= 0;
  sfr= 0;

  for(i=0; i<len; i++)
    {
      for(k=0; k<3; k++) {
	s[k]+= p[i].Mass * periodic(p[i].Pos[k] - p[0].Pos[k]);
	sv[k] += p[i].Mass * p[i].Vel[k];
      }

      mtot+= p[i].Mass;
      switch(p[i].Type)
	{
	case 0:
	  mgas+=    (p[i].Mass-p[i].Mfs-p[i].Mclouds);
	  mcold+=   p[i].Mclouds;
	  mstars+=  p[i].Mfs;
	  sfr+=  p[i].Sfr;
	  break;
	case 2:
	case 5:
	case 7:
	  mstars+= p[i].Mass; 
	  break;
	}
    }

  for(k=0; k<3; k++)
    {
      s[k]= s[k]/mtot;
      s[k]= periodic_wrap(s[k] + p[0].Pos[k]); 
      sv[k] = sv[k] / mtot;
    }

  for(k=0; k<3; k++) {
    pcm[k]= s[k];
    pcmv[k] = sv[k];
  }

  /* For groups, find the virial radius (r200) and calculate virial mass */

  // Sort by radius and search for an enclosed density of 200 times
  // the critical density.  Search outside-in.

  pindex = (int*) malloc(len*sizeof(int));

  if (subgroup == 0) {
    radius = (float*) malloc(len*sizeof(float));
    for (i = 0; i < len; i++) {
      radius[i] = 0;
      for (dim = 0; dim < 3; dim++) {
	del = periodic(p[i].Pos[dim] - pcm[dim]);
	radius[i] += del*del;
      }
      radius[i] = sqrt(radius[i]);
    }
    indexx(len, radius-1, pindex-1);

    // Convert one-based (from indexx) to zero-based.
    for (i = 0; i < len; i++)
      pindex[i]--;

    // Convert to Msun from 1e10 Msun and pre-compute the (4PI/3)
    // factor.  Rho will be in units of Msun / kpc^3, as is rho_crit.
    // Radius is comoving.
    factor = 1e10 / (4*M_PI/3.0);
    rho200 = 200 * RhoCritical0;

    menc = mtot;
    for (i = len-1; i >= 0; i--) {
      menc -= p[pindex[i]].Mass;
      r3 = radius[pindex[i]] * radius[pindex[i]] * radius[pindex[i]];
      rho = factor * menc / max(r3, 1e-20);
      if (rho > rho200) 
	break;
    }

    irvir = max(i, 0);
    rvir = radius[pindex[irvir]] * Time;  // comoving -> proper
    mvir = menc + p[pindex[irvir]].Mass;  // Add the last particle removed

    free(radius);

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
    L[dim] /= mvir * 1e3 / Time;
  vrms /= mvir;
  vrms = sqrt(vrms);

  // Convert to solar masses, 0->1 domain
  mvir *= 1e10;
  mtot *= 1e10;
  mstars *= 1e10;
  dwidth_inv = 1.0/(rightEdge[0] - leftEdge[0]);
  for (dim = 0; dim < 3; dim++) {
    pcm[dim] /= BoxSize * dwidth_inv;
    pcm[dim] += leftEdge[dim];
  }

  // Spin parameter
  float ang_mom = 0;
  float SpinUnits = (CM_PER_MPC * 1.0e5 * 1.0e5) / (GRAVITY * SOLAR_MASS);
  for (dim = 0; dim < 3; dim++)
    ang_mom += L[dim] * L[dim];
  ang_mom = sqrt(ang_mom);
  spin = SpinUnits * ang_mom * vrms / mvir;

  free(pindex);
  
  *pmtot=   mtot;
  *pmgas=   mgas;
  *pmstars= mstars;
  *pmsfr=   sfr;
  *pmcold=  mcold;
  *pvrms = vrms;
  *pspin = spin;
  *pmvir = mvir;
  *prvir = rvir;
  for (dim = 0; dim < 3; dim++)
    pL[dim] = L[dim];

}
