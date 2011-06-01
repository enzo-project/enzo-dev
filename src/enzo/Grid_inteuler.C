//=======================================================================
/////////////////////////  SUBROUTINE INTEULER  \\\\\\\\\\\\\\\\\\\\\\\\\
//
//      subroutine inteuler(
//     &            dslice, pslice, gravity, grslice, geslice,
//     &            uslice, vslice, wslice, dxi, flatten,
//     &            idim, jdim, i1, i2, j1, j2, idual, eta1, eta2,
//     &            isteep, iflatten, dt, gamma, ipresfree,
//     &            dls, drs, pls, prs, gels, gers,
//     &            uls, urs, vls, vrs, wls, wrs,
//     &            ncolor, colslice, colls, colrs
//     &                    )
//
//  COMPUTES LEFT AND RIGHT EULERIAN INTERFACE VALUES FOR RIEMANN SOLVER
//
//  written by: Jim Stone
//  date:       January,1991
//  modified1:  June, 1993 by Greg Bryan (changed to eulerian)
//  modified2:  July, 1994 by Greg Bryan (changed to slicewise)
//  modified3:  April, 1995 by GB (added gravity)
//  modified4:  January, 2010 by JHW (translated to C++)
//
//  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
//    and right interface values to be fed into Riemann solver during a
//    one dimensional sweeps.  This version computes the Eulerian corrections
//    to the left and right states described in section three of Colella &
//    Woodward (1984), JCP.  The routine works on one two dimensional
//    slice at a time.
//
//  INPUT:
//    dslice - extracted 2d slice of the density, d
//    dt     - timestep in problem time
//    dxi    - distance between Eulerian zone edges in sweep direction
//    eta1   - (dual) selection parameter for gas energy (typically ~0.001)
//    flatten - ammount of flattening (calculated in calcdiss)
//    gamma  - ideal gas law constant
//    gravity - gravity flag (0 = off)
//    grslice - acceleration in this direction in this slice
//    i1,i2  - starting and ending addresses for dimension 1
//    idim   - declared leading dimension of slices
//    idual  - dual energy formalism flag (0 = off)
//    iflatten - integer flag for flattener (eq. A1, A2) (0 = off)
//    isteep - integer flag for steepener (eq. 1.14,1.17,3.2) (0 = off)
//    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
//    j1,j2  - starting and ending addresses for dimension 2
//    jdim   - declared second dimension of slices
//    pslice - extracted 2d slice of the pressure, p
//    uslice - extracted 2d slice of the 1-velocity, u
//    vslice - extracted 2d slice of the 2-velocity, v
//    wslice - extracted 2d slice of the 3-velocity, w
//    
//  OUTPUT:
//    dl,rs  - density at left and right edges of each cell
//    pl,rs  - pressure at left and right edges of each cell
//    ul,rs  - 1-velocity at left and right edges of each cell
//    vl,rs  - 2-velocity at left and right edges of each cell
//    wl,rs  - 3-velocity at left and right edges of each cell
//
//  EXTERNALS:
//
//  LOCALS:
//
//  PARAMETERS:
//    ft     - a constant used in eq. 1.124 (=2*2/3)
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

void intvarC(float *qslice, int in, int is, int ie, int j, int isteep,
	     float *steepen, int iflatten, float *flatten, float *c1,
	     float *c2, float *c3, float *c4, float *c5, float *c6,
	     float *char1, float *char2, float *c0, float *dq, float *ql,
	     float *qr, float *q6, float *qla, float *qra, float *ql0,
	     float *qr0);

int grid::inteuler(int idim,
		   float *dslice, float *pslice, int gravity, float *grslice,
		   float *geslice, float *uslice, float *vslice, float *wslice, 
		   float *dxi, float *flatten, 
		   float *dls, float *drs, float *pls, float *prs, float *gels,
		   float *gers, float *uls, float *urs, float *vls, float *vrs,
		   float *wls, float *wrs, int ncolors, float *colslice,
		   float *colls, float *colrs)
{

  /* Local variables */

  const int NN = MAX_ANY_SINGLE_DIRECTION;
  const int NC = GridDimension[idim]*MAX_COLOR;
  const float ft = 4.0/3.0;

  int i, j, ic, index, cindex, scindex;
  float steepen[NN], tmp1[NN], tmp2[NN], tmp3[NN], tmp4[NN];
  float qa,qb,qc,qd,qe,s1,s2,f1; 
  float c1[NN], c2[NN], c3[NN], c4[NN], c5[NN], c6[NN], dxinv[NN],
    dp[NN],  pl[NN],  pr[NN],  p6[NN],  du[NN],  ul[NN],  ur[NN], u6[NN],
    dla[NN], dra[NN], pla[NN], pra[NN], ula[NN], ura[NN], vla[NN],
    vra[NN], wla[NN], wra[NN], plm[NN], prm[NN], ulm[NN], urm[NN],
    dl0[NN], dr0[NN], pl0[NN], pr0[NN], plp[NN], prp[NN], ulp[NN],
    urp[NN], ul0[NN], ur0[NN], vl0[NN], vr0[NN], wl0[NN], wr0[NN],
    cs[NN],  d2d[NN], dxb[NN], cm[NN],  c0[NN],  cp[NN],  char1[NN],
    char2[NN], betalp[NN], betalm[NN], betal0[NN], cla[NN], betarp[NN], 
    betarm[NN], betar0[NN], cra[NN], gela[NN], gera[NN], gel0[NN], ger0[NN],
    clainv[NN], crainv[NN], dx2i[NN], dxdt2[NN];
  float colla[NC], colra[NC], coll0[NC], colr0[NC];
  //float colla[MAX_COLOR][NN], colra[MAX_COLOR][NN], coll0[MAX_COLOR][NN],
  //  colr0[MAX_COLOR][NN];

  /* Shorthand for data bounds and dimensions */

  int jdim = (idim+1) % 3;
  int in, jn, is, ie, js, je;

  in = GridDimension[idim];
  jn = GridDimension[jdim];
  is = GridStartIndex[idim];
  ie = GridEndIndex[idim];
  js = 0;
  je = GridDimension[jdim]-1;

//
// Compute coefficients used in interpolation formulae (from eq. 1.6)
//

  for (i = is-2; i <= ie+2; i++) {
    qa = dxi[i] / (dxi[i-1] + dxi[i] + dxi[i+1]);
    c1[i] = qa*(2.0*dxi[i-1] + dxi[i])/(dxi[i+1] + dxi[i]);
    c2[i] = qa*(2.0*dxi[i+1] + dxi[i])/(dxi[i-1] + dxi[i]);
    dxinv[i] = 1.0f/dxi[i];
    dxdt2[i] = 0.5f*dtFixed*dxinv[i];
  }

  for (i = is-1; i <= ie+2; i++) {
    qa    = dxi[i-2] + dxi[i-1] + dxi[i] + dxi[i+1];
    qb    = dxi[i-1]/(dxi[i-1] + dxi[i]);
    qc    = (dxi[i-2] + dxi[i-1])/(2.0*dxi[i-1] + dxi[i  ]);
    qd    = (dxi[i+1] + dxi[i  ])/(2.0*dxi[i  ] + dxi[i-1]);
    qb    = qb + 2.0*dxi[i]*qb/qa*(qc-qd);
    c3[i] = 1.0 - qb;
    c4[i] = qb;
    c5[i] =  dxi[i  ]/qa*qd;
    c6[i] = -dxi[i-1]/qa*qc;
    dx2i[i] = 0.5f*dxinv[i];
  }

//
//    Loop over sweep lines (in this slice)
//
  for (j = js; j <= je; j++) {
//
//     Precompute steepening coefficients if needed (eqns 1.14-1.17, plus 3.2)
//
    if (PPMSteepeningParameter) {
      for (i = is-2, index = in*j+is-2; i <= ie+2; i++, index++) {
	qa     = dxi[i-1] + dxi[i] + dxi[i+1];
	d2d[i] = (dslice[index+1] - dslice[index])/(dxi[i+1] + dxi[i]);
	d2d[i] = (d2d[i] - (dslice[index]-dslice[index-1])
		  /(dxi[i]+dxi[i-1]))/qa;
	dxb[i] = 0.5*(dxi[i] + dxi[i+1]);
      } // ENDFOR i
      for (i = is-1, index = in*j+is-1; i <= ie+1; i++, index++) {
	qc = fabs(dslice[index+1] - dslice[index-1])
	  - 0.01*min(fabs(dslice[index+1]), fabs(dslice[index-1]));
	s1 = (d2d[i-1] - d2d[i+1])*(powf(dxb[i-1],3) + powf(dxb[i],3))
	  /((dxb[i] + dxb[i-1])*
	    (dslice[index+1] - dslice[index-1] + tiny_number));
	if (d2d[i+1]*d2d[i-1] > 0.0) s1 = 0.0;
	if (qc <= 0.0) s1 = 0.0;
	s2 = max(0.0, min(20.0*(s1-0.05), 1.0));
	qa = fabs(dslice[index+1] - dslice[index-1])/
	  min(dslice[index+1],  dslice[index-1]);
	qb = fabs(pslice[index+1] - pslice[index-1])/
	  min(pslice[index+1],  pslice[index-1]);
	steepen[i] = (0.1*Gamma*qa >= qb) ? s2 : 0.0;
      } // ENDFOR i
    } // ENDIF steepen

//
//     Precompute left and right characteristic distances
//

    for (i = 0, index = in*j; i < in; i++, index++) {
      cs[i] = sqrtf(Gamma * pslice[index] / dslice[index]);
      if (PressureFree) cs[i] = tiny_number;
      char1[i] = max(0.0f,  dtFixed * (uslice[index] + cs[i])) * dx2i[i];
      char2[i] = max(0.0f, -dtFixed * (uslice[index] - cs[i])) * dx2i[i];
    } // ENDFOR i

    for (i = is-2, index = in*j+is-2; i <= ie+2; i++, index++) {
      cm[i] = dtFixed*(uslice[index]-cs[i])*(0.5f*dxinv[i]);
      c0[i] = dtFixed*(uslice[index]      )*(0.5f*dxinv[i]);
      cp[i] = dtFixed*(uslice[index]+cs[i])*(0.5f*dxinv[i]);
    } // ENDFOR i    

//
//     Compute left and right states for each variable
//       (steepening, if requested, is only applied to density)
//
    intvarC(dslice, in, is, ie, j, PPMSteepeningParameter, steepen,
	    PPMFlatteningParameter, flatten, c1, c2, c3, c4, c5, c6, char1, 
	    char2, c0, tmp1, tmp2, tmp3, tmp4, dla, dra, dl0, dr0);

    intvarC(pslice, in, is, ie, j, 0, steepen,
	    PPMFlatteningParameter, flatten, c1, c2, c3, c4, c5, c6, char1, 
	    char2, c0, dp, pl, pr, p6, pla, pra, pl0, pr0);

    intvarC(uslice, in, is, ie, j, 0, steepen,
	    PPMFlatteningParameter, flatten, c1, c2, c3, c4, c5, c6, char1, 
	    char2, c0, du, ul, ur, u6, ula, ura, ul0, ur0);

    intvarC(vslice, in, is, ie, j, 0, steepen,
	    PPMFlatteningParameter, flatten, c1, c2, c3, c4, c5, c6, char1, 
	    char2, c0, tmp1, tmp2, tmp3, tmp4, vla, vra, vl0, vr0);

    intvarC(wslice, in, is, ie, j, 0, steepen,
	    PPMFlatteningParameter, flatten, c1, c2, c3, c4, c5, c6, char1, 
	    char2, c0, tmp1, tmp2, tmp3, tmp4, wla, wra, wl0, wr0);

    if (DualEnergyFormalism)
      intvarC(geslice, in, is, ie, j, 0, steepen,
	      PPMFlatteningParameter, flatten, c1, c2, c3, c4, c5, c6, char1, 
	      char2, c0, tmp1, tmp2, tmp3, tmp4, gela, gera, gel0, ger0);

    for (ic = 0, index = 0, cindex = 0; ic < ncolors; 
	 ic++, index += in*jn, cindex += in)
      intvarC(colslice+index, in, is, ie, j, 0, steepen,
	      PPMFlatteningParameter, flatten, c1, c2, c3, c4, c5, c6, char1, 
	      char2, c0, tmp1, tmp2, tmp3, tmp4, 
	      colla+cindex, colra+cindex, coll0+cindex, colr0+cindex);
//
//
// Correct the initial guess from the linearized gas equations
//
//     First, compute averge over characteristic domain of dependance (3.5)
//
//
    for (i = is; i <= ie+1; i++) {
      plm[i]= pr[i-1]-cm[i-1]*(dp[i-1]-(1.0f-ft*cm[i-1])*p6[i-1]);
      prm[i]= pl[i  ]-cm[i  ]*(dp[i  ]+(1.0f+ft*cm[i  ])*p6[i  ]);
      plp[i]= pr[i-1]-cp[i-1]*(dp[i-1]-(1.0f-ft*cp[i-1])*p6[i-1]);
      prp[i]= pl[i  ]-cp[i  ]*(dp[i  ]+(1.0f+ft*cp[i  ])*p6[i  ]);
    }

    for (i = is; i <= ie+1; i++) {
      ulm[i]= ur[i-1]-cm[i-1]*(du[i-1]-(1.0f-ft*cm[i-1])*u6[i-1]);
      urm[i]= ul[i  ]-cm[i  ]*(du[i  ]+(1.0f+ft*cm[i  ])*u6[i  ]);
      ulp[i]= ur[i-1]-cp[i-1]*(du[i-1]-(1.0f-ft*cp[i-1])*u6[i-1]);
      urp[i]= ul[i  ]-cp[i  ]*(du[i  ]+(1.0f+ft*cp[i  ])*u6[i  ]);
    }
//
//     Compute correction terms (3.7)
//
    for (i = is; i <= ie+1; i++) {
      cla[i] = sqrtf(max(Gamma*pla[i]*dla[i], 0.0f));
      clainv[i] = 1.0f/cla[i];
      cra[i] = sqrtf(max(Gamma*pra[i]*dra[i], 0.0f));
      crainv[i] = 1.0f/cra[i];
    }
//
//     a) left side
//
    for (i = is; i <= ie+1; i++) {
      betalp[i] = (ula[i]-ulp[i]) + (pla[i]-plp[i])*clainv[i];
      betalm[i] = (ula[i]-ulm[i]) - (pla[i]-plm[i])*clainv[i];
      betal0[i] = (pla[i]-pl0[i])*clainv[i]*clainv[i] + 
	1.0f/dla[i] - 1.0f/dl0[i];
    }
//
//     Add gravity component
//      
    if (gravity)
      for (i = is, index = j*in+is; i <= ie+1; i++, index++) {
	betalp[i] = betalp[i] - 0.25f*dtFixed*(grslice[index-1]+grslice[index]);
	betalm[i] = betalm[i] - 0.25f*dtFixed*(grslice[index-1]+grslice[index]);
      } // ENDFOR i

    for (i = is; i <= ie+1; i++) {
      betalp[i] = -betalp[i]*(0.5f*clainv[i]);
      betalm[i] = +betalm[i]*(0.5f*clainv[i]);
    }

    for (i = is; i <= ie+1; i++) {
      if (cp[i-1] <= 0.0) betalp[i] = 0.0;
      if (cm[i-1] <= 0.0) betalm[i] = 0.0;
      if (c0[i-1] <= 0.0) betal0[i] = 0.0;
    }
//
//     b) right side
//
    for (i = is; i <= ie+1; i++) {
      betarp[i] = (ura[i]-urp[i]) + (pra[i]-prp[i])*crainv[i];
      betarm[i] = (ura[i]-urm[i]) - (pra[i]-prm[i])*crainv[i];
      betar0[i] = (pra[i]-pr0[i])*crainv[i]*crainv[i] + 
	1.0f/dra[i] - 1.0f/dr0[i];
    }
    
    if (gravity)
      for (i = is, index = j*in+is; i <= ie+1; i++, index++) {
	betarp[i] = betarp[i] - 0.25f*dtFixed*(grslice[index-1]+grslice[index]);
	betarm[i] = betarm[i] - 0.25f*dtFixed*(grslice[index-1]+grslice[index]);
      }

    for (i = is; i <= ie+1; i++) {
      betarp[i] = -betarp[i]*(0.5f*crainv[i]);
      betarm[i] = +betarm[i]*(0.5f*crainv[i]);
    }

    for (i = is; i <= ie+1; i++) {
      if (cp[i] > 0.0) betarp[i] = 0.0;
      if (cm[i] > 0.0) betarm[i] = 0.0;
      if (c0[i] > 0.0) betar0[i] = 0.0;
    }
//
//     Finally, combine to create corrected left/right states (eq. 3.6)
//
    for (i = is, index = j*in+is; i <= ie+1; i++, index++) {
      pls[index] = pla[i] + (betalp[i]+betalm[i])*cla[i]*cla[i];
      prs[index] = pra[i] + (betarp[i]+betarm[i])*cra[i]*cra[i];

      uls[index] = ula[i] + (betalp[i]-betalm[i])*cla[i];
      urs[index] = ura[i] + (betarp[i]-betarm[i])*cra[i];

      dls[index] = 1.0f/(1.0f/dla[i] - (betal0[i]+betalp[i]+betalm[i]));
      drs[index] = 1.0f/(1.0f/dra[i] - (betar0[i]+betarp[i]+betarm[i]));
    } // ENDFOR i
//
//     Take the appropriate state from the advected variables
//
    for (i = is, index = j*in+is; i <= ie+1; i++, index++) {
      if (uslice[index-1] <= 0.0) {
	vls[index] = vla[i];
	wls[index] = wla[i];
	gels[index] = gela[i];
      } else {
	vls[index] = vl0[i];
	wls[index] = wl0[i];
	gels[index] = gel0[i];
      }
      
      if (uslice[index] >= 0.0) {
	vrs[index] = vra[i];
	wrs[index] = wra[i];
	gers[index] = gera[i];
      } else {
	vrs[index] = vr0[i];
	wrs[index] = wr0[i];
	gers[index] = ger0[i];
      }
    } // ENDFOR i

    for (ic = 0; ic < ncolors; ic++) {
      cindex = (jn*ic+j)*in+is;
      index = j*in+is;
      scindex = ic*in+is;
      for (i = is; i <= ie+1; i++, scindex++, cindex++, index++) {
	if (uslice[index-1] <= 0.0)
	  colls[cindex] = colla[scindex];
	else
	  colls[cindex] = coll0[scindex];
	
	if (uslice[index] > 0.0)
	  colrs[cindex] = colra[scindex];
	else
	  colrs[cindex] = colr0[scindex];
      } // ENDFOR i
    } // ENDFOR ic
//
//  Dual energy formalism: if sound speed squared is less than eta1*v^2 
//    then discard the corrections and use pla, ula, dla.  This amounts
//     to assuming that we are outside the shocked region but the flow is
//     hypersonic so this should be true.  This is inserted because the
//     corrections are inaccurate for hypersonic flows.
//
    if (DualEnergyFormalism)
      for (i = is, index = j*in+is; i <= ie+1; i++, index++) {

	if (Gamma*pla[i]/dla[i] < DualEnergyFormalismEta2*ula[i]*ula[i] ||
	    max(max(fabsf(cm[i-1]), fabsf(c0[i-1])), fabsf(cp[i-1])) < 1e-3 ||
	    dls[index]/dla[i] > 5.0) {
	  pls[index] = pla[i];
	  uls[index] = ula[i];
	  dls[index] = dla[i];
	}

	if (Gamma*pra[i]/dra[i] < DualEnergyFormalismEta2*ura[i]*ura[i] ||
	    max(max(fabs(cm[i]), fabs(c0[i])), fabs(cp[i])) < 1e-3 ||
	    drs[index]/dra[i] > 5.0) {
	  prs[index] = pra[i];
	  urs[index] = ura[i];
	  drs[index] = dra[i];
	}

      } // ENDFOR i
//
//     Enforce minimum values.
//
    for (i = is, index = j*in+is; i <= ie+1; i++, index++) {
      pls[index] = max(pls[index], tiny_number);
      prs[index] = max(prs[index], tiny_number);
      dls[index] = max(dls[index], tiny_number);
      drs[index] = max(drs[index], tiny_number);
    }
//
//     If approximating pressure free conditions, then the density should be
//       reset to the pre-corrected state.
//
    if (PressureFree)
      for (i = is, index = j*in+is; i <= ie+1; i++, index++) {
	dls[index] = dla[i];
	drs[index] = dra[i];
      }

  } // ENDFOR j

  return SUCCESS;

}
