//=======================================================================
/////////////////////////  SUBROUTINE INTVAR  \\\\\\\\\\\\\\\\\\\\\\\\\\\
//
//      subroutine intvar(qslice, idim, i1, i2, isteep, steepen, iflatten, 
//     &                  flatten, c1, c2, c3, c4, c5, c6, char1, char2,
//     &                  c0, dq, ql, qr, q6, qla, qra, ql0, qr0)
//
//  COMPUTES LEFT AND RIGHT EULERIAN INTERFACE VALUES FOR RIEMANN SOLVER
//
//  written by: Greg Bryan
//  date:       March, 1996
//  modified1:  January, 2010 by JHW (translated to C++)
//
//  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
//    and right interface values to be fed into Riemann solver during a
//    one dimensional sweeps.  This version computes the Eulerian corrections
//    to the left and right states described in section three of Colella &
//    Woodward (1984), JCP.  The routine works on a single variable in
//    one dimension.
//
//  INPUT:
//    qslice   - one dimensional field of quantity q (one of d,e,u,v...)
//    idim     - declared dimension of 1D fields
//    i1, i2   - start and end indexes of active region
//    isteep   - steepening flag (1 = on, 0 = off); only apply to density!
//    steepen    - steepening coefficients
//    iflatten - flattening flag (1 = on, 0 = off)
//    flatten  - flattening coefficients
//    c1-6     - precomputed grid coefficients
//    char1,2  - characteristic distances for +/- waves (for average)
//    c0       - characteristic distance (for lagrangean cell face)
//    dq, ql, qr, q6 - 1D field temporaries
//    
//  OUTPUT:
//    qla, qra - left and right state values (from char1,2)
//    ql0, qr0 - left and right state values (from c0)
//
//  EXTERNALS:
//
//  LOCALS:
//
//  PARAMETERS:
//    ft     - a constant used in eq. 1.124 (=2*2/3)
//
//-----------------------------------------------------------------------

#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"

void intvarC(float *qslice, int in, int is, int ie, int j, int isteep,
	     float *steepen, int iflatten, float *flatten, float *c1,
	     float *c2, float *c3, float *c4, float *c5, float *c6,
	     float *char1, float *char2, float *c0, float *dq, float *ql,
	     float *qr, float *q6, float *qla, float *qra, float *ql0,
	     float *qr0)
{

  const float ft = 4.0/3.0;
  int i, index, indexm1, indexm2;
  float qplus, qmnus, qvanl, temp1, temp2, temp3, temp22, temp23;
  indexm1 = j*in+is-1;
  indexm2 = indexm1-1;
//
//     Compute average linear slopes (eqn 1.7)
//      Monotonize (eqn 1.8)
//
  for (i = is-2, index = indexm2; i <= ie+2; i++, index++) {
    qplus = qslice[index+1] - qslice[index];
    qmnus = qslice[index  ] - qslice[index-1];
    qvanl = 2.0 * qplus * qmnus / (qmnus + qplus);
    dq[i] = c1[i] * qplus + c2[i] * qmnus;
    temp1 = min(min(min(fabsf(dq[i]), 2.0f*fabsf(qmnus)), 2.0f*fabsf(qplus)),
		fabsf(qvanl));
    dq[i] = (qplus*qmnus > 0) ? temp1*sign(dq[i]) : 0.0;
  }
//     
//     Construct left and right values (eqn 1.6)
//
  for (i = is-1, index = indexm1; i <= ie+2; i++, index++) {
    ql[i] = c3[i] * qslice[index-1] + c4[i] * qslice[index] +
      c5[i] * dq[i-1] + c6[i] * dq[i];
    qr[i-1] = ql[i];
  }
//
//     Steepen if asked for (use precomputed steepening parameter)
//
  if (isteep)
    for (i = is-1, index = indexm1; i <= ie+1; i++, index++) {
      ql[i] = (1.0f-steepen[i])*ql[i] + 
	steepen[i]*(qslice[index-1] + 0.5f*dq[i-1]);
      qr[i] = (1.0f-steepen[i])*qr[i] +
	steepen[i]*(qslice[index+1] + 0.5f*dq[i+1]);
    }
//
//     Monotonize again (eqn 1.10)
//
  for (i = is-1, index = indexm1; i <= ie+1; i++, index++) {
    temp1 = (qr[i]-qslice[index])*(qslice[index]-ql[i]);
    temp2 = qr[i]-ql[i];
    temp3 = 6.0f*(qslice[index]-0.5f*(qr[i]+ql[i]));
    if (temp1 < 0.0) {
      ql[i] = qslice[index];
      qr[i] = qslice[index];
    }
    temp22 = temp2*temp2;
    temp23 = temp2*temp3;
    if (temp22 < temp23)  ql[i] = 3.0f*qslice[index] - 2.0f*qr[i];
    if (temp22 < -temp23) qr[i] = 3.0f*qslice[index] - 2.0f*ql[i];
  } // ENDFOR i
//
//     If requested, flatten slopes with flatteners calculated in calcdiss (4.1)
//
  if (iflatten)
    for (i = is-1, index = indexm1; i <= ie+1; i++, index++) {
      ql[i] = qslice[index]*flatten[i] + ql[i]*(1.0f-flatten[i]);
      qr[i] = qslice[index]*flatten[i] + qr[i]*(1.0f-flatten[i]);
    }
//
//     Ensure that the L/R values lie between neighboring cell-centered 
//     values (Taken from ATHENA, lr_states)
//
#define CHECK_LR
#ifdef CHECK_LR
  for (i = is-1, index = indexm1; i <= ie+1; i++, index++) {
    ql[i] = max(min(qslice[index], qslice[index-1]), ql[i]);
    ql[i] = min(max(qslice[index], qslice[index-1]), ql[i]);
    qr[i] = max(min(qslice[index], qslice[index+1]), qr[i]);
    qr[i] = min(max(qslice[index], qslice[index+1]), qr[i]);
  }
#endif
//
//    Now construct left and right interface values (eqn 1.12 and 3.3)
//
  for (i = is-1, index = indexm1; i <= ie+1; i++, index++) {
    q6[i] = 6.0f*(qslice[index]-0.5f*(ql[i]+qr[i]));
    dq[i] = qr[i] - ql[i];
  }

  for (i = is; i <= ie+1; i++) {
    qla[i]= qr[i-1]-char1[i-1]*(dq[i-1]-(1.0f-ft*char1[i-1])*q6[i-1]);
    qra[i]= ql[i  ]+char2[i  ]*(dq[i  ]+(1.0f-ft*char2[i  ])*q6[i  ]);
  }

  for (i = is; i <= ie+1; i++) {
    ql0[i] = qr[i-1]-c0[i-1]*(dq[i-1]-(1.0f-ft*c0[i-1])*q6[i-1]);
    qr0[i] = ql[i  ]-c0[i  ]*(dq[i  ]+(1.0f+ft*c0[i  ])*q6[i  ]);
  }

  return;

}
