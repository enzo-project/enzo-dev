/* declarations for gadget cooling.
   created by Brian O'Shea, june 2002
   modified by James Bordner 2003-10-30

   note that gadget cooling is turned on/off by UseGadgetCooling.

/  NOTE:  This routine is based on code in the GADGET code, written
/    by Volker Springel and Naoki Yoshida.  The subroutines are
/    essentially identical to the Gadget cooling routines, with an
/    Enzo-friendly wrapper around them.  As of the addition of this
/    comment to the code (Feb. 2004) the cooling code is private and
/    NOT FOR PUBLIC RELEASE until Volker makes the cooling code public
/    in Gadget.  See the Gadget web page for more details:
/
/      http://www.mpa-garching.mpg.de/~volker/gadget/
/
/      --B.W.O'Shea, Feb. 2004

*/

#include "macros_and_parameters.h"
#include "Gadget.h"

 float *BetaH0, *BetaHep, *Betaff;
 float *AlphaHp, *AlphaHep, *Alphad, *AlphaHepp;
 float *GammaeH0, *GammaeHe0, *GammaeHep;
 float nH0, nHp, nHep, nHe0, nHepp;
 float bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
 float J_UV, gJH0, gJHep, gJHe0, epsH0, epsHep, epsHe0;
 float a0, planck, ev, e0_H, e0_He, e0_Hep;
 float Tmin,Tmax,deltaT;
 float nHcgs, XH;
 float yhelium;
 float ne;
 float DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input;
 float inlogz[TABLESIZE];
 float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
 float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
 int nheattab;                /* length of table */
 float ethmin;                /* minimum internal energy for neutral gas */
 float mhboltz;               /* hydrogen mass over Boltzmann constant */





