/* definitions for gadget cooling.
   created by Brian O'Shea, june 2002

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

#ifndef __Gadget_h_
#define __Gadget_h_

#define  NCOOLTAB  2000
#define  IONIZETABLE       /* turns on ionization table - comment this out
			     to use the ionization function */
#define  JAMPL	1.0		/* amplitude factor relative to input table */
#define  TABLESIZE 200		/* Max # of lines in TREECOOL */
#define  MAXITER 200
#define  SMALLNUM 1.0e-60
#define  GAMMA         (5.0/3.0)
#define  GAMMA_MINUS1  (GAMMA-1.0)
#define  MINGASTEMP 0.1
#define  HYDROGEN_MASSFRAC 0.76
#define  PROTONMASS  1.6726e-24
#define  BOLTZMANN   1.3806e-16

#endif  /* #ifndef __Gadget_h_ */

extern float *BetaH0, *BetaHep, *Betaff;
extern float *AlphaHp, *AlphaHep, *Alphad, *AlphaHepp;
extern float *GammaeH0, *GammaeHe0, *GammaeHep;
extern float nH0, nHp, nHep, nHe0, nHepp;
extern float bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
extern float J_UV, gJH0, gJHep, gJHe0, epsH0, epsHep, epsHe0;
extern float a0, planck, ev, e0_H, e0_He, e0_Hep;
extern float Tmin,Tmax,deltaT;
extern float nHcgs, XH;
extern float yhelium;
extern float ne;
extern float DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input;
extern float inlogz[TABLESIZE];
extern float gH0[TABLESIZE], gHe[TABLESIZE], gHep[TABLESIZE];
extern float eH0[TABLESIZE], eHe[TABLESIZE], eHep[TABLESIZE];
extern int nheattab;		/* length of table */
extern float ethmin;		/* minimum internal energy for neutral gas */
extern float mhboltz;		/* hydrogen mass over Boltzmann constant */
