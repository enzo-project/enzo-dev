/***********************************************************************
/
/  GRID CLASS (SOLVE THE GADGET COOLING/HEATING RATE EQUATIONS)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/        given temp, density, guess at electron number density, solve for
/        the abundances of various species and their cooling/heating rates
/
/  RETURNS:
/    void - no return
/
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
/
************************************************************************/

#include <stdio.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "Gadget.h"

/*float necgs, nHcgs;
  float DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input;
  float bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep;
  float nH0,nHe0,nHp,nHep,nHepp,ne;
*/
void grid::Gadgetfind_abundances_and_rates(float logT, float rho, float *ne_guess)
{

  float neold, nenew;
  int j, niter;
  float Tlow, Thi, flow, fhi, t;
  float necgs;
  float logT_input, rho_input, ne_input;


  float gJH0ne, gJHe0ne, gJHepne;

  logT_input = logT;
  rho_input = rho;
  ne_input = *ne_guess;

  //if(debug) printf("In Gadget_find_abundances_and_rates:1\n");

  if(logT <= Tmin)		/* everything neutral */
    {
      nH0 = 1.0;
      nHe0 = yhelium;
      nHp = 0;
      nHep = 0;
      nHepp = 0;
      ne = 0;
      *ne_guess = 0;
      //if(debug) printf("In Gadget_find_abundances_and_rates:2\n");
      return;
    }
  //if(debug) printf("In Gadget_find_abundances_and_rates:3\n");
  if(logT >= Tmax)		/* everything is ionized */
    {
      nH0 = 0;
      nHe0 = 0;
      nHp = 1.0;
      nHep = 0;
      nHepp = yhelium;
      ne = nHp + 2.0 * nHepp;
      *ne_guess = ne;		/* note: in units of the hydrogen number density */
      //if(debug) printf("In Gadget_find_abundances_and_rates:4\n");
      return;
    }
  //if(debug) printf("In Gadget_find_abundances_and_rates:5\n");
  //if(debug) printf("logT:  %lf  Tmin: %lf  deltaT:  %lf\n",logT,Tmin,deltaT);
  t = (logT - Tmin) / deltaT;
  j = (int) t;
  //if(debug) printf("t:  %lf   j: %i\n",t,j);
  Tlow = Tmin + deltaT * j;
  Thi = Tlow + deltaT;
  fhi = t - j;
  flow = 1 - fhi;

  if(*ne_guess == 0)
    *ne_guess = 1.0;

  //if(debug) printf("In Gadget_find_abundances_and_rates:6\n");

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

  ne = *ne_guess;
  neold = ne;
  niter = 0;
  necgs = ne * nHcgs;

  /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
  do
    {
      niter++;
      //if(debug) printf("In Gadget_find_abundances_and_rates:7\n");
      //if(debug) printf("value of j is %i\n",j);
      aHp = flow * AlphaHp[j] + fhi * AlphaHp[j + 1];
      //if(debug) printf("In Gadget_find_abundances_and_rates:7.05\n");
      aHep = flow * AlphaHep[j] + fhi * AlphaHep[j + 1];
      aHepp = flow * AlphaHepp[j] + fhi * AlphaHepp[j + 1];
      ad = flow * Alphad[j] + fhi * Alphad[j + 1];
      geH0 = flow * GammaeH0[j] + fhi * GammaeH0[j + 1];
      geHe0 = flow * GammaeHe0[j] + fhi * GammaeHe0[j + 1];
      geHep = flow * GammaeHep[j] + fhi * GammaeHep[j + 1];
      //if(debug) printf("In Gadget_find_abundances_and_rates:7.1\n");
      if(necgs <= 1.e-25 || J_UV == 0)
	{
	  gJH0ne = gJHe0ne = gJHepne = 0;
	}
      else
	{
	  gJH0ne = gJH0 / necgs;
	  gJHe0ne = gJHe0 / necgs;
	  gJHepne = gJHep / necgs;
	}
      //if(debug) printf("In Gadget_find_abundances_and_rates:7.2\n");
      nH0 = aHp / (aHp + geH0 + gJH0ne);	/* eqn (33) */
      nHp = 1.0 - nH0;		/* eqn (34) */

      if((gJHe0ne + geHe0) <= SMALLNUM)	/* no ionization at all */
	{
	  nHep = 0.0;
	  nHepp = 0.0;
	  nHe0 = yhelium;
	}
      else
	{
	  nHep = yhelium / (1.0 + (aHep + ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	/* eqn (35) */
	  nHe0 = nHep * (aHep + ad) / (geHe0 + gJHe0ne);	/* eqn (36) */
	  nHepp = nHep * (geHep + gJHepne) / aHepp;	/* eqn (37) */
	}
      //if(debug) printf("In Gadget_find_abundances_and_rates:7.3\n");
      neold = ne;

      //if(debug) printf("In Gadget_find_abundances_and_rates:8\n");

      ne = nHp + nHep + 2 * nHepp;	/* eqn (38) */
      necgs = ne * nHcgs;

      if(J_UV == 0)
	break;

      nenew = 0.5 * (ne + neold);
      ne = nenew;
      necgs = ne * nHcgs;

      if(fabs(ne - neold) < 1.0e-4)
	break;

      if(niter > (MAXITER - 10))
	printf("ne= %e  niter=%"ISYM"\n", ne, niter);
      //if(debug) printf("In Gadget_find_abundances_and_rates:9\n");
    }
  while(niter < MAXITER);

  if(niter >= MAXITER)
    {
      printf("no convergence reached in find_abundances_and_rates()\n");
      printf("logT_input= %e  rho_input= %e  ne_input= %e\n", logT_input, rho_input, ne_input);
      printf("DoCool_u_old_input=%e\nDoCool_rho_input= %e\nDoCool_dt_input= %e\nDoCool_ne_guess_input= %e\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
      /* endrun(13); */
    }

  bH0 = flow * BetaH0[j] + fhi * BetaH0[j + 1];
  bHep = flow * BetaHep[j] + fhi * BetaHep[j + 1];
  bff = flow * Betaff[j] + fhi * Betaff[j + 1];

  *ne_guess = ne;

}
