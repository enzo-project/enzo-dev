/***********************************************************************
/
/  GRID CLASS (SOLVE THE GADGET COOLING/HEATING RATE EQUATIONS)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/           given an internal energy density, return a temperature.
/  RETURNS:
/           temperature, in kelvin
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

void Gadgetfind_abundances_and_rates(float logT, float rho, float *ne_guess);

float grid::Gadgetconvert_u_to_temp(float u, float rho, float *ne_guess)
{

  float temp, temp_old, temp_new, max = 0, ne_old;
  float mu;
  int iter = 0;

  float u_input, rho_input, ne_input;

  //if(debug) printf("In Gadgetconvert_u_to_temp:1\n");


  u_input = u;
  rho_input = rho;
  ne_input = *ne_guess;

  mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);


  temp = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

  //if(debug) printf("GAMMA_MINUS1: %lf  BOLTZMANN:  %e  u:  %e\n",GAMMA_MINUS1, BOLTZMANN, u);
  //if(debug) printf("PROTONMASS:  %e  mu:  %e\n",PROTONMASS,mu);
  //if(debug) printf("temp:  %e log10(temp):  %lf\n",temp,log10(temp));
  //if(debug) printf("In Gadgetconvert_u_to_temp:2\n");

  do
    {
      ne_old = *ne_guess;

      //if(debug) printf("In Gadgetconvert_u_to_temp:3\n");
      Gadgetfind_abundances_and_rates(log10(temp), rho, ne_guess);
      temp_old = temp;
      //if(debug) printf("In Gadgetconvert_u_to_temp:4\n");
      mu = (1 + 4 * yhelium) / (1 + yhelium + *ne_guess);

      temp_new = GAMMA_MINUS1 / BOLTZMANN * u * PROTONMASS * mu;

      max =
	max(max,
	     temp_new / (1 + yhelium + *ne_guess) * fabs((*ne_guess - ne_old) / (temp_new - temp_old + 1.0)));

      temp = temp_old + (temp_new - temp_old) / (1 + max);
      iter++;
      //if(debug) printf("In Gadgetconvert_u_to_temp:5\n");
      if(iter > (MAXITER - 10))
	printf("-> temp= %e ne=%e\n", temp, *ne_guess);
    }
  while(fabs(temp - temp_old) > 1.0e-3 * temp && iter < MAXITER);

  if(iter >= MAXITER)   /* if it breaks! */
    {
      printf("failed to converge in convert_u_to_temp()\n");
      printf("u_input= %e\nrho_input=%e\n ne_input=%e\n", u_input, rho_input, ne_input);
      printf("DoCool_u_old_input=%e\nDoCool_rho_input= %e\nDoCool_dt_input= %e\nDoCool_ne_guess_input= %e\n",
	     DoCool_u_old_input, DoCool_rho_input, DoCool_dt_input, DoCool_ne_guess_input);
    }  /* need to think about how to do error-checking here! */

  return temp;
}
