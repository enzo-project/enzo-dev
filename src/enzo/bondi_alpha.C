#include "preincludes.h"
#include "macros_and_parameters.h"


/* Routine to return alpha, defined as rho/rho_inf, for a critical
   Bondi accretion solution.  The argument is x = r / rBondiHoyle. 
   Adapted from Orion, courtesy of Mark Krumholz */

float bondi_alpha(float x) {

#define XMIN 0.01
#define XMAX 2.0
#define NTABLE 51

  float lambda_c, xtable, xtablep1, alpha_exp;
  int idx;

  /* This is a precomputed table of alpha values.  These correspond to x values
     that run from 0.01 to 2.0 with uniform logarithmic spacing.  The reason for
     this choice of range is that the asymptotic expressions are accurate to
     better than 2% outside this range */

  float alphatable[NTABLE] = {820.254, 701.882, 600.752, 514.341, 440.497, 377.381, 323.427,
			      277.295, 237.845, 204.1, 175.23, 150.524, 129.377, 111.27, 95.7613,
			      82.4745, 71.0869, 61.3237, 52.9498, 45.7644, 39.5963, 34.2989,
			      29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, 13.0714,
			      11.4717, 10.0903, 8.89675, 7.86467, 6.97159, 6.19825, 5.52812,
			      4.94699, 4.44279, 4.00497, 3.6246, 3.29395, 3.00637, 2.75612,
			      2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, 1.70804,
			      1.62439};

  // A constant that appears in the following formulae.  This hardcoded value is
  // valid for an isothermal gas.
  lambda_c = 0.25*exp(1.5);

  // deal with the off-the-table cases
  if (x < XMIN) 
    return lambda_c / sqrt(2.*x*x);
  else if (x >= XMAX)
    return exp(1./x);
  else {
    // we are on the table
    
    idx = floor((NTABLE-1)*log(x/XMIN)/log(XMAX/XMIN));
    xtable = exp(log(XMIN) + idx*log(XMAX/XMIN)/(NTABLE-1));
    xtablep1 = exp(log(XMIN) + (idx+1)*log(XMAX/XMIN)/(NTABLE-1));
    alpha_exp = log(x/xtable) / log(xtablep1/xtable);

    return alphatable[idx] * POW(alphatable[idx+1]/alphatable[idx],alpha_exp);
  }

#undef NTABLE
#undef XMIN
#undef XMAX

}
