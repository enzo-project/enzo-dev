/***********************************************************************
 *
 *  STOCHASTIC FORCING CLASS: GaussDeviate, RandUni
 *
 *  written by: Wolfram Schmidt
 *  date:       May, 2005
 *  modified1: Oct, 2014: updated to support Enzo 2.4 // P. Grete
 *
 *  PURPOSE: random number generator producing Gaussian deviates
 *           (adapted from Numerical Recipes) 
 *
 ***********************************************************************/

#include "preincludes.h"
#include "StochasticForcing.h"

unsigned_long_int mt_random();

//
// generate couple of normally distributed random deviates (Box-Muller-Algorithm)
//
void StochasticForcing::GaussDeviate(float amplt, float *x, float *y)
{
	float v_sqr, v1, v2;
	float norm;

	do {
	    v1 = 2.0* (float)(mt_random()%2147483563)/(2147483563.0) - 1.0;
	    v2 = 2.0* (float)(mt_random()%2147483563)/(2147483563.0) - 1.0;
	    v_sqr = v1*v1+v2*v2;
	} while (v_sqr >= 1.0 || v_sqr == 0.0);
	
	norm = amplt * sqrt(-2.0*log(v_sqr)/v_sqr);

	*x = norm * v1; *y = norm * v2;
}

