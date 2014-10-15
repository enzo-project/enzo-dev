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

//
// generate couple of normally distributed random deviates (Box-Muller-Algorithm)
//
void StochasticForcing::GaussDeviate(float amplt, float *x, float *y)
{
	double v_sqr, v1, v2;
	float norm;

	do {
	    v1 = 2.0*RandUni(seed) - 1.0;
	    v2 = 2.0*RandUni(seed) - 1.0;
	    v_sqr = v1*v1+v2*v2;
	} while (v_sqr >= 1.0 || v_sqr == 0.0);
	
	norm = amplt * sqrt(-2.0*log(v_sqr)/v_sqr);

	*x = norm * v1; *y = norm * v2;
}

//
// random deviate of uniform distribution (L'Ecuyer with Bays-Durham shuffle)
//
double StochasticForcing::RandUni(int &idum)
{
	const int IM1=2147483563,IM2=2147483399;
	const int IA1=40014,IA2=40692,IQ1=53668,IQ2=52774;
	const int IR1=12211,IR2=3791,NTAB=32,IMM1=IM1-1;
	const int NDIV=1+IMM1/NTAB;
	const double EPS=3.0e-16,RNMX=1.0-EPS,AM=1.0/double(IM1);
	int j,k;
	double temp;

	if (idum <= 0) {
		idum=(idum==0 ? 1 : -idum);
		idum2=idum;
		for (j=NTAB+7;j>=0;j--) {
			k=idum/IQ1;
			idum=IA1*(idum-k*IQ1)-k*IR1;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}
	k=idum/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;
	if (idum < 0) idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
