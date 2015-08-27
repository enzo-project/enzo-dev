/***********************************************************************
 *
 *  STOCHASTIC FORCING CLASS: Evolve, Inject, WriteSpectrum
 *
 *  written by: Wolfram Schmidt
 *  date:       May, 2005
 *  modified1: Oct, 2014: updated to support Enzo 2.4 // P. Grete
 *
 *  PURPOSE: evolves the random forcing spectrum in the fashion of
 *           a multi-dimensional Ornstein-Uhlenbeck process
 *           
 *           Parameters:
 *           dt -- time step (small compared to AutoCorrlTime)
 *
 ***********************************************************************/

#include "preincludes.h"
#include "StochasticForcing.h"
#include "global_data.h"


void StochasticForcing::Evolve(float dt)
{
    if (MyProcessorNumber == ROOT_PROCESSOR) {

	float DriftCoeff[MAX_DIMENSION], DiffCoeff[MAX_DIMENSION];

	if (decay == 0) {

	    Inject();
		    
	    /* Increment forcing spectrum (drift and random diffusion) 
         * For general properties of Ornstein-Uhlenbeck process, see e.g.
         * Turbulent Flows by Pope (2000) Appendix J with 
         * drift and diffusion coefficients given eq (J.41)
         * */

	    for (int dim = 0; dim < SpectralRank; dim++) {
		DriftCoeff[dim] = exp(-dt/AutoCorrlTime[dim]);
		DiffCoeff[dim]  = sqrt(1 - DriftCoeff[dim]*DriftCoeff[dim]);
		for (int n = 0, m = 0; n < NumModes; n++)
		    if (mask[n]) {
			SpectrumEven[dim][m] = DriftCoeff[dim] * SpectrumEven[dim][m] +
			                       DiffCoeff [dim] * InjectionEven[dim][n];
			SpectrumOdd [dim][m] = DriftCoeff[dim] * SpectrumOdd [dim][m] + 
			                       DiffCoeff [dim] * InjectionOdd [dim][n];
			++m;
		    }
	    }

	} else {

	    /* increment forcing spectrum (drift only) */

	    for (int dim = 0; dim < SpectralRank; dim++) {
		DriftCoeff[dim] = exp(-dt/AutoCorrlTime[dim]);
		for (int m = 0; m < NumNonZeroModes; m++) {
		    SpectrumEven[dim][m] = DriftCoeff[dim] * SpectrumEven[dim][m];
		    SpectrumOdd [dim][m] = DriftCoeff[dim] * SpectrumOdd [dim][m];
		}
	    }
	}
    }

    CommunicationBroadcastSpectrum();

}

//
// Compute new random injection
//
void StochasticForcing::Inject(void)
{
    if (MyProcessorNumber == ROOT_PROCESSOR) {

	int i, j, k, n, dim;
	float a, b, contr, div;

	/* compute Gaussian deviates */

	for (dim = 0; dim < SpectralRank; dim++)
	    for (n = 0; n < NumModes; n++) {
		if (mask[n]) {
		    GaussDeviate(Amplitude[dim][n], &a, &b);
		} else {
		    a = 0.0; b = 0.0;
		}
		InjectionEven[dim][n] = a;
		InjectionOdd[dim][n]  = b;
	    }

	/* project modes 
     * see eq (8) in Schmidt et al., A&A (2009)
     * http://dx.doi.org/10.1051/0004-6361:200809967 */

	for (i = 0; i < i2; i++) { // wave vectors in positive x-direction
	    InjectionEven[0][i] = (1.0 - SolenoidalWeight) * InjectionEven[0][i];
	    InjectionOdd[0][i]  = (1.0 - SolenoidalWeight) * InjectionOdd[0][i];
	}

	if (SpectralRank > 1) {
	    
	    for (n = 0; n < i2; n++) { // wave vectors in positive x-direction
		InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n];
		InjectionOdd [1][n] = SolenoidalWeight * InjectionOdd [1][n];
	    }
	    
	    n = i2;
	    for (j = 1; j <= j2; j++) { // wave vectors in xy-plane
		for (i = i1; i <= i2; i++) {
		    contr = (1.0 - 2.0 * SolenoidalWeight) * 
			(i*InjectionEven[0][n] + 
			 j*InjectionEven[1][n]) / float(i*i + j*j);
		    InjectionEven[0][n] = SolenoidalWeight * InjectionEven[0][n] + i*contr;
		    InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n] + j*contr;
		    contr = (1.0 - 2.0 * SolenoidalWeight) * 
			(i*InjectionOdd[0][n] + 
			 j*InjectionOdd[1][n]) / float(i*i + j*j);
		    InjectionOdd[0][n] = SolenoidalWeight * InjectionOdd[0][n] + i*contr;
		    InjectionOdd[1][n] = SolenoidalWeight * InjectionOdd[1][n] + j*contr;
		    ++n;
		}
	    }
	    
	    if (SpectralRank > 2) {
		
		for (n = 0; n < i2 + j2*(i2-i1+1); n++) { // wave vectors in xy-plane
		    InjectionEven[2][n] = SolenoidalWeight * InjectionEven[2][n];
		    InjectionOdd[2][n]  = SolenoidalWeight * InjectionOdd [2][n];
		}
		
		for (k = 1; k <= k2; k++) { // wave vectors not aligned to xy-plane
		    for (j = j1; j <= j2; j++) {
			for (i = i1; i <= i2; i++) {
			    contr = (1.0 - 2.0 * SolenoidalWeight) * 
				(i*InjectionEven[0][n] + 
				 j*InjectionEven[1][n] + 
				 k*InjectionEven[2][n] ) / float(i*i + j*j + k*k);
			    InjectionEven[0][n] = SolenoidalWeight * InjectionEven[0][n] + i*contr;
			    InjectionEven[1][n] = SolenoidalWeight * InjectionEven[1][n] + j*contr;
			    InjectionEven[2][n] = SolenoidalWeight * InjectionEven[2][n] + k*contr;
			    contr = (1.0 - 2.0 * SolenoidalWeight) * 
				(i*InjectionOdd[0][n] + 
				 j*InjectionOdd[1][n] +
				 k*InjectionOdd[2][n]) / float(i*i + j*j + k*k);
			    InjectionOdd[0][n] = SolenoidalWeight * InjectionOdd[0][n] + i*contr;
			    InjectionOdd[1][n] = SolenoidalWeight * InjectionOdd[1][n] + j*contr;
			    InjectionOdd[2][n] = SolenoidalWeight * InjectionOdd[2][n] + k*contr;
			    ++n;

			}
		    }
		}
	    }
	}

    }
}

float StochasticForcing::RMS(void)
{
    int m;
    float sum_even = 0.0, sum_odd = 0.0;

    for (int dim = 0; dim < SpectralRank; dim++) {
        for (m = 0; m < NumNonZeroModes; m++)
            sum_even += SpectrumEven[dim][m] * SpectrumEven[dim][m];
        for (m = 0; m < NumNonZeroModes; m++)
            sum_odd  += SpectrumOdd[dim][m]  * SpectrumOdd[dim][m];
    }

    return sqrt(sum_even + sum_odd);
}
