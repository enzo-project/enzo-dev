/***********************************************************************
 *  
 *  STOCHASTIC FORCING CLASS: constructor
 *
 *  written by: Wolfram Schmidt
 *  date:       May, 2005
 *  modified1:  Oct. 2014: updated to support Enzo 2.4 // P. Grete
 *
 ***********************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "StochasticForcing.h"
#define DEBUG

//
//  Default constructor (Set all data to null/default state).
//
StochasticForcing::StochasticForcing() 
{
    i1 = i2 = 0;
    j1 = j2 = 0;
    k1 = k2 = 0;
    NumModes = 0;
    decay = 0;    
    
    idum2=123456789;
    iy=0;

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
	alpha[dim]         = 0;         
	IntgrVelocity[dim] = 0.0; 
	IntgrLength[dim]   = 0.0;
	WaveNumber[dim]    = 0.0;
	IntgrTime[dim]     = 0.0;
	AutoCorrlTime[dim] = 0.0;

	mask               = NULL;
	Amplitude[dim]     = NULL;
	InjectionEven[dim] = NULL;
	InjectionOdd[dim]  = NULL;
	SpectrumEven[dim]  = NULL;
	SpectrumOdd[dim]   = NULL;
    }

    SpectProfile = Peak;
    SolenoidalWeight = 1.0;
}

