/***********************************************************************
 *  
 *  STOCHASTIC FORCING CLASS: destructor
 *
 *  written by: Wolfram Schmidt
 *  date:       May, 2005
 *  modified1:
 *
 ***********************************************************************/
#include "preincludes.h"
#include "StochasticForcing.h"
//
//  Default constructor (Set all data to null/default state).
//
StochasticForcing::~StochasticForcing() 
{
    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
	delete [] Amplitude[dim];
	delete [] InjectionEven[dim];
	delete [] InjectionOdd[dim];
	delete [] SpectrumEven[dim];
	delete [] SpectrumOdd[dim];
    }
	delete [] mask;
}
