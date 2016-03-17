#ifndef __INDIVIDUAL_STAR_PROPERTIES_H
#define __INDIVIDUAL_STAR_PROPERTIES_H

#include "typedefs.h"

#ifdef DEFINE_STORAGE
# define ISEXTERN
#else /*Define storage*/
# define ISEXTERN extern
#endif
                                         // INPUT - OUTPUT
// mass input always in solar masses, everything else always cgs (in and out)
float IndividualStarLifetime(float *mp); // solar - cgs
float IndividualStarLuminosity(float *mp); // solar - cgs
float IndividualStarLuminosity(float *mp, float *lifetime); // solar, cgs - cgs
float IndividualStarRadius(float *mp);                 // solar - cgs
float IndividualStarTeff(float *mp, float *lifetime);       // solar, cgs - cgs
float IndividualStarTeff(float *mp, float *L, float *R); // solar, cgs, cgs - cgs
float IndividualStarSurfaceGravity(float *mp); //solar - cgs
float IndividualStarSurfaceGravity(float *mp, float *R); // solar - cgs

int   IndividualStarComputeIonizingRates(float *q0, float *q1,
                                             float *Teff, float *g, float *metallicity);

int   IndividualStarInterpolateRadData(float *q0, float *q1,
                                       float *Teff, float *g, float *metallicity);
int PhotonRadianceBlackBody(float *q, float x);

int ComputeAverageEnergy(float *energy, float *e_i, float *Teff);
int AverageEnergyBlackBody(float *energy, float x);

float GaussianRandomVariable(void);

#endif
