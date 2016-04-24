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
float IndividualStarLifetime(const float &mp, const float &metallicity); // solar, fraction - cgs
float IndividualStarLuminosity(const float &mp, const float &lifetime); // solar, cgs - cgs
float IndividualStarSurfaceGravity(const float &mp, const float &R);


/* Functions for stellar properties data */
int IndividualStarInterpolateLuminosity(float &L, const float &M, const float &metallicity);
int IndividualStarInterpolateProperties(float *Teff, float *R,
                                        const float &M, const float &metallicity);


/* Functions for radiation data */
int IndividualStarComputeIonizingRates(float *q0, float *q1,
                                       float *Teff, float *g, float *metallicity);

int IndividualStarInterpolateRadData(float *q0, float *q1,
                                     float *Teff, float *g, float *metallicity);
int PhotonRadianceBlackBody(float *q, float x);

int ComputeAverageEnergy(float *energy, float *e_i, float *Teff);
int AverageEnergyBlackBody(float *energy, float x);


/* General helper functions */
float GaussianRandomVariable(void);

#endif
