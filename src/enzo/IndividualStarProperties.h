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
float IndividualStarSurfaceGravity(const float &mp, const float &R);

/* Functions for stellar properties data */
int IndividualStarInterpolateLuminosity(float &L, const float &M, const float &metallicity);
int IndividualStarInterpolateProperties(float &Teff, float &R,
                                        const float &M, const float &metallicity);

int IndividualStarInterpolateLifetime(float & tau, const float &M, const float &metallicity, const int &mode);


/* Functions for radiation data */
int IndividualStarComputeIonizingRates(float &q0, float &q1,
                                       const float &Teff, const float &g, const float &metallicity);

int IndividualStarInterpolateRadData(float &q0, float &q1,
                                     const float &Teff, const float &g, const float &metallicity);

int PhotonRadianceBlackBody(float &q, const float &x);

int ComputeAverageEnergy(float *energy, float *e_i, float *Teff);
int AverageEnergyBlackBody(float *energy, float x);
int ComputeBlackBodyFlux(float &flux, const float &Teff, const float &e_min, const float &e_max);


int BlackBodyFlux(float &F, const float &x1, const float &x2);
int BlackBodyFlux(float &F, const float &x);



/* General helper functions */
float GaussianRandomVariable(void);

#endif
