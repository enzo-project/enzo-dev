/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE SPHERICAL INFALL TEST
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define SIEXTERN
#else /* DEFINE_STORAGE */
# define SIEXTERN extern
#endif /* DEFINE_STORAGE */

/* Flag indicating whether to use fixed perturbation mass. */

SIEXTERN int SphericalInfallFixedAcceleration;

/* Perturbation mass. */

SIEXTERN float SphericalInfallFixedMass;

/* Center of peturbation. */

SIEXTERN FLOAT SphericalInfallCenter[MAX_DIMENSION];
