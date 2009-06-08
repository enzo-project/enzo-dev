/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE TEST GRAVITY SPHERE PROBLEM
/
/  written by: Greg Bryan
/  date:       September, 1995
/  modified1:
/
/  PURPOSE:
/    This is global data that pertains only to the TestGravitySphere problem.
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define TGEXTERN
#else /* DEFINE_STORAGE */
# define TGEXTERN extern
#endif /* DEFINE_STORAGE */

/* Inside and outside densities. */

TGEXTERN float TestGravitySphereInteriorDensity;
TGEXTERN float TestGravitySphereExteriorDensity;

/* Sphere radius. */

TGEXTERN float TestGravitySphereRadius;

/* Type of sphere:
   0 - uniform inside and outside
   1 - r^-2 power law inside and uniform outside
   1 - r^-2.25 power law inside and uniform outside */

TGEXTERN int TestGravitySphereType;
