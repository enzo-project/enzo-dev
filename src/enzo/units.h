/* Define units in CGS */

#ifndef UNITS_DEFINED__
#define UNITS_DEFINED__


#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* Units of length is centimetres */

EXTERN float GlobalLengthUnits;

/* Units of mass in grams */

EXTERN double GlobalMassUnits;

/* Units of density in g/cm^3 */

EXTERN float GlobalDensityUnits;

/* Units of time in seconds */

EXTERN float GlobalTimeUnits;

#endif
