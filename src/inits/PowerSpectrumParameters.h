/* Power Spectrum Parameters */

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* Setable parameters. */

EXTERN FLOAT sigma8;
EXTERN FLOAT PrimordialIndex;
EXTERN FLOAT Gamma;
EXTERN FLOAT kcutoff;
EXTERN int   RandomSeed;
EXTERN FLOAT kmin;
EXTERN FLOAT kmax;
EXTERN int   NumberOfkPoints;
EXTERN int   PowerSpectrumType;

EXTERN char *PSFileName[2];

EXTERN float WDMPartMass;  // WDM particle mass

EXTERN float WDMg_x;  // WDM degrees of freedom

/* Space for pre-computed look-up table. */

EXTERN FLOAT *PSLookUpTable[MAX_SPECIES];

/* Internal parameters. */

EXTERN FLOAT Normalization;
EXTERN FLOAT GrowthFactor;
EXTERN FLOAT Redshift;
EXTERN FLOAT TophatRadius;
