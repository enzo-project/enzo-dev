/* Cosmology Parameters */

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* The Hubble constant at z=0 (now) in units of 100 km/s/Mpc. */

EXTERN float HubbleConstantNow;

/* The value of Omega due to non-relativistic particles at z=0. */

EXTERN float OmegaMatterNow;

/* The value of Omega due to lamba (the cosmological constant) at z=0. */

EXTERN float OmegaLambdaNow;

/* The comoving size of the simulation box (along the x-dir) in h^{-1} Mpc. */

EXTERN float ComovingBoxSize;

/* The maximum allowed fractional increase of the expansion. */

EXTERN float MaxExpansionRate;

/* The time in code units at the initial redshift. */

EXTERN FLOAT InitialTimeInCodeUnits;

/* The initial redshift. */

EXTERN FLOAT InitialRedshift;

/* Redshift output information: */

EXTERN FLOAT CosmologyOutputRedshift[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
EXTERN char *CosmologyOutputRedshiftName[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
EXTERN FLOAT CosmologyOutputRedshiftTime[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
