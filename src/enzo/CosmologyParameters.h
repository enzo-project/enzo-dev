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

/* The value of Omega due to non-relativistic collisionless particles at z=0. */

EXTERN float OmegaDarkMatterNow;

/* The value of Omega due to lamba (the cosmological constant) at z=0. */

EXTERN float OmegaLambdaNow;

/* The value of Omega due to relativistic particles at z=0.
   rad = photon + neutrino; neutrinos are considered relativistic here.
   At z>~100, all neutrinos are close to relativistic, and later their
   impact on the Friedmann equation diminishes.                         */

EXTERN float OmegaRadiationNow;

/* The comoving size of the simulation box (along the x-dir) in h^{-1} Mpc. */

EXTERN float ComovingBoxSize;

/* The maximum allowed fractional increase of the expansion. */

EXTERN float MaxExpansionRate;

/* The time in code units at the initial redshift. */

EXTERN FLOAT InitialTimeInCodeUnits;

/* The initial redshift. */

EXTERN FLOAT InitialRedshift;

/* The final redshift. */

EXTERN FLOAT FinalRedshift;

/* Redshift output information: */

EXTERN FLOAT CosmologyOutputRedshift[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
EXTERN char *CosmologyOutputRedshiftName[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];
EXTERN FLOAT CosmologyOutputRedshiftTime[MAX_NUMBER_OF_OUTPUT_REDSHIFTS];

/* Tabulated a(t). */

EXTERN int CosmologyTableNumberOfBins;
EXTERN FLOAT CosmologyTableLogaInitial;
EXTERN FLOAT CosmologyTableLogaFinal;
EXTERN FLOAT *CosmologyTableLoga;
EXTERN FLOAT *CosmologyTableLogt;
