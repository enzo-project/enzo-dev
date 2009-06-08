/* Cosmology Parameters */

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* The Hubble constant at z=0 (now) in units of 100 km/s/Mpc. */

EXTERN FLOAT HubbleConstantNow;

/* The value of Omega due to non-relativistic particles at z=0. */

EXTERN FLOAT OmegaMatterNow;

/* The value of Omega due to lamba (the cosmological constant) at z=0. */

EXTERN FLOAT OmegaLambdaNow;

/* The value of Omega due to WDM at z=0 */

EXTERN float OmegaWDMNow;

/* The value of Omega due to HDM (neutrinos) at z=0. */

EXTERN FLOAT OmegaHDMNow;

/* The value of Omega due to baryons at z=0. */

EXTERN FLOAT OmegaBaryonNow;

/* The comoving size of the simulation box (along the x-dir) in h^{-1} Mpc. */

EXTERN FLOAT ComovingBoxSize;

/* The initial redshift. */

EXTERN FLOAT InitialRedshift;

