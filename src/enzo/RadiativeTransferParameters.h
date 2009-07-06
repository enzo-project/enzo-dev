/* Radiative Transfer Parameters */

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* Radiation emitted not from a single point but from the surface of a sphere */
/* Radius of that sphere in cell-widths of the finest grid on which the source is found  */
/* Radius = 0 : pointsource */

EXTERN float RadiativeTransferSourceRadius;

/* Fraction of speed of Light to evolve photons on */

EXTERN float RadiativeTransferPropagationSpeedFraction;

/* Fraction of box length packages move in one time step */

EXTERN FLOAT RadiativeTransferPropagationDistance;

/* Only split photon packages within this radius (kpc) */

EXTERN float RadiativeTransferSplitPhotonRadius;

/* Criteron for splitting (rays per cell) */

EXTERN float RadiativeTransferRaysPerCell;

/* Initial HEALPix level from point source */

EXTERN int RadiativeTransferInitialHEALPixLevel;

/* Base radius from which to measure photon escape fractions (kpc) */

EXTERN float RadiativeTransferPhotonEscapeRadius;

/* Whether to use the interpolated densities (i.e. ray segment
   centers) for calculating the optical depth */

EXTERN int RadiativeTransferInterpolateField;

/* Associated velocity in km/s to limit the photon timestep.  
   min(dt) = dx_min / v_limit */

EXTERN float RadiativeTransferTimestepVelocityLimit;

/* Flag whether to cluster sources on a binary tree and merge rays at
   a certain radii with associated SuperSources */

EXTERN int RadiativeTransferSourceClustering;

EXTERN float RadiativeTransferPhotonMergeRadius;

EXTERN int RadiationPressure;

EXTERN int RadiativeTransferOpticallyThinH2;

EXTERN int RadiativeTransferPeriodicBoundary;

EXTERN int RadiativeTransferFLD;
