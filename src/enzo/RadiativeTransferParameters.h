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

/* Radius to merge rays in units of separation of the two sources
   associated with a super source. */

EXTERN float RadiativeTransferPhotonMergeRadius;

/* Radiative pressure flag and scale factor */

EXTERN int RadiationPressure;
EXTERN float RadiationPressureScale;

/* Flag to turn on a 1/r^2 Lyman-Werner radiation field */

EXTERN int RadiativeTransferOpticallyThinH2;

/* Periodic boundary conditions for the photon packages */

EXTERN int RadiativeTransferPeriodicBoundary;

/* Which level, i.e. the frequency it's called, we call the FLD
   solver */

EXTERN int RadiativeTransferFLDCallOnLevel;

/* Flag to use timestepping to restrict HII fraction change to 50% */

EXTERN int RadiativeTransferHIIRestrictedTimestep;

/* Flag to use adaptive timestepping for ray tracing.  However this
   removes the time-derivative from the radiative transfer
   equation. */

EXTERN int RadiativeTransferAdaptiveTimestep;

EXTERN float GlobalMaximumkphIfront;

/* Angle (in degrees) to collimate radiation for PhotonSourceType = 2 */

EXTERN float RadiativeTransferSourceBeamAngle;

/* Flag to trace the spectrum in ray tracing */

EXTERN int RadiativeTransferTraceSpectrum;

EXTERN char *RadiativeTransferTraceSpectrumTable;

/* Flag for temporary load balancing for the ray tracing */

EXTERN int RadiativeTransferLoadBalance;

/* Flux threshold when rays are deleted in units of the UV background
   flux (RadiationFieldType > 0) */

EXTERN float RadiativeTransferFluxBackgroundLimit;
