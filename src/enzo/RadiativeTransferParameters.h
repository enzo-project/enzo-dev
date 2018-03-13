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
   min(dt) = dx_min / v_limit on the specified level */

EXTERN float RadiativeTransferTimestepVelocityLimit;
EXTERN int RadiativeTransferTimestepVelocityLevel;

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

/* The maximum distance a ray can travel. By default this will 
 * be set to sqrt(3.0) meaning that the ray can travel at most from one 
 * corner of the box to another. To avoid rays hitting the same cell 
 * twice this should be set to be <= 0.5.
 */
EXTERN float RadiativeTransferRayMaximumLength;

/* 
 * If set to TRUE use the shielding fitting function to 
 * calculate the H2 dissociation due to radiation in the 
 * Lyman-Werner band. If set to false the dissociation
 * rate is based purely on the H2I cross section for
 * radiation in the Lyman-Werner band. 
 * Default = False
 */
EXTERN int RadiativeTransferUseH2Shielding;

/* This flag controls which shielding function is used in the 
 * calculation of the H2 shielding. Current options are the 
 * standard equation due to Draine & Beltoldi (1996) or the 
 * updated function due to Wolcott-Green (2011).
 */
EXTERN int RadiativeTransferH2Shield;

/* 
 * This flag controls whether the photo-dissociation of 
 * the H2II ion is included as part of the IR and LW 
 * radiation. By default it is included.
 */
EXTERN int RadiativeTransferH2IIDiss;
/* Parameter to control flux threshold when rays are deleted when the
   photo-ionization rate becomes as long as the Hubble time times this
   constant. */

EXTERN float RadiativeTransferHubbleTimeFraction;
