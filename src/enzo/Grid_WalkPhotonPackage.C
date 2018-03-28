#define DEBUG 0
#define MYPROC MyProcessorNumber == ProcessorNumber
/***********************************************************************
/
/  GRID CLASS (WALK PHOTON PACKAGES ACROSS GRID)
/
/  written by: Tom Abel
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: This is the heart of the radiative transfer algorithm.
/    All the work is done here. Trace particles, split them, compute
/    photo and heating rates, communicate PhotonPackages if they are on
/    the same processor
/
/  RETURNS: FAIL or SUCCESS
/
/
/ N.B.
/ In the following file the convention for radiation is as follows:
/ 0 = HI ionising radiation
/ 1 = HeI ionising radiation
/ 2 = HeII ionising radiation
/ 3 = H2I ionising radiation (LW)
/ 4 = HM photodetachment radiation (IR)
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"
#include "RadiativeTransferHealpixRoutines64.h"
#define MAX_HEALPIX_LEVEL 29
#define MAX_COLUMN_DENSITY 1e25
#define MIN_TAU_IFRONT 0.1
#define TAU_DELETE_PHOTON 10.0
#define GEO_CORRECTION
#define H_SPECIES           3         //Includes HI, HeI, HeII
#define ALLSPECIES          6         //includes HI, HeI, HeII, H2I, H2II and HM
#define HIField             0
#define HeIField            1
#define HeIIField           2
#define H2IField            3
#define HIIField            4
#define HMField             5
#define H2IIField           6

int SplitPhotonPackage(PhotonPackageEntry *PP);
FLOAT FindCrossSection(int type, float energy);
float ReturnValuesFromSpectrumTable(float ColumnDensity, float dColumnDensity, int mode);
static void CalculateCrossSection(PhotonPackageEntry **PP, FLOAT *sigma, float LengthUnits, 
				  float &nSecondaryHII, float &nSecondaryHeII);
static void ResetdPi(FLOAT *dPi);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::WalkPhotonPackage(PhotonPackageEntry **PP, 
			    grid **MoveToGrid, grid *ParentGrid, grid *CurrentGrid, 
			    grid **Grids0, int nGrids0, int &DeleteMe, 
			    int &PauseMe, int &DeltaLevel, float LightCrossingTime,
			    float LightSpeed, int level, float MinimumPhotonFlux) {

  const float c_cgs = 2.99792e10;
  const float EnergyThresholds[] = {13.6, 24.6, 54.4, 11.2, 0.755, 100.0};
  const float PopulationFractions[] = {1.0, 0.25, 0.25, 1.0, 1.0, 1.0, 1.0}; //Matches Fields
  const float EscapeRadiusFractions[] = {0.5, 1.0, 2.0};
  const int offset[] = {1, GridDimension[0], GridDimension[0]*GridDimension[1]};

  int i, index, dim, splitMe, direction;
  int keep_walking, count, H2Thin, TemperatureField; 
  int type = (*PP)->Type;
  int g[3], celli[3], u_dir[3], u_sign[3];
  int cindex;
  float m[3], slice_factor, slice_factor2, sangle_inv;
  float MinTauIfront, PhotonEscapeRadius[3], c, c_inv, tau, taua, adj_taua;
  float DomainWidth[3], dx, dx2, dxhalf, fraction, dColumnDensity;
  float shield1, shield2, solid_angle, midpoint, nearest_edge;
  float tau_delete, flux_floor;
  double dN;
  FLOAT radius, oldr, cdt, dr;
  FLOAT CellVolume = 1, Volume_inv, Area_inv, SplitCriteron, SplitWithinRadius;
  FLOAT SplitCriteronIonized, PauseRadius, r_merge, d_ss, d2_ss, u_dot_d, sqrt_term;
  FLOAT dir_vec[3], sigma[H2II + 1]; //Accounts for all of the cross sections needed
  FLOAT ddr, dP, dP1, dp2,EndTime;
  FLOAT dPi[H_SPECIES + 1], dPXray[H_SPECIES + 1];  //+ 1 is to account for Compton
  FLOAT thisDensity, min_dr, adj_thisDensity;
  FLOAT ce[3], nce[3];
  FLOAT s[3], u[3], f[3], u_inv[3], r[3], dri[3];
  static int secondary_flag = 1, compton_flag = 1;
  static int photoncounter = 0;
#ifdef GEO_CORRECTION
  float absm[3], one_minus_s2;
  int closest_m, adj_index;
  FLOAT adj_dPi[H_SPECIES + 1], dP2;
#endif
  /* Check for early termination */

  if ((*PP)->Photons <= 0) {
    DeleteMe = TRUE;
    return SUCCESS;
  }

  if ((*PP) == NULL || (*PP)->PreviousPackage == NULL ||
      (*PP)->PreviousPackage->NextPackage != (*PP)) {
    ENZO_VFAIL("Called grid::WalkPhotonPackage with an invalid pointer.\n"
	    "\t %p %p %p\n",
	    (*PP), (*PP)->PreviousPackage, PhotonPackages)
  }
  /* Get units. */
  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);
  /* This controls the splitting condition, where this many rays must
     exist in each cell */

  float RaysPerCell = RadiativeTransferRaysPerCell;
  float ConvertToProperNumberDensity = DensityUnits/1.673e-24f;
  // Only split photons within this radius if specified
  SplitWithinRadius = (RadiativeTransferSplitPhotonRadius > 0) ?
    RadiativeTransferSplitPhotonRadius * (3.086e21 / LengthUnits) : 2.0;

  /* Convert escape fraction radius into code units */

  for (i = 0; i < 3; i++)
    PhotonEscapeRadius[i] = EscapeRadiusFractions[i] * 
      RadiativeTransferPhotonEscapeRadius * (3.086e21f / LengthUnits);

  // speed of light in code units. note this one is independent of a(t)
  c = LightSpeed;
  c_inv = 1.0 / LightSpeed;

  /* Calculate the normal direction (HEALPix) */
  pix2vec_nest64((int64_t) (1 << (*PP)->level), (*PP)->ipix, dir_vec);

  if (DEBUG) 
    fprintf(stderr,"grid::WalkPhotonPackage: %"GSYM" %"GSYM" %"GSYM". \n", 
	    dir_vec[0], dir_vec[1], dir_vec[2]);

  /***********************************************************************/
  /* Quantities that help finding the cell index and ray characteristics */
  /***********************************************************************/

  for (dim = 0; dim < GridRank; dim++) {
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    CellVolume *= CellWidth[dim][0];
    s[dim] = (*PP)->SourcePosition[dim];
    u[dim] = dir_vec[dim];
    u_sign[dim] = sign(u[dim]);
    u_dir[dim] = (u_sign[dim]+1) / 2;  // Ray direction

    // Zeros in y&z directions possible
    if (dim > 0)
      if (fabs(u[dim]) < PFLOAT_EPSILON)
	u[dim] = u_sign[dim]*PFLOAT_EPSILON;

    u_inv[dim] = 1.0 / u[dim];

    r[dim] = s[dim] + (*PP)->Radius * u[dim]; // Ray position

    // Current cell in integer and floating point
    g[dim] = GridStartIndex[dim] + 
      nint(floor((r[dim] - GridLeftEdge[dim]) / CellWidth[dim][0]));

    /* Check whether photon is inside the ghost zones.  They can exist
       outside after splitting in tightly packed grids, close to the
       source.  Automatically move to the parent grid. */

    if (g[dim] < 0 || g[dim] >= GridDimension[dim]) {
      DeltaLevel = -1;
      *MoveToGrid = ParentGrid;
      return SUCCESS;
    }

    f[dim] = CellLeftEdge[dim][g[dim]];

    // On cell boundaries, the index will change in negative directions
    if (r[dim] == f[dim])
      g[dim] += (u_sign[dim]-1)/2;

  }

  /* Before we do anything else, check if the photon belongs in this
     grid, or needs to be moved */

  cindex = GRIDINDEX_NOGHOST(g[0],g[1],g[2]);
  if (SubgridMarker[cindex] != this) {
    FindPhotonNewGrid(cindex, r, u, g, *PP, *MoveToGrid,
		      DeltaLevel, DomainWidth, DeleteMe, 
		      ParentGrid);
    return SUCCESS;
  }

  /* Compute the photon distance that corresponds to a distance =
     R_merge away from the super source. 
     | (PauseRadius * dir_vec) + (vector_between_source_and_super) | = R_merge
     solve for PauseRadius, which results in a quadratic equation.

     PauseRadius = -u1 +/- sqrt(u1^2 - d^2 + R_merge^2), where
     u1 = dir_vec \dot d
     d := vector between current and super source
  */

  if (RadiativeTransferSourceClustering && (*PP)->CurrentSource != NULL) {
    r_merge = RadiativeTransferPhotonMergeRadius *
      (*PP)->CurrentSource->ClusteringRadius;
    d2_ss = 0.0;
    u_dot_d = 0.0;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      d_ss = (*PP)->SourcePosition[dim] - (*PP)->CurrentSource->Position[dim];
      d2_ss += d_ss * d_ss;
      u_dot_d += dir_vec[dim] * d_ss;
    }
    sqrt_term = sqrt(u_dot_d*u_dot_d - d2_ss + r_merge*r_merge);
    if (sqrt_term > u_dot_d)
      PauseRadius = -u_dot_d + sqrt_term;
    else
      PauseRadius = -u_dot_d - sqrt_term;
    if (PauseRadius < 0)
      PauseRadius = huge_number;
  } else
    PauseRadius = huge_number;

  /* find relevant cross-section and number of secondary ionizations
     for X-rays */

  float nSecondaryHII = 1, nSecondaryHeII = 1;
  float xx, heat_factor = 1.0;
  float ion2_factor[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
 
  CalculateCrossSection(PP, sigma, LengthUnits, nSecondaryHII, nSecondaryHeII);  //sigma in cm^2 * LengthUnits
  
  MinTauIfront = MIN_TAU_IFRONT / (sigma[type] / LengthUnits);  // absorb sigma
  tau_delete = TAU_DELETE_PHOTON / (sigma[type] / LengthUnits); // and make dimensions cm^-2

  // solid angle associated with package (= 4 Pi/N_package[on this level])
  uint64_t Hlevel = (*PP)->level, res = 2L, BRP = 12L;
  uint64_t Nlevel = (BRP * (res << (res*Hlevel-1)));
  double n_on_this_level = Nlevel;
  double omega_package=4*M_PI/(n_on_this_level);
  double dtheta = sqrt(omega_package);
  if(n_on_this_level <= 0.0  || omega_package == INFINITY)  {
    fprintf(stdout, "%s: level = %llu\n", __FUNCTION__, level);
    fprintf(stdout, "%s: 2*(*PP)->level-1 = %llu\n", __FUNCTION__, res*(*PP)->level-1);
    fprintf(stdout, "%s: 12 * (2 << (2*(*PP)->level-1)) = %llu\n", __FUNCTION__, 
	    BRP * (res << (res*(*PP)->level-1)));
    fprintf(stdout, "%s: Hlevel = %llu\n", __FUNCTION__, Hlevel);
    fprintf(stdout, "%s: Nlevel = %llu\n", __FUNCTION__, Nlevel);
    fprintf(stdout, "%s: n_on_this_level = %lf\n", __FUNCTION__, n_on_this_level);
    ENZO_VFAIL("Computation of the number of pixels failed.\n"
	       "\t level = %llu\n", (*PP)->level);
  }
 
  if (RadiativeTransferAdaptiveTimestep)
    EndTime = PhotonTime + LightCrossingTime;
  else
    EndTime = PhotonTime + dtPhoton;
  
  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
  /* Find the species fields */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;
  IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
			HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  
  /* Find radiative transfer fields. */
  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				  kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum);
  const int kphNum[] = {kphHINum, kphHeINum, kphHeIINum};  //MultiSpecies = 1
  /* Find the pressure fields */
  int RPresNum1, RPresNum2, RPresNum3;
  if (RadiationPressure)
    IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3);
  /* Get the correct baryon fields (make it pretty) */

  type = (*PP)->Type;
  float *density = BaryonField[DensNum];
  float *fields[7] = { BaryonField[HINum],                                  //HI   = 0
		       BaryonField[HeINum],                                 //HeI  = 1
		       BaryonField[HeIINum],                                //HeII = 2
		       (MultiSpecies > 1) ? BaryonField[H2INum]  : NULL,    //H2I  = 3
		       BaryonField[HIINum],                                 //HII  = 4
		       (MultiSpecies > 1) ? BaryonField[HMNum]   : NULL,    //HM   = 5
		       (MultiSpecies > 1) ? BaryonField[H2IINum] : NULL};   //H2II = 6

  /* Pre-compute some quantities to speed things up */

  /* Pre-compute some quantities to speed things up */
  if(RadiativeTransferH2Shield == 1) {
    TemperatureField = this->GetTemperatureFieldNumberForH2Shield();
  }

  dx = CellWidth[0][0];
  dx2 = dx*dx;
  dxhalf = 0.5f * dx;
  SplitCriteron = dx2 / RaysPerCell;
  SplitCriteronIonized = dx2;
  Volume_inv = 1.0 / CellVolume;
  Area_inv = 1.0 / dx2;

  FLOAT emission_dt_inv = 1.0 / (*PP)->EmissionTimeInterval;
  FLOAT factor1 = emission_dt_inv;
  FLOAT ExcessEnergyfactor[5]; //Includes iHI, iHeI, iHeII, LW & IR
   
  /* For X-ray photons, we do heating and ionization for HI/HeI/HeII
     in one shot; see Table 2 of Shull & van Steenberg (1985) */
  if ((*PP)->Type == XRAYS) {
    for (i = 0; i < 3; i++) //Loop over IONISING_HI -> IONISING_HeII 
      if (RadiationXRaySecondaryIon)
	ExcessEnergyfactor[i] = factor1 * (*PP)->Energy;
      else
	ExcessEnergyfactor[i] = factor1 * ((*PP)->Energy - EnergyThresholds[i]);
    if (RadiationXRayComptonHeating) 
      TemperatureField = this->GetTemperatureFieldNumberForComptonHeating();
  }
  else if ((*PP)->Type != TRACINGSPECTRUM) {
    for (i = 0; i <= (*PP)->Type; i++) {
      if((*PP)->Type == LW)
	ExcessEnergyfactor[LW] = 0.0;
      else
	ExcessEnergyfactor[i] = factor1 * ((*PP)->Energy - EnergyThresholds[i]);
    }
  }

  /* Calculate minimum photo-ionization rate (*dOmega) before ray
     termination with a radiation background.  Probably only accurate
     with ray merging. */

  if (RadiationFieldType > 0 && RadiativeTransferSourceClustering > 0) {
    flux_floor = dtPhoton * RadiativeTransferFluxBackgroundLimit / sigma[type];
    switch (type) {
    case iHI:
      flux_floor *= RateData.k24;
      break;
    case iHeI:
      flux_floor *= RateData.k25;
      break;
    case iHeII:
      flux_floor *= RateData.k26;
      break;
    default:
      flux_floor = 0;
    } // ENDSWITCH
  } // ENDIF

  /* Calculate conversion factor for radiation pressure.  In cgs, we
     actually calculate acceleration due to radiation pressure, (all
     unit conversions are in []),

     dA = dMomentum * emission_dt_inv / Mass 
        = (N_absorbed * Energy / c) * emission_dt_inv / (CellVolume * Density) * r_hat
     ---  N_absorbed = dP * [L^3] / [t]
     dA = (dP * [L^3] / [t] * Energy / c) * emission_dt_inv / (CellVolume * Density) * r_hat

     LHS = [~] * [v] / [t]
     RHS = [L^3] / [t] * (erg_eV) / (c_cgs) * emission_dt_inv / [L^3] / [rho]

     (unit conversion for acceleration will be [~])
     [~] = RHS / ([v] / [t])
         = (erg_eV / c_cgs) * emission_dt_inv / ([t] * [rho]) / ([v] / [t])
         = (erg_eV / c_cgs) * emission_dt_inv / [rho] / [v]

     Therefore,
     dA = [~] * dP * Energy / (CellVolume * Density) * r_hat

     Since CellVolume is constant for the grid, incorporate it into [~].
     dA = [~] * dP * Energy / Density * r_hat
  */

  double RadiationPressureConversion = 0.0;
  if (RadiationPressure)
    RadiationPressureConversion =
      erg_eV * emission_dt_inv * Volume_inv / (DensityUnits * VelocityUnits * clight);

  // Mark that this grid has radiation (mainly for the coupled rate solver)
  HasRadiation = TRUE;

  DeltaLevel = 0;

  /************************************************************************/
  /*                       MAIN RAY TRACING LOOP                          */
  /************************************************************************/

  count = 0;
  keep_walking = 1;
  //cindex = GRIDINDEX_NOGHOST(g[0],g[1],g[2]);
  while (keep_walking) {

    /* If the photon has left the grid, determine MoveToGrid,
       DeltaLevel, and DeleteMe, and exit the loop. */

    if (SubgridMarker[cindex] != this) {
      FindPhotonNewGrid(cindex, r, u, g, *PP, *MoveToGrid,
			DeltaLevel, DomainWidth, DeleteMe, 
			ParentGrid);
      break;
    }

    /* Check for photons that have left the domain (only for root
       grids without periodic radiation boundaries).  We also adjust
       the photon coordinates if it needs wrapping around the periodic
       boundary.  Only for root grids that are not split because the
       ray will just be wrapped around the grid, and not moved to
       another grid.  */

    if (GravityBoundaryType != SubGridIsolated && nGrids0 == 1)
      if (this->PhotonPeriodicBoundary(cindex, r, g, s, *PP, *MoveToGrid,
    				       DomainWidth, DeleteMe) == FALSE)
    	break;

    oldr = (*PP)->Radius;
    min_dr = 1e20;
       
    /****** next cell edge crossing radii ******/

    // This and the next cell edge
    for (dim = 0; dim < 3; dim++) {
      ce[dim] = CellLeftEdge[dim][g[dim]];
      nce[dim] = CellLeftEdge[dim][g[dim] + u_dir[dim]];
    }

    // Radius of the next edge crossing in each dimension
    for (dim = 0; dim < 3; dim++)
      dri[dim] = u_inv[dim] * (nce[dim] - s[dim]);

    // the closest one is the one we want
    for (dim = 1, direction = 0, min_dr = dri[0]; dim < 3; dim++)
      if (dri[dim] < min_dr) {
	direction = dim;
	min_dr = dri[dim];
      }

    radius = min_dr + PFLOAT_EPSILON;
    dr = radius - oldr;

    if (dr < 0) {
      printf("dr < 0:   %"GSYM" %"GSYM" %"GSYM"\n", dr, min_dr, oldr);
      (*PP)->Photons = -1;
      DeleteMe = TRUE;
      return SUCCESS;
    }

    // My Position in coordinates [0..1]
    for (dim = 0; dim < 3; dim++)
      r[dim] = s[dim] + radius*u[dim];

    // splitting condition split package if its associated area is
    // larger than dx^2/RaysPerCell but only as long as the radius is
    // less than one box length.  Don't split beyond SplitWithinRadius
    // if the radiation is optically thin (Xray, LW)

    // One ray per cell in ionized cells
//    if (HII[cindex] > 0.5*fH*density[cindex])
//      splitMe = omega_package*r*r > SplitCriteronIonized;
//    else
    //float total_r = r + (*PP)->SourcePositionDiff;
    solid_angle = radius * radius * omega_package;
    splitMe = (solid_angle > SplitCriteron);

    if (splitMe && radius < SplitWithinRadius && 
	(*PP)->level < MAX_HEALPIX_LEVEL) {

      // split the package
      int return_value = SplitPhotonPackage((*PP));

      // discontinue parent ray 
      (*PP)->Photons = -1;

      DeleteMe = TRUE;
      NumberOfPhotonPackages += 4;
      return return_value;

    }  // if (splitting condition)
    
    if (DEBUG > 1) 
      fprintf(stdout, "%x %"ISYM" %"ISYM" %"ISYM" %"GSYM" %"GSYM"\t|\n",
	      (*PP), g[0], g[1], g[2], (*PP)->Radius, dr);

    index = cindex;
    ddr    = dr;

    // nor do we want transport longer than the grid timestep
    ddr    = min(ddr, c*(EndTime-(*PP)->CurrentTime));
    cdt = ddr * c_inv;

    // Check for ray merging, only consider a fraction of the ray to
    // make r=PauseRadius and return.
    if ((*PP)->Radius+ddr > PauseRadius) {
      fraction = (PauseRadius-(*PP)->Radius) / ddr;
      fraction = max(fraction, PFLOAT_EPSILON);
      //fraction = min(fraction,0.1);
      //fraction = 1.0;
      ddr *= fraction;
      cdt *= fraction;
      if (DEBUG > 1) 
	fprintf(stderr, "PAUSE: PP->Photons: %"GSYM"  PP->Radius: %"GSYM"\n",
		(*PP)->Photons, (*PP)->Radius);
      PauseMe = TRUE;
    }

    /* If requested, keep track of hydrogen ionizing photons passing
       certain radii (for photon escape fractions) */

    if (RadiativeTransferPhotonEscapeRadius > 0 && (*PP)->Type == iHI) {
      for (i = 0; i < 3; i++) {
	if (radius > PhotonEscapeRadius[i] && oldr < PhotonEscapeRadius[i])
	  EscapedPhotonCount[i+1] += (*PP)->Photons;
      } // ENDFOR i
    } // ENDIF PhotonEscapeRadius > 0

    /* Geometric correction factor because the ray's solid angle could
       not completely cover the cell */

#ifdef GEO_CORRECTION
    midpoint = oldr + 0.5f*ddr - PFLOAT_EPSILON;
    for (dim = 0; dim < 3; dim++)
      m[dim] = fabs(s[dim] + midpoint * u[dim] - (ce[dim] + dxhalf));
    for (dim = 1, nearest_edge = m[0]; dim < 3; dim++)
      if (m[dim] > nearest_edge) nearest_edge = m[dim];
    sangle_inv = 1.0 / (dtheta*radius);
    slice_factor = min(0.5f + (dxhalf-nearest_edge) * sangle_inv, 1.0f);
    slice_factor2 = slice_factor * slice_factor;
#else
    slice_factor2 = 1.0;
#endif

    // Adjust length and energy due to cosmological expansion
    // assumes that da/a=dl/l=2/3 dt/t which is strictly only true for OmegaM=1
    // note that we only uses this locally though. 

    /*
    if (ComovingCoordinates) {
      FLOAT daovera = 2/3*cdt/(*PP)->CurrentTime ;
      LengthUnits  = (LengthUnits + (1+daovera) * LengthUnits)/2;
      (*PP)->Energy  -= daovera * (*PP)->Energy;
    }
      
    */

    /* Calculate the absorption.  There are three different cases:
       single ionizations, H2 dissociation, and X-rays. */
/* Calculate the absorption. The different cases (currently) are:
     *
     * 0/1/2. Single ionizations (HI, HeI and HeII)
     * 3. H2 dissociation in the Lyman Werner band
     * 4. Infrared Radiation
     * 5. X-rays. 
     * 6. A Tracing Spectrum
     * 7. A Full Spectrum including IR, LW, UV Ionising radiation and XRAYS
     *
     */
#if DEBUG
    printf("%s: Radiation Type = %d\t PhotonEnergy = %f\n", __FUNCTION__, type, (*PP)->Energy);
#endif
    switch (type) {

      /************************************************************/
      /* HI or HeI or HeII */
      /************************************************************/
    case iHI:
    case iHeI:
    case iHeII:
      dP = dN = tau = taua = 0.0;
      ResetdPi(dPi);
      /* Loop over absorbers - type = 0 or 1 or 2 */
      for (i = 0; i <= type; i++) {
	thisDensity = PopulationFractions[i] * fields[i][index] * 
	  ConvertToProperNumberDensity; //[cm^-3] for species i
	taua = thisDensity * ddr * sigma[i];  //in cgs

#ifdef GEO_CORRECTION
	adj_thisDensity = PopulationFractions[i] * fields[i][adj_index] * 
	  ConvertToProperNumberDensity;
	adj_taua = adj_thisDensity * ddr * sigma[i];
#endif


      if(FAIL == RadiativeTransferIonization(PP, dPi, index, i, taua, factor1, 
					     ExcessEnergyfactor, slice_factor2, kphNum, 
					     gammaNum))
	{
	  fprintf(stderr, "Failed to calculate the ionizing radiation");
	  return FAIL;
	}
      // Exit the loop if everything's been absorbed
      if (taua > 20.0) break;
	
    } 

    for (i = 0; i <= type; i++) dP += dPi[i];
    (*PP)->ColumnDensity += thisDensity * ddr * LengthUnits; //in cgs

    break;

      /************************************************************/
      /* Lyman-Werner radiation  - Type 3 */
      /************************************************************/

    case LW:
      if(MultiSpecies < 2)
	break;
      dP = 0;
      thisDensity = PopulationFractions[H2IField] * fields[H2IField][index] * 
	  ConvertToProperNumberDensity;
      //Use either shielding due to Draine & Bertoldi (1996) or Wolcott-Green et al. (2011)
      if(RadiativeTransferUseH2Shielding) { 
	if(FAIL == RadiativeTransferLWShielding(PP, dP, thisDensity, ddr, index, 
						LengthUnits, kdissH2INum, TemperatureField,
						slice_factor2)) {
	  fprintf(stderr, "Failed to calculate the LW radiation");
	  return FAIL;
	}
      }
      else { //Calculate the dissociation rate based on the cross section
	dN = thisDensity * ddr;
	tau = dN*sigma[LW];  //[dimensionless]

	if(FAIL == RadiativeTransferLW(PP, dP, index, tau, factor1, 
				       slice_factor2, kdissH2INum)) {
	  fprintf(stderr, "Failed to calculate the LW radiation");
	  return FAIL;
	}
	
	(*PP)->ColumnDensity += dN * LengthUnits; //update Column Density

      }



     if(RadiativeTransferH2IIDiss == TRUE) {
	/* 
	 * A further consideration in this waveband is the dissociation of H2II.
	 */
	thisDensity = PopulationFractions[H2IIField] * fields[H2IIField][index] * 
	  ConvertToProperNumberDensity;
	dN = thisDensity * ddr;

	tau = dN * sigma[H2II];  //[dimensionless]
	
	if(FAIL == RadiativeTransferH2II(PP, index, tau, factor1, 
					 slice_factor2, kdissH2IINum)) {
	  fprintf(stderr, "Failed to calculate the LW radiation");
	  return FAIL;
	}
      }
   
      break;
      
      /************************************************************/
      /* IR radiation  Type 4 */
      /************************************************************/
    case IR:
      if(MultiSpecies < 2)
	break;
      thisDensity = PopulationFractions[HMField] * fields[HMField][index] * 
	ConvertToProperNumberDensity;
      dN = thisDensity * ddr;
      tau = dN*sigma[IR];  //[dimensionless]

     
      if(FAIL == RadiativeTransferIR(PP,dP, index, tau, factor1, 
				     ExcessEnergyfactor, slice_factor2, 
				     kphHMNum, gammaNum)) {
	fprintf(stderr, "Failed to calculate the IR radiation");
	return FAIL;
      }
      
      if(RadiativeTransferH2IIDiss == TRUE) {
	/* 
	 * A further consideration in this waveband is the dissociation of H2II.
	 */
	thisDensity = PopulationFractions[H2IIField] * fields[H2IIField][index] * 
	  ConvertToProperNumberDensity;
	dN = thisDensity * ddr;

	tau = dN * sigma[H2II];  //[dimensionless]
	if(FAIL == RadiativeTransferH2II(PP, index, tau, factor1, 
					 slice_factor2, kdissH2IINum)) {
	  fprintf(stderr, "Failed to calculate the IR radiation");
	  return FAIL;
	}
      }
      break;
	
      /************************************************************/
      /* X-rays (HI/HeI/HeII all in one!) - Type 5*/
      /************************************************************/
    case XRAYS:

      // Shull & van Steenberg (1985)
      if (RadiationXRaySecondaryIon) {
	xx = max(fields[HIIField][index] / 
		 (fields[HIIField][index] + fields[HIIField][index]), 1e-4);
	heat_factor    = 0.9971 * (1 - powf(1 - powf(xx, 0.2663f), 1.3163));

	ion2_factor[HIField] = 0.3908 * nSecondaryHII * 
	  powf(1 - powf(xx, 0.4092f), 1.7592f);
	ion2_factor[HeIField] = 0.0554 * nSecondaryHeII * 
	  powf(1 - powf(xx, 0.4614f), 1.6660f);
      }
      else {
	if(secondary_flag) {
	  fprintf(stderr, "%s: WARNING - Secondary ionisations NOT turned on\n", __FUNCTION__);
	  secondary_flag = 0;
	}
      }
      dP = 0.0; 
      ResetdPi(dPi);

      /* Loop over absorbers */
      for (i = 0; i < 3; i++) {   //##### for TraceSpectrum test 3 -> 1

	thisDensity = PopulationFractions[i] * fields[i][index] *
	  ConvertToProperNumberDensity;
	
	// optical depth of ray segment
	dN = thisDensity * ddr;
	tau = dN*sigma[i];

	if(FAIL == RadiativeTransferXRays(PP, dPi, index, i, ddr, tau, 
					  slice_factor2, factor1, ExcessEnergyfactor, 
					  ion2_factor, heat_factor, kphNum, gammaNum))
	  {
	     fprintf(stderr, "Failed to calculate the LW radiation\n");
	     return FAIL;
	  }
      } // ENDFOR absorber
      /* Update Column Densities from all species */
      (*PP)->ColumnDensity += dN * LengthUnits;

      if (RadiationXRayComptonHeating) {  

	thisDensity = BaryonField[DeNum][index] * ConvertToProperNumberDensity;
	dN = thisDensity * ddr;
	if(FAIL == RadiativeTransferComptonHeating(PP, dPi, index, LengthUnits, factor1, 
						   TemperatureField, ddr, dN, slice_factor2, 
						   gammaNum))
	  {
	     fprintf(stderr, "Failed to calculate the Compton Heating\n");
	     return FAIL;
	  }
      }
      else {
	if(compton_flag) {
	  fprintf(stderr, "%s: WARNING - Compton Heating NOT turned on\n", __FUNCTION__);
	  compton_flag = 0;
	}
      }
      
      // find the total absorbed number of photons including Compton heating
      for (i = 0; i < H_SPECIES + 1; i++) 
	dP += dPi[i] * slice_factor2;

      break;

      /************************************************************/
      /* tracing spectrum (HI/HeI/HeII all in one!) - Type 6 */
      /************************************************************/
      // instead of calculating the optical depth (tau), ColumnDensity will 
      // give all the necessary information of the photon spectrum: dP, min E
      // at the moment, secondary ionization is ignored (Alvarez & Kim 2010)
    case TRACINGSPECTRUM:

     dP = dN = 0.0;
      for (i = 0; i < 4; i++) dPXray[i] = 0.0;

      // calculate dColumnDensity of this ray segment (only for the main absorber, HI)     
      thisDensity = PopulationFractions[HIField] * fields[HIField][index] * 
	ConvertToProperNumberDensity; 
      dColumnDensity = thisDensity * ddr * LengthUnits;
      
      /* Loop over absorbers */
      for (i = 0; i < 3; i++) {   //##### for TraceSpectrum test 3 -> 1

	// the spectrum table returns the fraction of photons absorbed at this column density
	dPXray[i] = (*PP)->Photons * 
	  ReturnValuesFromSpectrumTable((*PP)->ColumnDensity, dColumnDensity, i);
	dP1 = dPXray[i] * slice_factor2;

	// units are 1/s *TimeUnits
	BaryonField[kphNum[i]][index] += dP1 * factor1; 
	
	// units are eV/s *TimeUnits;
	// the spectrum table returns the mean energy of the spectrum at this column density
	BaryonField[gammaNum][index] += dP1 * factor1 * 
	  ( ReturnValuesFromSpectrumTable((*PP)->ColumnDensity, dColumnDensity, 3) - 
	    EnergyThresholds[i] );

      }

      // update column density
      (*PP)->ColumnDensity += dColumnDensity;

      // find the total absorbed number of photons 
      for (i = 0; i < 4; i++) dP += dPXray[i] * slice_factor2;

      break;
    default:
      printf("Photon type = %d, radius = %g, pos = %"FSYM" %"FSYM" %"FSYM"\n",
	     type, radius, r[0], r[1], r[2]);
      ENZO_FAIL("Bad photon type.");

    } // ENDSWITCH type


    /* Keep track of the maximum hydrogen photo-ionization rate in the
       I-front, so we can calculate the maximum ionization timescale
       for timestepping purposes. */

    if (RadiativeTransferHIIRestrictedTimestep)
      if (type == iHI || type == XRAYS) {
	if ((*PP)->ColumnDensity > MinTauIfront) {
	  if (BaryonField[kphNum[iHI]][index] > this->MaximumkphIfront) {
	    this->MaximumkphIfront = BaryonField[kphNum[iHI]][index];
	    this->IndexOfMaximumkph = index;
	  } // ENDIF max
	} // ENDIF tau > min_tau (I-front)
      } // ENDIF type==iHI || Xrays
      
    /* Acceleration due to radiation pressure */

    // Remember:  dA = [~] * dP * Energy / Density * r_hat
    if (RadiationPressure && 
	(*PP)->Radius >= (*PP)->SourcePositionDiff)
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	BaryonField[RPresNum1+dim][index] += 
	  RadiationPressureConversion * RadiationPressureScale * dP * (*PP)->Energy / 
	  density[index] * dir_vec[dim];
    
    (*PP)->CurrentTime += cdt;
    (*PP)->Photons     -= dP;
    (*PP)->Radius      += ddr;
    
    if (RadiativeTransferLoadBalance) {
      int RaySegNum = FindField(RaySegments, FieldType, NumberOfBaryonFields);
      BaryonField[RaySegNum][index] += 1.0;
    }

    // return in case we're pausing to merge
    if (PauseMe)
      return SUCCESS;

    // return in case we're out of photons
    if ((*PP)->Photons < MinimumPhotonFlux*(solid_angle*Area_inv) || 
	(*PP)->ColumnDensity > tau_delete) {
      if (DEBUG > 1) {
	fprintf(stderr, "PP-Photons: %"GSYM" (%"GSYM"), PP->Radius: %"GSYM
		"PP->CurrentTime: %"FSYM"\n",
		(*PP)->Photons, MinimumPhotonFlux*solid_angle*Area_inv, 
		(*PP)->Radius, (*PP)->CurrentTime);
	fprintf(stderr, "\tdP: %"GSYM"\tddr: %"GSYM"\t cdt: %"GSYM"\t tau: %"GSYM"\n", 
		dP, ddr, cdt, tau);
      }
      (*PP)->Photons = -1;
      DeleteMe = TRUE;
      return SUCCESS;
    }

    // If we are using a radiation background, check if the
    // optically-thin flux has dropped below 1% of the background.
    // Only can be used accurately with ray merging.
    // Ray photon flux < dOmega * 0.01 * photo-ionization(background)
    //   dOmega = 4*pi*r^2/N_pixels
    //   flux_floor = 0.01 * photon-ionization(background) / cross-section
    if (RadiationFieldType > 0 && RadiativeTransferSourceClustering > 0)
      if ((*PP)->Photons < flux_floor * solid_angle) {
	if (DEBUG)
	  printf("Deleting photon %p: r=%g, P=%g, limit=%g\n", *PP, radius,
		 (*PP)->Photons, flux_floor*solid_angle);
	(*PP)->Photons = -1;
	DeleteMe = TRUE;
	return SUCCESS;
      }
    
    // are we done ? 
    if (((*PP)->CurrentTime) >= EndTime) {
      if (RadiativeTransferAdaptiveTimestep) {
	(*PP)->Photons = -1;
	DeleteMe = TRUE;
      }
      return SUCCESS;
    }

    count++;
    
    g[direction] += u_sign[direction];
    cindex += u_sign[direction] * offset[direction];
    
  } // while keep walking

  return SUCCESS;
}
/*
 * Set cross sections for the photo rates. 
 * The H2I cross section is taken from a fit to the 
 * Stancil et al. (1994) data. The data is a function of both 
 * temperature and frequency and so I fit the data to a 
 * temperarure close to 5000K. The fit is shown below.
 * A significantly more accurate treatment of the cross section 
 * can be obtained by using a look-up table (provided) by setting USE_H2II_LOOKUP 
 * to be true in RadiativeTransferH2IILoopUpTable.h
 * This comes with a (substantial) computational cost. Currently Disabled.
 */
static void CalculateCrossSection(PhotonPackageEntry **PP, 
				  FLOAT *sigma, float LengthUnits, 
				  float &nSecondaryHII, float &nSecondaryHeII)
{
  int i = 0; 
  float X = 0.0;
  float A = 3.486e-16;
  float a = 3.35485518, b = 0.93891875, c = -0.01176537;
  /* If Radiation is ionising radiation (type = 0 or 1 or 2)*/
  if ((*PP)->Type == iHI || (*PP)->Type == iHeI || 
      (*PP)->Type == iHeII) {
    for (i = 0; i <= (*PP)->Type; i++)
      if (i == (*PP)->Type)
	sigma[i] = (*PP)->CrossSection * LengthUnits;
      else
	sigma[i] = FindCrossSection(i, (*PP)->Energy) * LengthUnits; //fits from Verner et al. 1996
  } /* If Radiation is H2 dissociating radiation (type = 3)*/
  else if ((*PP)->Type == LW) {
    sigma[LW] = 3.71e-18 * LengthUnits; // H2I average cross-section
     /* Fit created from Stancil et al. 1994 */
    X = log((*PP)->Energy/8.0);
    sigma[H2II] = a*pow(10.0, -b*X*X)*exp(-c*X) * 1e-18 * LengthUnits;
    //sigma[H2II] = 1.495509e-18 * LengthUnits; // H2II average cross-section
    } /* If Radiation is in XRays (type = 5) plus ionising radiation */
  else if ((*PP)->Type == XRAYS) {
    for (i = 0; i < 3; i++)
      sigma[i] = FindCrossSection(i, (*PP)->Energy) * LengthUnits; //fits from Verner et al. 1996
    nSecondaryHII = (*PP)->Energy / 13.6;
    nSecondaryHeII = (*PP)->Energy / 24.6; 
  } /* If radiation is in the IR (type = 4) */
  else if (((*PP)->Type == IR)) {
    /* Fit taken from Tegmark et al. (1997) */
    X = (*PP)->Energy / 0.74;
    sigma[IR] = (A*pow((X - 1), 1.5)/pow(X, 3.11)) * LengthUnits;
    /* Fit created from Stancil et al. 1994 */
    X = log((*PP)->Energy/8.0);
    sigma[H2II] = a*pow(10.0, -b*X*X)*exp(-c*X) * 1e-18 * LengthUnits;
  }
  return;
}

static void ResetdPi(FLOAT *dPi)
{
  int i = 0;
  for (i = 0; i < H_SPECIES + 1; i++) {
    dPi[i] = 0.0;
  }
  return;
}

