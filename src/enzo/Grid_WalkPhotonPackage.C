#define DEBUG 0
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
#include "RadiativeTransferHealpixRoutines.h"

#ifdef CONFIG_BFLOAT_4
#define ROUNDOFF 1e-6f
#endif
#ifdef CONFIG_BFLOAT_8
#define ROUNDOFF 1e-12
#endif
#ifdef CONFIG_BFLOAT_16
#define ROUNDOFF 1e-16
#endif

#define MAX_HEALPIX_LEVEL 13
#define MAX_COLUMN_DENSITY 1e25
#define MIN_TAU_IFRONT 0.1

int SplitPhotonPackage(PhotonPackageEntry *PP);
FLOAT FindCrossSection(int type, float energy);

enum species { iHI, iHeI, iHeII, iH2I, iHII };
int grid::WalkPhotonPackage(PhotonPackageEntry **PP, 
			    grid **MoveToGrid, grid *ParentGrid, grid *CurrentGrid, 
			    grid **Grids0, int nGrids0, int DensNum, int DeNum,
			    int HINum, int HeINum, int HeIINum, int H2INum, 
			    int kphHINum, int gammaNum, int kphHeINum, 
			    int kphHeIINum, 
			    int kdissH2INum, int RPresNum1, int RPresNum2, 
			    int RPresNum3, int &DeleteMe, int &PauseMe, 
			    int &DeltaLevel, float LightCrossingTime,
			    float DensityUnits, float TemperatureUnits,
			    float VelocityUnits, float LengthUnits,
			    float TimeUnits) {

  const float erg_eV = 1.602176e-12;
  const float c_cgs = 2.99792e10;
  const float EnergyThresholds[] = {13.6, 24.6, 54.4, 11.2};
  const float PopulationFractions[] = {1.0, 0.25, 0.25, 1.0};
  const float EscapeRadiusFractions[] = {0.5, 1.0, 2.0};
  const int kphNum[] = {kphHINum, kphHeINum, kphHeIINum};

  float ConvertToProperNumberDensity = DensityUnits/1.673e-24f;

  int i, index, dim, splitMe, direction;
  int keep_walking, count, H2Thin, type;
  int g[3], celli[3], u_dir[3], u_sign[3];
  long cindex;
  float m[3], slice_factor, slice_factor2, sangle_inv;
  float MinTauIfront, PhotonEscapeRadius[3], c, c_inv, tau;
  float DomainWidth[3], dx, dx2, dxhalf, fraction;
  float shield1, shield2, solid_angle, midpoint, nearest_edge;
  double dN;
  FLOAT radius, oldr, cdt, dr;
  FLOAT CellVolume = 1, Volume_inv, Area_inv, SplitCriteron, SplitWithinRadius;
  FLOAT SplitCriteronIonized, PauseRadius, r_merge, d_ss, d2_ss, u_dot_d, sqrt_term;
  FLOAT dir_vec[3], sigma[4]; 
  FLOAT ddr, dP, dP1, EndTime;
  FLOAT xE, dPXray[4];  
  FLOAT thisDensity, min_dr;
  FLOAT ce[3], nce[3];
  FLOAT s[3], u[3], f[3], u_inv[3], r[3], dri[3];

  /* Check for early termination */

  if ((*PP)->Photons <= 0) {
    DeleteMe = TRUE;
    return SUCCESS;
  }

  if ((*PP) == NULL || (*PP)->PreviousPackage->NextPackage != (*PP)) {
    fprintf(stderr, "Called grid::WalkPhotonPackage with an invalid pointer.\n"
	    "\t %x %x %x %x\n",
	    (*PP), (*PP)->PreviousPackage, (*PP)->PreviousPackage->NextPackage,
	    PhotonPackages);
    ENZO_FAIL("");
  }

  /* This controls the splitting condition, where this many rays must
     exist in each cell */

  float RaysPerCell = RadiativeTransferRaysPerCell;

  // Only split photons within this radius if specified
  SplitWithinRadius = (RadiativeTransferSplitPhotonRadius > 0) ?
    RadiativeTransferSplitPhotonRadius * (3.086e21 / LengthUnits) : 2.0;

  /* Convert escape fraction radius into code units */

  for (i = 0; i < 3; i++)
    PhotonEscapeRadius[i] = EscapeRadiusFractions[i] * 
      RadiativeTransferPhotonEscapeRadius * (3.086e21f / LengthUnits);

  // speed of light in code units. note this one is independent of a(t)
  c = c_cgs/VelocityUnits;

  // Modify the photon propagation speed by this parameter
  c *= RadiativeTransferPropagationSpeedFraction;
  c_inv = 1.0 / c;

  /* Calculate the normal direction (HEALPix) */

  if (pix2vec_nest((long) (1 << (*PP)->level), (*PP)->ipix, dir_vec)==FAIL) {
    fprintf(stdout,"WalkPhotonPackage: pix2vec_nest outor %ld %ld %g %x\n",
	    (long) (1 << (*PP)->level), (*PP)->ipix, (*PP)->Photons, (*PP) );
    ENZO_FAIL("");
  }

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
      if (fabs(u[dim]) < ROUNDOFF)
	u[dim] = u_sign[dim]*ROUNDOFF;

    u_inv[dim] = 1.0 / u[dim];

    r[dim] = s[dim] + (*PP)->Radius * u[dim]; // Ray position

    // Current cell in integer and floating point
    g[dim] = GridStartIndex[dim] + 
      (int) ((r[dim] - GridLeftEdge[dim]) / CellWidth[dim][0]);
    f[dim] = CellLeftEdge[dim][g[dim]];

    // On cell boundaries, the index will change in negative directions
    if (r[dim] == f[dim])
      g[dim] += (u_sign[dim]-1)/2;

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
    r_merge = 2*RadiativeTransferPhotonMergeRadius *
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
  } else
    PauseRadius = huge_number;

  /* find relevant cross-section and number of secondary ionizations
     for X-rays */

  int nSecondaryHII = 1, nSecondaryHeIII = 1;
  float xx, heat_factor = 1.0;
  float ion2_factor[] = {1.0, 1.0, 1.0};

  if ((*PP)->Type == iHI || (*PP)->Type == iHeI || (*PP)->Type == iHeII) {
    sigma[0] = (*PP)->CrossSection * LengthUnits;
  }
  else if ((*PP)->Type == iH2I) {
    sigma[0] = 3.71e-18 * LengthUnits; // H2I average cross-section
  }
  else if ((*PP)->Type == 4) {
    for (i = 0; i < 3; i++)
      sigma[i] = FindCrossSection(i, (*PP)->Energy) * LengthUnits;
    nSecondaryHII = (int) floor((*PP)->Energy / 13.6 - 1);
    nSecondaryHeIII = (int) floor((*PP)->Energy / 54.4 - 1);
  }
  
  MinTauIfront = MIN_TAU_IFRONT / sigma[0];  // absorb sigma

  // solid angle associated with package (= 4 Pi/N_package[on this level]) 
  float n_on_this_level = (float) (12 * (2 << (2*(*PP)->level-1)));
  FLOAT omega_package=4*M_PI/(n_on_this_level);
  float dtheta = sqrt(omega_package);
 
  if (RadiativeTransferAdaptiveTimestep)
    EndTime = PhotonTime + LightCrossingTime;
  else
    EndTime = PhotonTime + dtPhoton;

  /* Get the correct baryon fields (make it pretty) */

  type = (*PP)->Type;
  float *density = BaryonField[DensNum];
  float *fields[5] = { BaryonField[HINum],
		       BaryonField[HeINum],
		       BaryonField[HeIINum],
		       (MultiSpecies > 1) ? BaryonField[H2INum] : NULL,
		       BaryonField[HINum+1] };

  /* Pre-compute some quantities to speed things up */

  dx = CellWidth[0][0];
  dx2 = dx*dx;
  dxhalf = 0.5f * dx;
  SplitCriteron = dx2 / RaysPerCell;
  SplitCriteronIonized = dx2;
  Volume_inv = 1.0 / CellVolume;
  Area_inv = 1.0 / dx2;

  FLOAT emission_dt_inv = 1.0 / (*PP)->EmissionTimeInterval;
  FLOAT factor1 = emission_dt_inv;
  FLOAT factor2[3];
  FLOAT factor3 = Area_inv*emission_dt_inv;

  /* For X-ray photons, we do heating and ionization for HI/HeI/HeII
     in one shot */

  if ((*PP)->Type == 4) 
    for (i = 0; i < 3; i++)
      factor2[i] = factor1 * ((*PP)->Energy - EnergyThresholds[i]);
  else 
    factor2[0] = factor1 * ((*PP)->Energy - EnergyThresholds[type]);


  /* Calculate conversion factor for radiation pressure.  In cgs, we
     actually calculate acceleration due to radiation pressure, (all
     unit conversions are in []),

     dA = dMomentum / Mass 
        = (N_absorbed * Energy / c) / (CellVolume * Density) * r_hat
     ---  N_absorbed = dP * [L^3] / [t]
     dA = (dP * [L^3] / [t] * Energy / c) / (CellVolume * Density) * r_hat

     LHS = [~] * [v] / [t]
     RHS = [L^3] / [t] * (erg_eV) / (c_cgs) / [L^3] / [rho]

     (unit conversion for acceleration will be [~])
     [~] = RHS / ([v] / [t])
         = (erg_eV / c_cgs) / ([t] * [rho]) / ([v] / [t])
         = (erg_eV / c_cgs) / [rho] / [v]

     Therefore,
     dA = [~] * dP * Energy / (CellVolume * Density) * r_hat

     Since CellVolume is constant for the grid, incorporate it into [~].
     dA = [~] * dP * Energy / Density * r_hat

  */

  double RadiationPressureConversion =
    erg_eV / c_cgs / DensityUnits / VelocityUnits * Volume_inv;

  // Mark that this grid has radiation (mainly for the coupled rate solver)
  HasRadiation = TRUE;

  DeltaLevel = 0;

  /************************************************************************/
  /*                       MAIN RAY TRACING LOOP                          */
  /************************************************************************/

  count = 0;
  keep_walking = 1;
  while (keep_walking) {

    /* If the photon has left the grid, determine MoveToGrid,
       DeltaLevel, and DeleteMe, and return. */

    if (PointInGridNB(r) == FALSE)
      if (FindPhotonNewGrid(Grids0, nGrids0, r, u, *PP, *MoveToGrid,
			    DeltaLevel, DomainWidth, DeleteMe, 
			    ParentGrid) == FALSE)
	return SUCCESS;

    cindex = GRIDINDEX_NOGHOST(g[0],g[1],g[2]);
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
    for (dim = 0; dim < 3; dim++)
      if (dri[dim] < min_dr) {
	direction = dim;
	min_dr = dri[dim];
      }

    radius = min_dr + ROUNDOFF;
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

    if (SubgridMarker[cindex] != CurrentGrid) {
      if (SubgridMarker[cindex] == NULL) {
	(*MoveToGrid) = ParentGrid;
	DeltaLevel = -1;
      } else {
	(*MoveToGrid) = SubgridMarker[cindex];
	DeltaLevel = 1;
	if (DEBUG) 
	  printf("different grid subgrid marker %x %x %ld %"ISYM" %"ISYM
		 " %"ISYM" %"FSYM" %"FSYM" %"FSYM"\n",
		 SubgridMarker[cindex], CurrentGrid, cindex, 
		 g[0], g[1], g[2], r[0], r[1], r[2]);
      }
      // move it at least a tiny fraction of the grid cell to not have
      // to worry about round off errors shifting photons back and
      // forth between two grids without doing anything.
      (*PP)->Radius += ROUNDOFF;
      return SUCCESS;
    }

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

    midpoint = oldr + 0.5f*ddr - ROUNDOFF;
    nearest_edge = -1e20;
    for (dim = 0; dim < 3; dim++)
      m[dim] = fabs(s[dim] + midpoint * u[dim] - (ce[dim] + dxhalf));
    nearest_edge = max(max(m[0], m[1]), m[2]);
    sangle_inv = 1.0 / (dtheta*radius);
    slice_factor = min(0.5f + (dxhalf-nearest_edge) * sangle_inv, 1.0f);
    slice_factor2 = slice_factor * slice_factor;

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

    switch (type) {

      /************************************************************/
      /* HI or HeI or HeII */
      /************************************************************/
    case iHI:
    case iHeI:
    case iHeII:
      thisDensity = PopulationFractions[type] * fields[type][index] * 
	ConvertToProperNumberDensity; 

      // optical depth of ray segment
      dN = thisDensity * ddr;
      tau = dN*sigma[0];

      // at most use all photons for photo-ionizations
      if (tau > 2.e1) dP = (1.0+ROUNDOFF) * (*PP)->Photons;
      else if (tau > 1.e-4) 
	dP = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
      else
	dP = min((*PP)->Photons*tau, (*PP)->Photons);
      dP1 = dP * slice_factor2;

      // contributions to the photoionization rate is over whole timestep
      BaryonField[kphNum[type]][index] += dP1*factor1;
	
      // the heating rate is just the number of photo ionizations
      // times the excess energy units here are eV/s *TimeUnits.
      BaryonField[gammaNum][index] += dP1*factor2[0];

      break;

      /************************************************************/
      /* Lyman-Werner radiation */
      /************************************************************/
    case iH2I:
      if (MultiSpecies > 1)
	thisDensity = PopulationFractions[type] * fields[type][index] * 
	  ConvertToProperNumberDensity;

      /* We treat H2 dissociation with the shielding function from
	 Draine & Bertoldi (1996) */

      if ((*PP)->ColumnDensity < 1e14) {
	shield1 = 1;
	H2Thin = TRUE;
      } else {
	shield1 = pow((*PP)->ColumnDensity / 1e14, -0.75);
	H2Thin = FALSE;
      }

      (*PP)->ColumnDensity += thisDensity * ddr * LengthUnits;
      if ((*PP)->ColumnDensity < 1e14) {
	shield2 = 1;
      } else {
	shield2 = pow((*PP)->ColumnDensity / 1e14, -0.75);
	H2Thin = FALSE;
      }

      if (H2Thin == TRUE)
	dP = 0;
      else
	dP = (*PP)->Photons * (1 - shield2/shield1);
      BaryonField[kdissH2INum][index] += (*PP)->Photons * sigma[0] * factor3 * 
	(ddr * omega_package * radius * radius * Volume_inv);

      break;

      /************************************************************/
      /* X-rays (HI/HeI/HeII all in one!) */
      /************************************************************/
    case 4:

      if (RadiationXRaySecondaryIon) {
	xx = max(fields[iHII][index] / 
		 (fields[iHI][index] + fields[iHII][index]), 1e-4);
	heat_factor    = 0.9971 * (1 - powf(1 - powf(xx, 0.2663f), 1.3163));
	ion2_factor[0] = 0.3908 * nSecondaryHII * 
	  powf(1 - powf(xx, 0.4092f), 1.7592f);
	ion2_factor[2] = 0.0554 * nSecondaryHeIII * 
	  powf(1 - powf(xx, 0.4614f), 1.6660f);
      }

      dP = 0.0; 
      for (i = 0; i < 4; i++) dPXray[i] = 0.0; //#####

      /* Loop over absorbers */
      for (i = 0; i < 3; i++) {

	thisDensity = PopulationFractions[i] * fields[i][index] *
	  ConvertToProperNumberDensity;
	
	// optical depth of ray segment
	dN = thisDensity * ddr;
	tau = dN*sigma[i];

	// at most use all photons for photo-ionizations
	if (tau > 2.e1) dPXray[i] = (1.0+ROUNDOFF) * (*PP)->Photons;
	else if (tau > 1.e-4) 
	  dPXray[i] = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
	else
	  dPXray[i] = min((*PP)->Photons*tau, (*PP)->Photons);
	dP1 = dPXray[i] * slice_factor2;

	// contributions to the photoionization rate is over whole timestep
	BaryonField[kphNum[i]][index] += dP1 * factor1 * ion2_factor[i];
	
	// the heating rate is just the number of photo ionizations times
	// the excess energy units here are  eV/s/cm^3 *TimeUnits.  
	BaryonField[gammaNum][index] += dP1 * factor2[i] * heat_factor;

      } // ENDFOR absorber

      // assume photon energy is much less than the electron rest mass energy 
      if (RadiationXRayComptonHeating) {  //#####

	thisDensity = BaryonField[DeNum][index] * ConvertToProperNumberDensity;

	// nonrelativistic Klein-Nishina cross section and optical depth
	// Ribicki & Lightman (1979)
	xE = (*PP)->Energy/5.11e5;  // mc^2 = 0.511 MeV
	sigma[3] = 6.65e-25 * (1 - 2.*xE + 26./5.*xE*xE) * LengthUnits;

	dN = thisDensity * ddr;
	tau = dN*sigma[3];

	// at most use all photons for Compton scattering
	if (tau > 2.e1) dPXray[3] = (1.0+ROUNDOFF) * (*PP)->Photons;
	else if (tau > 1.e-4) 
	  dPXray[3] = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
	else
	  dPXray[3] = min((*PP)->Photons*tau, (*PP)->Photons);
	dP1 = dPXray[3] * slice_factor2;

	BaryonField[gammaNum][index] += dP1 * factor1 * (*PP)->Energy * xE;

      }
      
      // find the total absorbed number of photons including Compton heating
      for (i = 0; i < 4; i++) dP += dPXray[i];

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
      if (type == iHI || type == 4) {
	(*PP)->ColumnDensity += dN;
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
	(*PP)->Radius > (*PP)->SourcePositionDiff)
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	BaryonField[RPresNum1+dim][index] += 
	  RadiationPressureConversion * dP * (*PP)->Energy / 
	  density[index] * dir_vec[dim];
    
    (*PP)->CurrentTime += cdt;
    (*PP)->Photons     -= dP;
    (*PP)->Radius      += ddr;

    //BaryonField[kphHeIINum][index] += 1;

    // return in case we're pausing to merge
    if (PauseMe)
      return SUCCESS;

    // return in case we're out of photons
    if ((*PP)->Photons < tiny_number) {
      if (DEBUG>1) 
	fprintf(stderr, "PP-Photons: %"GSYM"  PP->Radius: %"GSYM
		"PP->CurrentTime: %"FSYM"\n",
		(*PP)->Photons, (*PP)->Radius, (*PP)->CurrentTime);
      if (DEBUG>1) 
	fprintf(stderr, "\tdP: %"GSYM"\tddr: %"GSYM"\t cdt: %"GSYM"\t tau: %"GSYM"\n", 
		dP, ddr, cdt, tau);
      (*PP)->Photons = -1;
      DeleteMe = TRUE;
      return SUCCESS;
    }
    
    // are we done ? 
    if (((*PP)->CurrentTime) >= EndTime) {
      (*PP)->Photons = -1;
      DeleteMe = TRUE;
      return SUCCESS;
    }

    count++;
    
    g[direction] += u_sign[direction];
    
  } // while keep walking

  return SUCCESS;
}
