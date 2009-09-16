#define DEBUG 1 //#####
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
#define ROUNDOFF 1e-6
#endif
#ifdef CONFIG_BFLOAT_8
#define ROUNDOFF 1e-12
#endif
#ifdef CONFIG_BFLOAT_16
#define ROUNDOFF 1e-16
#endif

#define MAX_HEALPIX_LEVEL 13
#define MAX_COLUMN_DENSITY 1e25

int FindRootGrid(int &dummy, grid **Grids0, int nGrids0, 
		 FLOAT rx, FLOAT ry, FLOAT rz, FLOAT ux, FLOAT uy, FLOAT uz);
int SplitPhotonPackage(PhotonPackageEntry *PP);
FLOAT FindCrossSection(int type, float energy);
float ComputeInterpolatedValue(float *vc, int vci, int vcj, int vck, 
			       float mx, float my, float mz);

int grid::WalkPhotonPackage(PhotonPackageEntry **PP, 
			    grid **MoveToGrid, grid *ParentGrid, grid *CurrentGrid, 
			    grid **Grids0, int nGrids0, int DensNum, 
			    int HINum, int HeINum, int HeIINum, int H2INum, 
			    int kphHINum, int gammaHINum, int kphHeINum, 
			    int gammaHeINum, int kphHeIINum, int gammaHeIINum, 
			    int kdissH2INum, int RPresNum1, int RPresNum2, 
			    int RPresNum3, int &DeleteMe, int &PauseMe, 
			    int &DeltaLevel,
			    float DensityUnits, float TemperatureUnits,
			    float VelocityUnits, float LengthUnits,
			    float TimeUnits) {

  float mx, my, mz, slice_factor, slice_factor2;
  FLOAT rx, ry, rz, r, oldr, drx, dry, drz, prev_radius;
  FLOAT CellVolume=1, Volume_inv, Area_inv, SplitCriteron, SplitWithinRadius;
  FLOAT SplitCriteronIonized, PauseRadius, r_merge, d_ss, d2_ss, u_dot_d, sqrt_term;

  int i, j, k, index, dim, splitMe;
  int kphNum, gammaNum;
  i=j=k=index=0;
#define HIRecombinationRate 3.6e-13

  // absorbs the LengthUnits for tau computation here
//  FLOAT        sigmaHI   = 6.3e-18 * pow(((*PP)->Energy/13.6), -3) * LengthUnits;
//  FLOAT	sigmaHeI  = 6.98e-18 * LengthUnits;
//  FLOAT	sigmaHeII = 5.36e-20 * LengthUnits;
//  FLOAT	sigmaH2I  = 3.71e-18 * LengthUnits;

  FLOAT sigma, sigmaHI, sigmaHeI, sigmaHeII;
  /* find relevant crossection */
  if ((*PP)->Type >= 0 && (*PP)->Type <= 2)
    sigma = (*PP)->CrossSection * LengthUnits;
  else if ((*PP)->Type == 4) {
    sigmaHI   = FindCrossSection(0, (*PP)->Energy) * LengthUnits;
    sigmaHeI  = FindCrossSection(1, (*PP)->Energy) * LengthUnits;
    sigmaHeII = FindCrossSection(2, (*PP)->Energy) * LengthUnits;
  }
  const float erg_eV = 1.602176e-12;

  /* Secondary ionizations from X-rays */

  int nSecondaryHII = 1, nSecondaryHeIII = 1;
  float xx, heat_factor, kphHI_factor, kphHeII_factor;

  if ((*PP)->Type == 4) {
    nSecondaryHII = (int) floor((*PP)->Energy / 13.6 - 1);
    nSecondaryHeIII = (int) floor((*PP)->Energy / 54.4 - 1);
  }

  /* This controls the splitting condition, where this many rays must
     exist in each cell */

  float RaysPerCell = RadiativeTransferRaysPerCell;

  // Only split photons within this radius if specified
  if (RadiativeTransferSplitPhotonRadius > 0) // && (*PP)->Type == 4)
    SplitWithinRadius = RadiativeTransferSplitPhotonRadius * 
      (3.086e21 / LengthUnits);
  else
    SplitWithinRadius = 2.0;

  /* Convert escape fraction radius into code units */

  float PhotonEscapeRadius[3];
  PhotonEscapeRadius[0] = 0.5 * RadiativeTransferPhotonEscapeRadius * 
    (3.086e21 / LengthUnits);
  PhotonEscapeRadius[1] = RadiativeTransferPhotonEscapeRadius * 
    (3.086e21 / LengthUnits);
  PhotonEscapeRadius[2] = 2.0 * RadiativeTransferPhotonEscapeRadius * 
    (3.086e21 / LengthUnits);

  // speed of light in code units. note this one is independent of a(t)
  float c_cgs = 2.99792e10;
  double c = c_cgs/VelocityUnits, c_inv;

  // Modify the photon propagation speed by this parameter
  c *= RadiativeTransferPropagationSpeedFraction;
  c_inv = 1.0 / c;

  float tau, tauHI, tauHeI, tauHeII;

  // Something to do ?  
  if ((*PP)->Photons <= 0) {
    (*PP)->Photons=-1;
    DeleteMe = TRUE;
    if (DEBUG)
      fprintf(stdout, "called WalkPhotonPackge with empty PhotonPackage "
	      "%x %x %x %x %"GSYM"\n",  (*PP), 
	       (*PP)->PreviousPackage, 
	       (*PP)->NextPackage,  PhotonPackages, 
	      (*PP)->Photons);
    return SUCCESS;
  }
  if (((*PP) == NULL) || ((*PP)->PreviousPackage->NextPackage != (*PP))) {
    fprintf(stdout, "called WalkPhotonPackge with invalid pointer "
	    "%x %x %x %x\n",  (*PP), 
	     (*PP)->PreviousPackage, 
	     (*PP)->PreviousPackage->NextPackage, 
	     PhotonPackages);
    ENZO_FAIL("");
  }

  FLOAT dir_vec[3];
  if (pix2vec_nest((long) (1 << (*PP)->level), (*PP)->ipix, dir_vec)==FAIL) {
    fprintf(stdout,"grid::WalkPhotonPackage:  pix2vec_nest outor %"ISYM" %"ISYM" %"GSYM" %x\n",
	    (long) (1 << (*PP)->level), (*PP)->ipix, (*PP)->Photons, 
	     (*PP) );
    (*PP)->Photons=0;
    ENZO_FAIL("");
  }

  if (DEBUG) fprintf(stderr,"grid::WalkPhotonPackage: %"GSYM" %"GSYM" %"GSYM". \n", dir_vec[0],dir_vec[1],dir_vec[2]);

  // Quantities that help finding which cell index am I in ? 

  FLOAT DomainWidth[3];
  FLOAT dx, dx2;
  for (dim=0; dim<GridRank; dim++) {
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    CellVolume *= CellWidth[i][0];
  }

  dx = CellWidth[0][0];
  dx2 = dx*dx;
  SplitCriteron = dx2 / RaysPerCell;
  SplitCriteronIonized = dx2;
  Volume_inv = 1.0 / CellVolume;
  Area_inv = 1.0 / dx2;

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

  // solid angle associated with package (= 4 Pi/N_package[on this level]) 
  float n_on_this_level = (12.0*pow(4.0,(*PP)->level));
  FLOAT omega_package=4*M_PI/(n_on_this_level);
  float dtheta = sqrt(omega_package);
 
  // Do the computations

  int keep_walking=1;
  FLOAT ddr, dP, EndTime;
  FLOAT nH, nHI, nHeI, nHeII, fH, nH2I, nHI_inv, nHeI_inv, nHeII_inv, xe;
  FLOAT thisDensity, inverse_rho;
  float shield1, shield2, solid_angle, filling_factor, midpoint;
  EndTime = PhotonTime+dtPhoton;

  int dummy;
  FLOAT dr;
  long cindex;
  int  ci;
  int  cj;
  int  ck;
  int gi, gj, gk;

  fH = CoolData.HydrogenFractionByMass;
  DeltaLevel = 0;

  int H2Thin;

  float ConvertToProperNumberDensity = DensityUnits/1.673e-24f;

  /* Get the correct baryon fields (make it pretty) */

  float *density, *HI, *HII, *HeI, *HeII, *H2I;

  density = BaryonField[DensNum];
  HI      = BaryonField[HINum];
  HII     = BaryonField[HINum+1];
  HeI     = BaryonField[HeINum];
  HeII    = BaryonField[HeIINum];
  if (MultiSpecies > 1)
    H2I   = BaryonField[H2INum];

  FLOAT sx = (*PP)->SourcePosition[0];
  FLOAT sy = (*PP)->SourcePosition[1];
  FLOAT sz = (*PP)->SourcePosition[2];

  float EnergyThreshold;
  switch ((*PP)->Type) {
  case 0: EnergyThreshold = 13.6; break;
  case 1: EnergyThreshold = 24.6; break;
  case 2: EnergyThreshold = 54.4; break;
  case 3: EnergyThreshold = 11.2; break;
  default: EnergyThreshold = 0.0;
  }

  FLOAT emission_dt_inv = 1.0 / (*PP)->EmissionTimeInterval;
  FLOAT heat_energy = (*PP)->Energy - EnergyThreshold;
//  FLOAT factor1 = Volume_inv*emission_dt_inv;
//  FLOAT factor2 = factor1 * heat_energy;
  FLOAT factor1 = emission_dt_inv;
  FLOAT factor2 = factor1 * heat_energy;
  FLOAT factor3 = Area_inv*emission_dt_inv;
  FLOAT factor4 = factor3 * heat_energy;

  // For X-ray photons, we do heating and ionization for HI/HeI/HeII in one shot.
  FLOAT factor2_HI, factor2_HeI, factor2_HeII;
  if ((*PP)->Type == 4) {
    factor2_HI   = factor1 * ((*PP)->Energy - 13.6);
    factor2_HeI  = factor1 * ((*PP)->Energy - 24.6);
    factor2_HeII = factor1 * ((*PP)->Energy - 54.4);
  }

  int count = 0;
  FLOAT ux = dir_vec[0];
  FLOAT uy = dir_vec[1];
  FLOAT uz = dir_vec[2];
  //if (ux == 0) ux = ROUNDOFF;
  if (fabs(uy) < ROUNDOFF) uy = sign(uy)*ROUNDOFF; // zeros in y direction possible
  if (fabs(uz) < ROUNDOFF) uz = sign(uz)*ROUNDOFF; // zeros in z direction possible

  FLOAT ux_inv, uy_inv, uz_inv, dr_temp[MAX_DIMENSION], min_dr;
  int ux_dir, uy_dir, uz_dir;
  int direction;

  ux_inv = 1.0 / ux;
  uy_inv = 1.0 / uy;
  uz_inv = 1.0 / uz;

  // Ray direction
  ux_dir = (sign(ux)+1)/2;
  uy_dir = (sign(uy)+1)/2;
  uz_dir = (sign(uz)+1)/2;

  // My Position in coordinates [0..1] 		     
  rx = sx + (*PP)->Radius * ux;
  ry = sy + (*PP)->Radius * uy;
  rz = sz + (*PP)->Radius * uz;

  if ((rx <= GridLeftEdge[0]) || (GridRightEdge[0] <= rx)  ||
      (ry <= GridLeftEdge[1]) || (GridRightEdge[1] <= ry)  ||
      (rz <= GridLeftEdge[2]) || (GridRightEdge[2] <= rz)) {

    switch (GravityBoundaryType) {
    case TopGridPeriodic: 
      FindRootGrid(dummy, Grids0, nGrids0, rx, ry, rz, ux, uy, uz);
      if (dummy >= 0) {
	(*MoveToGrid) = Grids0[dummy];
	DeltaLevel = 0;
	(*PP)->Radius += ROUNDOFF;
	return SUCCESS;
      } else
	// Wrap the photon around the boundary
	if (RadiativeTransferPeriodicBoundary) {
	  if (rx < DomainLeftEdge[0]) {
	    (*PP)->SourcePosition[0] += DomainWidth[0];
	    rx += DomainWidth[0];
	  } else if (rx > DomainRightEdge[0]) {
	    (*PP)->SourcePosition[0] -= DomainWidth[0];
	    rx -= DomainWidth[0];
	  }
	  if (ry < DomainLeftEdge[1]) {
	    (*PP)->SourcePosition[1] += DomainWidth[1];
	    ry += DomainWidth[1];
	  } else if (ry > DomainRightEdge[1]) {
	    (*PP)->SourcePosition[1] -= DomainWidth[1];
	    ry -= DomainWidth[1];
	  }
	  if (rz < DomainLeftEdge[2]) {
	    (*PP)->SourcePosition[2] += DomainWidth[2];
	    rz += DomainWidth[2];
	  } else if (rz > DomainRightEdge[2]) {
	    (*PP)->SourcePosition[2] -= DomainWidth[2];
	    rz -= DomainWidth[2];
	  }
	  FindRootGrid(dummy, Grids0, nGrids0, rx, ry, rz, ux, uy, uz);
	  (*MoveToGrid) = Grids0[dummy];
	  DeltaLevel = 0;
	  (*PP)->Radius += ROUNDOFF;
	  return SUCCESS;
	} else {
	  // PhotonPackage left the box
	  (*PP)->Photons=-1;
	  DeleteMe = TRUE;
	  return SUCCESS;
	}
    case TopGridIsolated: 
      FindRootGrid(dummy, Grids0, nGrids0, rx, ry, rz, ux, uy, uz);
      if (dummy >= 0) {
	(*MoveToGrid) = Grids0[dummy];
	DeltaLevel = 0;
	(*PP)->Radius += ROUNDOFF;
      } else {
	// PhotonPackage left the box
	(*PP)->Photons=-1;
	DeleteMe = TRUE;
      }
      return SUCCESS;
    case SubGridIsolated:
      (*MoveToGrid) = ParentGrid;
      DeltaLevel = -1;
      //      (*PP)->Radius += 0.01*CellWidth[0][0];
      (*PP)->Radius += ROUNDOFF;
      if (DEBUG) 
	fprintf(stdout, "Walk: left grid: sent photon to grid %x\n", 
		ParentGrid);
      return SUCCESS;
    case GravityUndefined:
    default:
      fprintf(stdout, "grid::WalkPhotonPackage: "
	      "GravityBoundaryType = RadiationBoundary undefined %"ISYM".\n",
	      GravityBoundaryType);
      ENZO_FAIL("");
    }
  }

  // current cell ?
  gi = ci = int((rx-GridLeftEdge[0])/CellWidth[0][0]);
  gj = cj = int((ry-GridLeftEdge[1])/CellWidth[1][0]);
  gk = ck = int((rz-GridLeftEdge[2])/CellWidth[2][0]);

  // in floating point:
  FLOAT fx = GridLeftEdge[0] + (FLOAT) (ci) * CellWidth[0][0];
  FLOAT fy = GridLeftEdge[1] + (FLOAT) (cj) * CellWidth[1][0];
  FLOAT fz = GridLeftEdge[2] + (FLOAT) (ck) * CellWidth[2][0];
  
  // on cell boundaries the index will change in negative directions
  if (rx == fx) gi = ci + (sign(ux)-1)/2;
  if (ry == fy) gj = cj + (sign(uy)-1)/2;
  if (rz == fz) gk = ck + (sign(uz)-1)/2; 

  // Mark that this grid has radiation (mainly for the coupled rate solver)
  HasRadiation = TRUE;

  while (keep_walking) {

    cindex = GRIDINDEX(gi,gj,gk);
    oldr = (*PP)->Radius;
    min_dr = 1e20;
       
    // next cell edge crossing radii:
    drx = dr_temp[0] = 
      (GridLeftEdge[0] + 
       ((FLOAT) (gi + ux_dir)) * CellWidth[0][0] - sx) * ux_inv;
    dry = dr_temp[1] = 
      (GridLeftEdge[1] + 
       ((FLOAT) (gj + uy_dir)) * CellWidth[1][0] - sy) * uy_inv;
    drz = dr_temp[2] = 
      (GridLeftEdge[2] + 
       ((FLOAT) (gk + uz_dir)) * CellWidth[2][0] - sz) * uz_inv;

    // the closest one is the one we want
//    if (drx <= min(dry, drz))  r = drx;
//    if (dry <= min(drx, drz))  r = dry;
//    if (drz <= min(drx, dry))  r = drz;

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      if (dr_temp[dim] < min_dr) {
	direction = dim;
	min_dr = dr_temp[dim];
      }

    r = min_dr + ROUNDOFF;
    dr = r - oldr;

    if (dr < 0) {
      printf("dr < 0:   %"GSYM" %"GSYM" %"GSYM"\n"
	     "vec(dr) = %"GSYM" %"GSYM" %"GSYM"\n", dr, min_dr, oldr,
	     dr_temp[0], dr_temp[1], dr_temp[2]);
      (*PP)->Photons = -1;
      DeleteMe = TRUE;
      return SUCCESS;
    }

    // My Position in coordinates [0..1]	     
    rx = sx + r*ux;
    ry = sy + r*uy;
    rz = sz + r*uz;

    if (SubgridMarker[cindex] != CurrentGrid) {
      if (SubgridMarker[cindex] == NULL) {
	(*MoveToGrid) = ParentGrid;
	DeltaLevel = -1;
      } else {
	(*MoveToGrid) = SubgridMarker[cindex];
	DeltaLevel = 1;
	if (DEBUG) 
	  printf("different grid subgrid marker %x %x %"ISYM" %"ISYM" %"ISYM
		 " %"ISYM" %"ISYM" %"FSYM" %"FSYM" %"FSYM"\n",
		 SubgridMarker[cindex], CurrentGrid, cindex, count,
		 gi, gj, gk, rx, ry, rz);
      }
      // move it at least a tiny fraction of the grid cell to not have
      // to worry about round off errors shifting photons back and
      // forth between two grids without doing anything (*PP)->Radius
      // *= (1. + 0.001 * CellWidth[0][0]);
      //      (*PP)->Radius += 0.01*CellWidth[0][0];
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
    solid_angle = r * r * omega_package;
    splitMe = solid_angle > SplitCriteron;

    if (splitMe && 
	r < SplitWithinRadius && 
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
	      (*PP), gi, gj, gk, (*PP)->Radius, dr);

    index = cindex;
    ddr    = dr;

    // nor do we want transport longer than the grid timestep
    ddr    = min(ddr, c*(EndTime-(*PP)->CurrentTime));
    FLOAT cdt;
    cdt = ddr * c_inv;

    // Check for ray merging, only consider a fraction of the ray to
    // make r=PauseRadius and return.
    float fraction;
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

    if (RadiativeTransferPhotonEscapeRadius > 0 && (*PP)->Type == 0) {
      for (i = 0; i < 3; i++) {
	if (r > PhotonEscapeRadius[i] && oldr < PhotonEscapeRadius[i])
	  EscapedPhotonCount[i+1] += (*PP)->Photons;
      } // ENDFOR i
    } // ENDIF PhotonEscapeRadius > 0

    /* Get absorber density and its inverse */

    /* No interpolation */

    switch ((*PP)->Type) {
    case 0:  // HI
      thisDensity = HI[index]*ConvertToProperNumberDensity; 
      break;
    case 1:  // HeI
      thisDensity = 0.25*HeI[index]*ConvertToProperNumberDensity; 
      break;
    case 2:  // HeII
      thisDensity = 0.25*HeII[index]*ConvertToProperNumberDensity; 
      break;
    case 3:  // Lyman-Werner
      if (MultiSpecies > 1)
	thisDensity = H2I[index]*ConvertToProperNumberDensity;
      break;
    case 4:  // X-rays
      nHI   = HI[index]*ConvertToProperNumberDensity; 
      nHeI  = 0.25*HeI[index]*ConvertToProperNumberDensity; 
      nHeII = 0.25*HeII[index]*ConvertToProperNumberDensity; 
      nHI_inv   = 1.0 / nHI;
      nHeI_inv  = 1.0 / nHeI;
      nHeII_inv = 1.0 / nHeII;
      break;
    } // ENDSWITCH type

    switch ((*PP)->Type) {
    case 0:  // HI
      kphNum = kphHINum;
      gammaNum = gammaHINum;
      break;
    case 1:  // HeI
      kphNum = kphHeINum;
      gammaNum = gammaHeINum;
      break;
    case 2:  // HeII
      kphNum = kphHeIINum;
      gammaNum = gammaHeIINum;
      break;
    } // ENDSWITCH

//    //filling_factor = 0.5*(oldr*oldr + r*r) * ddr * omega_package * Volume_inv;
//    //filling_factor = r*r * ddr * omega_package * Volume_inv;
    mx = fabs(sx + (oldr + 0.5*ddr - ROUNDOFF) * ux -
	      (GridLeftEdge[0] + (gi+0.5)*CellWidth[0][0]));
    my = fabs(sy + (oldr + 0.5*ddr - ROUNDOFF) * uy -
	      (GridLeftEdge[1] + (gj+0.5)*CellWidth[1][0]));
    mz = fabs(sz + (oldr + 0.5*ddr - ROUNDOFF) * uz -
	      (GridLeftEdge[2] + (gk+0.5)*CellWidth[2][0]));
    float nearest_edge = max(max(mx, my), mz);
    slice_factor = min(0.5 + (0.5*dx-nearest_edge) / (dtheta*r), 1);
    slice_factor2 = slice_factor * slice_factor;
//    printf("mx, my, mz = %g %g %g, dtheta = %g, dtheta*r = %g, slice_factor = %g\n",
//	   mx, my, mz, dtheta, dtheta*r, slice_factor);
//    filling_factor = 0.5*(oldr*oldr + r*r) * ddr * omega_package * Volume_inv *
//      slice_factor;
//    filling_factor = (1.0/3.0) * (r*r*r - oldr*oldr*oldr) * omega_package *
//      Volume_inv * slice_factor;
//    filling_factor = slice_factor*slice_factor*r*r*dtheta*dtheta * ddr * Volume_inv;

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

    /* H2I dissociation */

    if ((*PP)->Type == 3) {

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
      BaryonField[kdissH2INum][index] += (*PP)->Photons * sigma * factor3 * 
	(ddr * omega_package * r * r * Volume_inv);

    } else if ((*PP)->Type >= 0 && (*PP)->Type <= 2) {
	
      /* HI/HeI/HeII ionization and heating */
      
      // optical depth of ray segment
      tau = thisDensity*sigma*ddr;
      // at most use all photons for photo-ionizations
      if (tau > 2.e1) dP = (1.0+ROUNDOFF) * (*PP)->Photons;
      else if (tau > 1.e-4) 
	dP = min((*PP)->Photons*(1-expf(-tau)), (*PP)->Photons);
      else
	dP = min((*PP)->Photons*tau, (*PP)->Photons);

      // contributions to the photoionization rate is over whole timestep
      BaryonField[kphNum][index] += dP*factor1*slice_factor2;
	
      // the heating rate is just the number of photo ionizations times
      // the excess energy units here are  eV/s/cm^3 *TimeUnits. 
      BaryonField[gammaNum][index] += dP*factor2*slice_factor2;

    }
    else if ((*PP)->Type == 4) {

      /* X-rays (HI/HeI/HeII all in one!) */

      // optical depth of ray segment
      tauHI = nHI*sigmaHI*ddr;
      tauHeI = nHeI*sigmaHeI*ddr;
      tauHeII = nHeII*sigmaHeII*ddr;

      // Secondary ionizations (minus heating)
      
      if (RadiationXRaySecondaryIon) {
	xx = max(HII[index] / (HI[index] + HII[index]), 1e-4);
	heat_factor    = 0.9971 * (1 - powf(1 - powf(xx, 0.2663f), 1.3163));
	kphHI_factor   = 0.3908 * powf(1 - powf(xx, 0.4092f), 1.7592f);
	kphHeII_factor = 0.0554 * powf(1 - powf(xx, 0.4614f), 1.6660f);
	kphHI_factor   *= nSecondaryHII;
	kphHeII_factor *= nSecondaryHeIII;
      } else {
	heat_factor    = 1;
	kphHI_factor   = 1;
	kphHeII_factor = 1;
      }

      // HI
      if (tauHI > 2.e1) dP = 1.001f * (*PP)->Photons;
      else if (tauHI > 1.e-3) dP = min((*PP)->Photons*(1-expf(-tauHI)), 
				       (*PP)->Photons);
      else dP = min((*PP)->Photons*tauHI, (*PP)->Photons);
      BaryonField[kphHINum][index] += kphHI_factor*dP*factor1;
      BaryonField[gammaHINum][index] += heat_factor*dP*factor2_HI;
      (*PP)->Photons -= dP;

      // HeI
      if (tauHeI > 2.e1) dP = 1.001f * (*PP)->Photons;
      else if (tauHeI > 1.e-3) dP = min((*PP)->Photons*(1-expf(-tauHeI)), 
				       (*PP)->Photons);
      else dP = min((*PP)->Photons*tauHeI, (*PP)->Photons);
      BaryonField[kphHeINum][index] += dP*factor1;
      BaryonField[gammaHeINum][index] += heat_factor*dP*factor2_HeI;
      (*PP)->Photons -= dP;

      // HeII
      if (tauHeII > 2.e1) dP = 1.001f * (*PP)->Photons;
      else if (tauHeII > 1.e-3) dP = min((*PP)->Photons*(1-expf(-tauHeII)), 
				       (*PP)->Photons);
      else dP = min((*PP)->Photons*tauHeII, (*PP)->Photons);
      BaryonField[kphHeIINum][index] += kphHeII_factor*dP*factor1;
      BaryonField[gammaHeIINum][index] += heat_factor*dP*factor2_HeII;

    } else
      dP = 0;

    // acceleration due to radiation pressure
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

    BaryonField[kphHeIINum][index] += 1;

    // return in case we're pausing to merge
    if (PauseMe)
      return SUCCESS;

    // return in case we're out of photons
    if ((*PP)->Photons < tiny_number) {
      if (DEBUG>1) 
	fprintf(stderr, "PP-Photons: %"GSYM"  PP->Radius: %"GSYM" "
		"PP->CurrentTime: \n",
		(*PP)->Photons, (*PP)->Radius, (*PP)->CurrentTime);
      if (DEBUG>1) 
	fprintf(stderr, "\tdP: %"GSYM"\tddr: %"GSYM"\t cdt: %"GSYM"\n", 
		dP, ddr, cdt);
      (*PP)->Photons = -1;
      DeleteMe = TRUE;
      return SUCCESS;
    }
    
    // are we done ? 
    if (((*PP)->CurrentTime) >= EndTime) {
      return SUCCESS;
    }

    //   if ((gi < 0) || ((GridEndIndex[0]-GridStartIndex[0]) < gi) ||
    //	(gj < 0) || ((GridEndIndex[1]-GridStartIndex[1]) < gj) ||
    //	(gk < 0) || ((GridEndIndex[2]-GridStartIndex[2]) < gk)) {
    if ((rx <= GridLeftEdge[0]) || (GridRightEdge[0] <= rx)  ||
	(ry <= GridLeftEdge[1]) || (GridRightEdge[1] <= ry)  ||
	(rz <= GridLeftEdge[2]) || (GridRightEdge[2] <= rz)) {
      switch (GravityBoundaryType) {
      case TopGridPeriodic: 
	FindRootGrid(dummy, Grids0, nGrids0, rx, ry, rz, ux, uy, uz);
	if (dummy >= 0) {
	  (*MoveToGrid) = Grids0[dummy];
	  DeltaLevel = 0;
	  (*PP)->Radius += ROUNDOFF;
	  return SUCCESS;
	} else
	  // Wrap the photon around the boundary
	  if (RadiativeTransferPeriodicBoundary) {
	    if (rx < DomainLeftEdge[0]) {
	      (*PP)->SourcePosition[0] += DomainWidth[0];
	      rx += DomainWidth[0];
	    } else if (rx > DomainRightEdge[0]) {
	      (*PP)->SourcePosition[0] -= DomainWidth[0];
	      rx -= DomainWidth[0];
	    }
	    if (ry < DomainLeftEdge[1]) {
	      (*PP)->SourcePosition[1] += DomainWidth[1];
	      ry += DomainWidth[1];
	    } else if (ry > DomainRightEdge[1]) {
	      (*PP)->SourcePosition[1] -= DomainWidth[1];
	      ry -= DomainWidth[1];
	    }
	    if (rz < DomainLeftEdge[2]) {
	      (*PP)->SourcePosition[2] += DomainWidth[2];
	      rz += DomainWidth[2];
	    } else if (rz > DomainRightEdge[2]) {
	      (*PP)->SourcePosition[2] -= DomainWidth[2];
	      rz -= DomainWidth[2];
	    }
	    FindRootGrid(dummy, Grids0, nGrids0, rx, ry, rz, ux, uy, uz);
	    (*MoveToGrid) = Grids0[dummy];
	    DeltaLevel = 0;
	    (*PP)->Radius += ROUNDOFF;
	    return SUCCESS;
	  } else {
	    // PhotonPackage left the box
	    (*PP)->Photons=-1;
	    DeleteMe = TRUE;
	    return SUCCESS;
	  }
      case TopGridIsolated: 
	FindRootGrid(dummy, Grids0, nGrids0, rx, ry, rz, ux, uy, uz);
	if (dummy >= 0) {
	  (*MoveToGrid) = Grids0[dummy];
	  DeltaLevel = 0;
	  (*PP)->Radius += ROUNDOFF;
	} else {
	  // PhotonPackage left the box
	  (*PP)->Photons=-1;
	  DeleteMe = TRUE;
	}
	return SUCCESS;
      case SubGridIsolated:
	(*MoveToGrid) = ParentGrid;
	DeltaLevel = -1;
	//	(*PP)->Radius += 0.01*CellWidth[0][0];
	(*PP)->Radius += ROUNDOFF;
	if (DEBUG) 
	  fprintf(stdout, "Walk: left grid: sent photon to grid %x\n", 
		  ParentGrid);
	return SUCCESS;
      case GravityUndefined:
      default:
	fprintf(stdout, "grid::WalkPhotonPackage: GravityBoundaryType = "
		"RadiationBoundary undefined %"ISYM".\n",
		GravityBoundaryType);
	ENZO_FAIL("");
      }
    } // ENDIF photon leaving grid

    count++;

    switch (direction) {
    case 0:
      gi += sign(ux);
      break;
    case 1:
      gj += sign(uy);
      break;
    case 2:
      gk += sign(uz);
      break;
    }
    
  } // while keep walking

  return SUCCESS;
}
