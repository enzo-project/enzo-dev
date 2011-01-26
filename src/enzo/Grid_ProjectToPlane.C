/***********************************************************************
/
/  GRID CLASS (PROJECT GRID VALUES TO A PLANE)
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:
/
/  PURPOSE:
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
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "CosmologyParameters.h"
#include "Grid.h"
 
extern "C" void FORTRAN_NAME(projplane)(
          float *grid1, float *grid2, float *flaggrid, int *iflag,
              int *ismooth,
          int *gdim1, int *gdim2, int *gdim3, FLOAT *gcellsize,
          float *plane, int *pdim1, int *pdim2, FLOAT *pcellsize,
          int *projdim, int *ifield, float *weight,
          FLOAT *gleft, FLOAT *gfarleft, FLOAT *gright, FLOAT *pleft,
              FLOAT *pright,
          int *npstart, int *npend, float *fracleft, float *fracright);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
 
int grid::ProjectToPlane(FLOAT ProjectedFieldLeftEdge[],
			 FLOAT ProjectedFieldRightEdge[],
			 int ProjectedFieldDims[], float *ProjectedField[],
			 int ProjectionDimension, int ProjectionSmooth,
			 int NumberOfProjectedFields, int level,
			 int XrayUseLookupTable, float XrayLowerCutoffkeV,
			 float XrayUpperCutoffkeV, char *XrayTableFileName)
{
 
  if (GridRank != 3)
    return SUCCESS;
 
  if (BaryonField[NumberOfBaryonFields] == NULL && level >= 0) {
    ENZO_FAIL("UNDER_SUBGRID_FLAG field not set.\n");
  }
 
  if (SelfGravity && GravityResolution != 1) {
    ENZO_FAIL("ProjectToPlane assumes GravityResolution == 1.\n");
  }
 
  /* Declarations */
 
  int i, j, k, dim, start, stop, Index[MAX_DIMENSION];
  float LeftCellFraction, RightCellFraction, ConversionFactor;
  float *temperature = NULL;
  float *xray_emissivity = NULL, *first_field, *second_field;
 
  /* Check To see if grid overlaps the projected field. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    if (GridLeftEdge[dim] > ProjectedFieldRightEdge[dim] ||
	GridRightEdge[dim] < ProjectedFieldLeftEdge[dim])
      return SUCCESS;
    size *= GridDimension[dim];
  }
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (NumberOfBaryonFields > 0)
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }
 
  /* Find metallicity field and set flag. */
 
  int MetallicityField = FALSE, MetalNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;
 
  /* Find the start and stop indicies in the ProjectionDimension of this
     grid for the projected region. */
 
  start = max(int((ProjectedFieldLeftEdge[ProjectionDimension] -
	 	   GridLeftEdge[ProjectionDimension]) /
	          CellWidth[ProjectionDimension][0]), 0);
  stop  = min(int((ProjectedFieldRightEdge[ProjectionDimension] -
	           GridLeftEdge[ProjectionDimension]) /
	          CellWidth[ProjectionDimension][0]),
	      GridEndIndex[ProjectionDimension] -
	      GridStartIndex[ProjectionDimension]);
 
  LeftCellFraction = min(1.0 - ((ProjectedFieldLeftEdge[ProjectionDimension] -
				 GridLeftEdge[ProjectionDimension]) /
				CellWidth[ProjectionDimension][0] - start), 1);
  RightCellFraction = min((ProjectedFieldRightEdge[ProjectionDimension] -
			   GridLeftEdge[ProjectionDimension]) /
			  CellWidth[ProjectionDimension][0] - stop, 1);
 
  start += GridStartIndex[ProjectionDimension];
  stop += GridStartIndex[ProjectionDimension];
  if (debug) printf("ProjectToGrid: start = %"ISYM"/%"ISYM" (%5.3"FSYM")  stop = %"ISYM"/%"ISYM" (%5.3"FSYM")  GridLeft/Right = %5.3"PSYM"/%5.3"PSYM"\n",
		start, GridStartIndex[ProjectionDimension], LeftCellFraction,
		stop, GridEndIndex[ProjectionDimension], RightCellFraction,
		    GridLeftEdge[ProjectionDimension],
		    GridRightEdge[ProjectionDimension]);
 
  /* Compute the volume factor for each projected cell (projected cell width
     for the plane dimensions and the grid cell width for the projected
     dimension). */
 
  float CellLength = CellWidth[ProjectionDimension][0];
 
  /* Set the Conversion factors for Density and X-rays.  If using comoving
     coordinates use solar masses and Mpc as the intrinsic units.
     Note: The X-ray units have been multiplied by 1.0e-20 to stop overflow.
     Note: The temperature field already has units of K. */
 
  float DensityConversion, XrayConversion, TempXrayConversion;
  DensityConversion = XrayConversion = TempXrayConversion = CellLength;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
  if (ComovingCoordinates) {
    const double SolarMass = 1.989e33, Mpc = 3.0824e24;
    DensityConversion *= float(double(DensityUnits)*double(LengthUnits)/
			       SolarMass*Mpc*Mpc);
    XrayConversion *= float(1.0e-20*POW(double(DensityUnits),2)
			      /POW(SolarMass,2)*POW(Mpc,6)
			    *double(LengthUnits)/Mpc);
    TempXrayConversion = XrayConversion;
  }
 
  /* Set some required numbers for the FORTRAN call. */
 
  int One = 1;
  FLOAT GridFarLeftEdge[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    GridFarLeftEdge[dim] = CellLeftEdge[dim][0];
  int adim = (ProjectionDimension == 0) ? 1 : 0;
  int bdim = (ProjectionDimension == 2) ? 1 : 2;
 
  /* Compute the projected field cell size. */
 
  FLOAT ProjectedFieldCellSize = (ProjectedFieldRightEdge[0] -
                                  ProjectedFieldLeftEdge[0])/
                                  FLOAT(ProjectedFieldDims[0]);
 
  /* If level < 0, then just do the star particle stuff. */
 
  if (level < 0)
    goto StarParticleLabel;
 
  /* Compute the temperature. */
 
  temperature = new float[size];
  if (NumberOfBaryonFields > 0) {
    if (this->ComputeTemperatureField(temperature) == FAIL) {
      ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
    }
 
    /* Set the temperature to zero wherever the baryon density is. */
 
    for (i = 0; i < size; i++)
      if (BaryonField[0][i] == 0) temperature[i] = 0;
  }
 
  /* 1) baryon density. */
 
  if (NumberOfBaryonFields > 0)
  FORTRAN_NAME(projplane)(BaryonField[0], NULL,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[0], ProjectedFieldDims+adim,
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &One, &DensityConversion,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
 
  /* 2) 'X-ray' luminosity */
 
  int ProjType;
  if (NumberOfBaryonFields > 0) {
 
    if (XrayUseLookupTable) {
 
    /* Compute the real X-ray emissivity.  Units will be in 10^-23 erg/cm^2/s
       (actually, the units depend on XrayTableFileName). */
 
      xray_emissivity = new float[size];
      this->ComputeXrayEmissivity(temperature, xray_emissivity,
				  XrayLowerCutoffkeV, XrayUpperCutoffkeV,
				  XrayTableFileName);
      first_field = xray_emissivity;
      second_field = NULL;
      ProjType = 1;
      ConversionFactor = LengthUnits*CellLength;
 
    } else {
 
    /* Compute rho^2*T^1/2, a measure of bolometric free-free emissivity. */
 
      first_field = temperature;
      second_field = BaryonField[0];
      ProjType = 2;
      ConversionFactor = XrayConversion;
 
    }
 
    FORTRAN_NAME(projplane)(first_field, second_field,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[1], ProjectedFieldDims+adim,
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &ProjType, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
  }
 
  /* 3) Dark matter density. */
 
  if (GravitatingMassFieldParticles != NULL) {
 
    float *temp = new float[size];
    for (i = 0; i < size; i++)
      temp[i] = 0.0;
    int Offset[MAX_DIMENSION], bfindex, dmindex;
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      Offset[dim] = -nint((GridFarLeftEdge[dim] -
                           GravitatingMassFieldParticlesLeftEdge[dim])/
                          CellWidth[dim][0]);
//      printf("Offset[%"ISYM"] = %"ISYM"  GFLE = %"GSYM"  GMFPLE = %"GSYM"\n", dim, Offset[dim],
//	     GridFarLeftEdge[dim],GravitatingMassFieldParticlesLeftEdge[dim]);
    }
 
    /* Copy grid-like region of GravitatingMassFieldParticles to temp field. */
 
    for (k = max(0,Offset[2]); k < GridDimension[2]-max(0,Offset[2]); k++)
      for (j = max(0,Offset[1]); j < GridDimension[1]-max(0,Offset[1]); j++) {
        bfindex = (k*GridDimension[1] + j)*GridDimension[0] + max(0,Offset[0]);
        dmindex = ((k - Offset[2])*GravitatingMassFieldParticlesDimension[1] +
                   (j - Offset[1]))*GravitatingMassFieldParticlesDimension[0] +
		   (max(0, Offset[0]) - Offset[0]);
	for (i = max(0,Offset[0]); i < GridDimension[0]-max(0,Offset[0]);
	     i++, bfindex++, dmindex++)
	  temp[bfindex] = GravitatingMassFieldParticles[dmindex];
      }
 
    FORTRAN_NAME(projplane)(temp, NULL,
                              BaryonField[NumberOfBaryonFields], &One,
                              &ProjectionSmooth,
                            GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                            ProjectedField[2], ProjectedFieldDims+adim,
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                            &ProjectionDimension, &One, &DensityConversion,
                            GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                            &start, &stop,
			    &LeftCellFraction, &RightCellFraction);
 
    delete temp;
  }
 
  /* 4) Temperature weighted by 'X-ray' luminosity */
 
  if (NumberOfBaryonFields > 0) {
 
    if (XrayUseLookupTable) {
      for (i = 0; i < size; i++)
	xray_emissivity[i] *= temperature[i];
      first_field = xray_emissivity;
      second_field = NULL;
      ProjType = 1;
      ConversionFactor = LengthUnits*CellLength;
    } else {
      first_field = temperature;
      second_field = BaryonField[0];
      ProjType = 3;
      ConversionFactor = TempXrayConversion;
    }
 
    FORTRAN_NAME(projplane)(first_field, second_field,
			    BaryonField[NumberOfBaryonFields], &One,
			       &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[3], ProjectedFieldDims+adim,
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &ProjType, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
 
    delete [] xray_emissivity;
 
  }
 
  /* 5) Level depth. */
 
  int PFStart[MAX_DIMENSION], PFStop[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
 
    PFStart[dim] = max(nint((GridLeftEdge[dim] - ProjectedFieldLeftEdge[dim])/
			     ProjectedFieldCellSize), 0);
    PFStop[dim]  = min(nint((GridRightEdge[dim] - ProjectedFieldLeftEdge[dim])/
		 	     ProjectedFieldCellSize),
		       ProjectedFieldDims[dim]) - 1;
  }
  /*  printf("PFStart = %"ISYM" %"ISYM"  PFStop = %"ISYM" %"ISYM", a/bdim = %"ISYM"/%"ISYM" level = %"GSYM"\n",
	 PFStart[adim], PFStart[bdim], PFStop[adim], PFStop[bdim],
	 adim, bdim, float(level)); */
 
  for (j = PFStart[bdim]; j <= PFStop[bdim]; j++)
    for (i = PFStart[adim]; i <= PFStop[adim]; i++)
      ProjectedField[4][j*ProjectedFieldDims[adim] + i] = float(level);
 
  if (NumberOfBaryonFields > 0 && ComovingCoordinates) {
 
    /* 6) y parameter for the SZ effect
       electron pressure (assumes fh=0.24, all ionized) */
 
    float *sz = new float[size];
    for (i = 0; i < size; i++)
      sz[i] = BaryonField[DensNum][i]*temperature[i];
 
    double sigma_thompson = 6.65e-25, mh = 1.67e-24, me = 9.11e-28,
           kboltz = 1.38e-16, clight = 3.00e10, csquared = 8.99e20;
 
    ConversionFactor = double(DensityUnits)*0.88/mh
                       *kboltz/(me*csquared)*sigma_thompson
		       *double(LengthUnits)*CellLength;
 
    if (NumberOfBaryonFields > 0)
      FORTRAN_NAME(projplane)(sz, NULL, BaryonField[NumberOfBaryonFields],
			      &One, &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[5], ProjectedFieldDims+adim,
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &One, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
 
 
    /* 7) Delta T/T for the bulk motion SZ effect
       v_par * rho (assumes fh=0.24, all ionized) */
 
    for (i = 0; i < size; i++)
      sz[i] = BaryonField[Vel1Num+ProjectionDimension][i] *
              BaryonField[DensNum][i];
 
    ConversionFactor = 1.0*double(VelocityUnits)*double(DensityUnits)*0.88
                       *sigma_thompson/mh/clight
		       *double(LengthUnits)*CellLength;
 
    if (NumberOfBaryonFields > 0)
      FORTRAN_NAME(projplane)(sz, NULL, BaryonField[NumberOfBaryonFields],
			         &One, &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[6], ProjectedFieldDims+adim,
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &One, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
 
 
    delete sz;
 
  } // end: SZ effects
 
  /* 8) metallicity. */
 
  if (MetallicityField)
    FORTRAN_NAME(projplane)(BaryonField[MetalNum], NULL,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[7], ProjectedFieldDims+adim,
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &One, &DensityConversion,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
 
  /* 10/11/12/13/14/15) OVII/OVIII/SiXIII/SiXIV/FeXXV/FeXXVI column density. */
 
  if (NumberOfBaryonFields > 0 && NumberOfProjectedFields > 9) {
    float *elemental_density = new float[size];
    ConversionFactor = LengthUnits*CellLength;
 
    for (i = 0; i < NumberOfProjectedFields-9; i++) {
//      if (this->ComputeElementalDensity(temperature, elemental_density, i)
      if (this->ComputeElementalDensity(temperature, elemental_density, 6)
	  == FAIL) {
	ENZO_FAIL("Error in grid->ComputeElementalDensity\n");
      }
      FORTRAN_NAME(projplane)(elemental_density, NULL,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[9+i], ProjectedFieldDims+adim,
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &One, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
    }
 
    delete [] elemental_density;
  }
 
  /* 9) Star particle field (NGP interpolation). */
 
 StarParticleLabel:
 
  float weight = DensityConversion*
    POW(CellWidth[0][0]/ProjectedFieldCellSize, FLOAT(2.0));
  if (StarParticleCreation && NumberOfParticleAttributes > 0)
    for (i = 0; i < NumberOfParticles; i++)
      if (ParticleAttribute[0][i] > 0) {
	for (dim = 0; dim < GridRank; dim++) {
	  Index[dim] = int((ParticlePosition[dim][i] -
			    ProjectedFieldLeftEdge[dim])/
			   ProjectedFieldCellSize + 1.0) - 1;
	  if (Index[dim] < 0 || Index[dim] >= ProjectedFieldDims[dim])
	    break;
	}
	if (dim == GridRank)

	  ProjectedField[8][Index[bdim]*ProjectedFieldDims[adim]+Index[adim]]
	    += ParticleMass[i]*weight;
      }
 
  /* Clean up. */
 
  delete [] temperature;
 
  return SUCCESS;
}
