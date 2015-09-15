/***********************************************************************
/
/  GRID CLASS (PROJECT GRID VALUES TO A PLANE)
/
/  written by: Greg Bryan
/  date:       February, 1996
/  modified1:  John Wise (June, 2008)
/              - changed fields and output to HDF5
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
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


int grid::ProjectToPlane2(FLOAT ProjectedFieldLeftEdge[], 
			  FLOAT ProjectedFieldRightEdge[],
			  int ProjectedFieldDims[], float *ProjectedField[], 
			  int ProjectionDimension, int ProjectionSmooth,
			  int NumberOfProjectedFields, int level,
			  int MetalLinesUseLookupTable, char *MetalLinesFilename)
{
  /* Return if the data is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Projection only allowed for 3D simulations */

  if (GridRank != 3) 
    return SUCCESS;
  
  if (BaryonField[NumberOfBaryonFields] == NULL && level >= 0)
    ENZO_FAIL("UNDER_SUBGRID_FLAG field not set.");

  if (SelfGravity && GravityResolution != 1)
    ENZO_FAIL("ProjectToPlane assumes GravityResolution == 1.");

  /* Declarations */

  const int MetalLumField = 17;
  const int NumberOfMetalLines = 11;
  const int NumberOfLuminosityFields = 18;
  int i, j, k, n, dim, start, stop, Index[MAX_DIMENSION];
  float LeftCellFraction, RightCellFraction, ConversionFactor;
  float *temp_field = NULL, *temperature = NULL;
  float *spin_temp = NULL, *bright_temp = NULL, *all_luminosities = NULL;
  float *first_field, *second_field;
  float kappa, gamma_p, gamma_e, C_e, C_H, C_p, y21, v_i, ri;
  

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
    IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
			       Vel3Num, TENum);

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;
  if (NumberOfBaryonFields > 0 && MultiSpecies)
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
			  HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  /* Find metallicity field and set flag. */

  int MetallicityField = FALSE, MetalNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) 
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

  /* Find SN Colour field and set flag */

  int SNColourField = FALSE, SNColourNum;
  if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields)) 
      != -1)
    SNColourField = TRUE;
  else
    SNColourNum = 0;

  if ((SNColourField || MetallicityField) && CoolData.metals != NULL)
    MetalCooling = JHW_METAL_COOLING;

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
  if (debug) 
    printf("ProjectToGrid: start = %d/%d (%5.3f)  stop = %d/%d (%5.3f)  "
	   "GridLeft/Right = %5.3"FSYM"/%5.3"FSYM"\n",
	   start, GridStartIndex[ProjectionDimension], LeftCellFraction,
	   stop, GridEndIndex[ProjectionDimension], RightCellFraction,
	   GridLeftEdge[ProjectionDimension],
	   GridRightEdge[ProjectionDimension]);

  /* Compute the volume factor for each projected cell (projected cell width
     for the plane dimensions and the grid cell width for the projected
     dimension). */

  float CellLength = CellWidth[ProjectionDimension][0];
  double CellVolume;

  /* Set the Conversion factors for Density and X-rays.  If using comoving
     coordinates use solar masses and Mpc as the intrinsic units. 
     Note: The X-ray units have been multiplied by 1.0e-20 to stop overflow.
     Note: The temperature field already has units of K. */

  float DensityConversion, XrayConversion, TempXrayConversion, 
    LuminosityConversion, dom;
  DensityConversion = XrayConversion = TempXrayConversion = 
    LuminosityConversion = CellLength;
  float TemperatureUnits, DensityUnits, LengthUnits, 
        VelocityUnits, TimeUnits;
  double sigma_thompson = 6.65e-25, mh = 1.67e-24, me = 9.11e-28,
    kboltz = 1.38e-16, clight = 3.00e10, csquared = 8.99e20;
  const double SolarMass = 1.989e33, Mpc = 3.0824e24;
  
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  FLOAT a=1, dadt, CurrentRedshift = 0.0;
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;
  }


//  DensityConversion *= float(double(DensityUnits)*double(LengthUnits)/
//			     SolarMass*Mpc*Mpc);
  dom = float(double(DensityUnits) / mh);
  DensityConversion *= dom;
  //DensityConversion *= float(double(dom)*double(LengthUnits));
  XrayConversion *= float(1.0e-20*pow(double(DensityUnits),2)
			  /pow(SolarMass,2)*pow(Mpc,6)
			  *double(LengthUnits)/Mpc);
  TempXrayConversion = XrayConversion;
  LuminosityConversion *= LengthUnits;
  CellVolume = pow(double(CellLength) * double(LengthUnits), 3);
  //LuminosityConversion = float(1e-40*CellVolume);

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

  temp_field = new float[size];
  spin_temp = new float[size];
  bright_temp = new float[size];
  all_luminosities = new float[NumberOfLuminosityFields*size];

  /* Compute the temperature. */

  temperature = new float[size];

  if (NumberOfBaryonFields > 0) {

    this->ComputeTemperatureField(temperature);

    /* Set the temperature to zero wherever the baryon density is. */

    for (i = 0; i < size; i++)
      if (BaryonField[0][i] == 0) temperature[i] = 0;

  } // ENDIF NumberOfBaryons > 0

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

  /* 2) Temperature weighted by rho^2. */

  int ProjType;
  if (NumberOfBaryonFields > 0) {

    first_field = temperature;
    second_field = BaryonField[0];
    ProjType = 4;
    ConversionFactor = DensityConversion;//*DensityConversion;

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

  /* 3) baryon density squared */

  if (NumberOfBaryonFields > 0) {

//    for (i = 0; i < size; i++)
//      temp_field[i] = BaryonField[0][i] * BaryonField[0][i];
    
    first_field = BaryonField[0];
    second_field = BaryonField[0];
    ProjType = 4;
    ConversionFactor = DensityConversion;//*DensityConversion;

    FORTRAN_NAME(projplane)(first_field, second_field,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[2], ProjectedFieldDims+adim, 
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &ProjType, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
  }

  /* 4) SN Colour weighted by density */

  if (SNColourField == TRUE) {
    
    for (i = 0; i < size; i++)
      temp_field[i] = BaryonField[SNColourNum][i] / BaryonField[DensNum][i];

    first_field = temp_field;
    second_field = BaryonField[0];
    //second_field = NULL;
    ProjType = 4;
    ConversionFactor = DensityConversion/CoolData.SolarMetalFractionByMass;
    //ConversionFactor = One;

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
  }

  /* 5) Smoothed Dark Matter Density */

//  if (SmoothedDarkMatterDensity != NULL) {
//
//    first_field = SmoothedDarkMatterDensity;
//    second_field = NULL;
//    ProjType = 1;
//    ConversionFactor = DensityConversion;
//
//    FORTRAN_NAME(projplane)(first_field, second_field,
//			    BaryonField[NumberOfBaryonFields], &One,
//			       &ProjectionSmooth,
//			  GridDimension, GridDimension+1,
//                              GridDimension+2, CellWidth[0],
//                          ProjectedField[4], ProjectedFieldDims+adim, 
//                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
//                          &ProjectionDimension, &ProjType, &ConversionFactor,
//                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
//                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
//                          &start, &stop, &LeftCellFraction,&RightCellFraction);
//
//  }

  /* 6) Electron Fraction weighted by density^2 */

  if (MultiSpecies) {

    for (i = 0; i < size; i++)
      temp_field[i] = BaryonField[DeNum][i] / BaryonField[DensNum][i];
    
    first_field = temp_field;
    second_field = BaryonField[0];
    ProjType = 4;
    ConversionFactor = DensityConversion;//*DensityConversion;

    FORTRAN_NAME(projplane)(first_field, second_field,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[5], ProjectedFieldDims+adim, 
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &ProjType, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
  }
  /* 7) H2 Fraction weighted by density^2 */

  if (MultiSpecies > 1) {

    for (i = 0; i < size; i++)
      temp_field[i] = (BaryonField[H2INum][i] + BaryonField[H2IINum][i]) / 
	BaryonField[DensNum][i];

    first_field = temp_field;
    second_field = BaryonField[0];
    ProjType = 4;
    ConversionFactor = DensityConversion;

    FORTRAN_NAME(projplane)(first_field, second_field,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[6], ProjectedFieldDims+adim, 
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &ProjType, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);

  }

  /* 8) Spin Temperature weighted by density */

  if (MultiSpecies && ComovingCoordinates) {

    const float t_star = 0.068;  // Energy difference between hyperfine levels
    const float A_em = 2.85e-15; // Spont. emission rate [/s]
    float nu0, phi, Delta_nu, Delta_nuD, hz, prefactor, Tcmb, high_dt, logT;

    // 21 cm redshifted to box redshift [Hz]
    nu0 = 1.4204e9 / (1 + CurrentRedshift);
    hz = 3.24044e-18 * HubbleConstantNow * 
      sqrt(OmegaMatterNow * pow(1+CurrentRedshift, 3) + OmegaLambdaNow);
    Tcmb = 2.723*(1+CurrentRedshift);
    high_dt = 1.14e5 /
      sqrt(OmegaMatterNow * HubbleConstantNow * HubbleConstantNow / 0.147) *
      pow((1 + CurrentRedshift)/11, -2.5) * dom;

    prefactor = t_star * CellLength * LengthUnits * dom *
      (3 * (clight*clight) * A_em) / (32 * M_PI * nu0 * nu0);

    for (i = 0; i < size; i++) {

      if (BaryonField[DeNum][i] < tiny_number) {
	spin_temp[i] = 0;
	bright_temp[i] = 0;
	continue;
      }

      kappa = 3.1e-11 * pow(temperature[i], 0.357) * exp(-32.0 / temperature[i]);
      gamma_e = -9.607 + 0.5 * log(temperature[i]) * 
	exp(-pow(log(temperature[i]), 4.5) / 1800.0);
      gamma_e = pow(10.0, gamma_e);
      gamma_p = 3.2 * kappa;

      C_H = BaryonField[HINum][i] * dom * kappa;
      C_e = BaryonField[DeNum][i] * dom * gamma_e;
      C_p = BaryonField[DeNum][i] * dom * gamma_p;
	
      y21 = (t_star / (A_em*temperature[i])) * (C_H + C_e + C_p);
      spin_temp[i] = (t_star + Tcmb + 
		      y21*temperature[i]) / (1+y21);
	
      /* Discretization as outlined in Kuhlen et al. (2005) */

//      ri = 0.0;  // 0 for now.  Should be proper distance from midplane.
//      v_i = VelocityUnits * BaryonField[Vel1Num+ProjectionDimension][i] + hz*ri;
//      Delta_nu = nu0 * v_i / clight;
//      Delta_nuD = nu0 * sqrt(2*kboltz*temperature[i]/(mh*clight*clight));
//      phi = exp(-( (Delta_nu*Delta_nu) / (Delta_nuD*Delta_nuD) )) / 
//	(1.77245*Delta_nuD); // line profile
//      Delta_tau[i] = prefactor * BaryonField[HINum][i] * phi / spin_temp[i];

      /* Instead, let's use an approximation for differential
	 brightness temperature */

      bright_temp[i] = high_dt * BaryonField[HINum][i] * 
	(1.0 - Tcmb / spin_temp[i]);

    } // ENDFOR size

    first_field = spin_temp;
    second_field = BaryonField[0];
    ProjType = 4;
    ConversionFactor = DensityConversion;

    FORTRAN_NAME(projplane)(first_field, second_field,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[7], ProjectedFieldDims+adim, 
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &ProjType, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);

    /* 9) Brightness Temperature */

    first_field = bright_temp;
    //second_field = NULL;
    second_field = BaryonField[0];
    ProjType = 4;
    ConversionFactor = DensityConversion;

    FORTRAN_NAME(projplane)(first_field, second_field,
                             BaryonField[NumberOfBaryonFields], &One,
                             &ProjectionSmooth,
			  GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
                          ProjectedField[8], ProjectedFieldDims+adim, 
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
                          &ProjectionDimension, &ProjType, &ConversionFactor,
                          GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
                          &start, &stop, &LeftCellFraction,&RightCellFraction);
  }

  /* 10-27) Cooling luminosities */

  int offset;

  if (MultiSpecies) {

    this->ComputeLuminosity(all_luminosities, NumberOfLuminosityFields);

    for (n = 0; n < NumberOfLuminosityFields; n++) {
      offset = n*size;
      for (i = 0; i < size; i++)
	temp_field[i] = all_luminosities[i+offset];

      first_field = temp_field;
      second_field = NULL;
      ProjType = 1;
      ConversionFactor = LuminosityConversion;
      //ConversionFactor = One;

      FORTRAN_NAME(projplane)(first_field, second_field,
			      BaryonField[NumberOfBaryonFields], &One,
			      &ProjectionSmooth,
			      GridDimension, GridDimension+1,
                              GridDimension+2, CellWidth[0],
			      ProjectedField[9+n], ProjectedFieldDims+adim, 
                              ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
			      &ProjectionDimension, &ProjType, &ConversionFactor,
			      GridLeftEdge, GridFarLeftEdge, GridRightEdge,
                              ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
			      &start, &stop, &LeftCellFraction,&RightCellFraction);

    } // ENDFOR n
  } // ENDIF MultiSpecies

  /* 28) Balmer emission */

  if (MultiSpecies) {

    // Based on JHW fit to Clegg et al. (1999).  
    // Good between n_elec = [1e2,1e9]
    // Returns H-alpha line emissivity in erg cm^3 s^-1

    float nelec, log_nelec, log_temp, log_nelec2, log_temp2, a0, a1, a2;

    for (i = 0; i < size; i++) {

      if (BaryonField[DeNum][i] < tiny_number) {
	temp_field[i] = 0;
	continue;
      }

      nelec = dom * BaryonField[DeNum][i];
      log_nelec = max(min(log10f(nelec), 9), 2);
      log_nelec2 = log_nelec * log_nelec;
      log_temp = log10f(temperature[i]);
      log_temp2 = log_temp * log_temp;
      a0 = -19.7929  - 1.28249*log_nelec   + 0.198715*log_nelec2;
      a1 = -1.34084  + 0.61503*log_nelec   - 0.097026*log_nelec2;
      a2 = +0.044186 - 0.0736903*log_nelec + 0.0118072*log_nelec2;

      temp_field[i] = powf(10.0f, a0 + a1*log_temp + a2*log_temp2) * nelec*nelec;
      if (isnan(temp_field[i]))
	printf("NaN: %d %g %g %g %g %g %g %g\n", 
	       i, a0, a1, a2, temperature[i], nelec,
	       log_temp, log_nelec);
    }

    first_field = temp_field;
    second_field = NULL;
    ProjType = 1;
    ConversionFactor = LuminosityConversion;
    //ConversionFactor = One;

    FORTRAN_NAME(projplane)(first_field, second_field,
			    BaryonField[NumberOfBaryonFields], &One,
			    &ProjectionSmooth,
			    GridDimension, GridDimension+1,
			    GridDimension+2, CellWidth[0],
			    ProjectedField[27], ProjectedFieldDims+adim, 
			    ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
			    &ProjectionDimension, &ProjType, &ConversionFactor,
			    GridLeftEdge, GridFarLeftEdge, GridRightEdge,
			    ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
			    &start, &stop, &LeftCellFraction,&RightCellFraction);

  } // ENDIF MultiSpecies

  /* Optional: 29-39) Metal line luminosities */

  float max_lum;
  float *all_metal_emis;

  if (MetalLinesUseLookupTable && MultiSpecies) {

    /* Get total metal line luminosity */

    max_lum = -1e20;
    offset = MetalLumField*size;
    for (i = 0; i < size; i++) {
      temp_field[i] = all_luminosities[i+offset];
      max_lum = max(max_lum, temp_field[i]);
    }

    /* Scale the total metal line luminosity to individual lines from
       the lookup table */

    if (max_lum > 0) {
      all_metal_emis = new float[NumberOfMetalLines * size];

      this->ComputeMetalLineLuminosity(temp_field, all_metal_emis, temperature);

      for (n = 0; n < NumberOfMetalLines; n++) {

	offset = n*size;
	for (i = 0; i < size; i++)
	  temp_field[i] = all_metal_emis[i+offset];

	first_field = temp_field;
	second_field = NULL;
	ProjType = 1;
	ConversionFactor = LuminosityConversion;
	//ConversionFactor = One;

	FORTRAN_NAME(projplane)(first_field, second_field,
				BaryonField[NumberOfBaryonFields], &One,
				&ProjectionSmooth,
				GridDimension, GridDimension+1,
				GridDimension+2, CellWidth[0],
				ProjectedField[28+n],
				ProjectedFieldDims+adim, 
				ProjectedFieldDims+bdim, &ProjectedFieldCellSize,
				&ProjectionDimension, &ProjType, &ConversionFactor,
				GridLeftEdge, GridFarLeftEdge, GridRightEdge,
				ProjectedFieldLeftEdge, ProjectedFieldRightEdge,
				&start, &stop, &LeftCellFraction,&RightCellFraction);

      } // ENDFOR metal lines

      delete [] all_metal_emis;
    } // ENDIF max_lum > 0
  } // ENDIF Metal lines

 StarParticleLabel:

  /* Clean up. */

  delete [] temperature;
  delete [] temp_field;
  delete [] spin_temp;
  delete [] bright_temp;
  delete [] all_luminosities;

  return SUCCESS;
}
