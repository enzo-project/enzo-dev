/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COLLAPSE TEST WITH PHOTONS)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Oct, 2003 by Tom Abel -- include photons
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

/********************* PROTOTYPES *********************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
float ph_gasdev();

// Used to compute Bonner-Ebert density profile
double ph_BE(double r);
double ph_q(double r);
double ph_Ang(double a1, double a2, double R, double r);

// Returns random velocity from Maxwellian distribution
double ph_Maxwellian(double c_tilda, double vel_unit, double mu, double gamma);

/*******************************************************/

static int PhotonTestParticleCount = 0;

int grid::PhotonTestInitializeGrid(int NumberOfSpheres,
			     float SphereRadius[MAX_SPHERES],
			     float SphereCoreRadius[MAX_SPHERES],
			     float SphereDensity[MAX_SPHERES],
			     float SphereTemperature[MAX_SPHERES],
			     FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
			     float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
                             float SphereFracKeplerianRot[MAX_SPHERES],
                             float SphereTurbulence[MAX_SPHERES],
                             float SphereCutOff[MAX_SPHERES],
                             float SphereAng1[MAX_SPHERES],
                             float SphereAng2[MAX_SPHERES],
                             int   SphereNumShells[MAX_SPHERES],
			     int   SphereType[MAX_SPHERES],
			     int   SphereConstantPressure[MAX_SPHERES],
			     int   SphereSmoothSurface[MAX_SPHERES],
			     float SphereSmoothRadius[MAX_SPHERES],
			     float SphereHII[MAX_SPHERES],
			     float SphereHeII[MAX_SPHERES],
			     float SphereHeIII[MAX_SPHERES],
			     float SphereH2I[MAX_SPHERES],
			     int   SphereUseParticles,
			     float UniformVelocity[MAX_DIMENSION],
			     int   SphereUseColour,
			     float InitialTemperature, int level, 
			     float PhotonTestInitialFractionHII, 
			     float PhotonTestInitialFractionHeII,
			     float PhotonTestInitialFractionHeIII, 
			     float PhotonTestInitialFractionHM,
			     float PhotonTestInitialFractionH2I, 
			     float PhotonTestInitialFractionH2II,
			     int RefineByOpticalDepth,
			     int TotalRefinement,
			     char *DensityFilename,
			     char *HIIFractionFilename,
			     char *HeIIFractionFilename,
			     char *HeIIIFractionFilename,
			     char *TemperatureFilename)
{
  /* declarations */

  int dim, i, j, k, m, field, sphere, size, active_size, index, cindex;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum,  kphHINum, gammaNum, kphHeINum,
    kphHeIINum, kdissH2INum, RPresNum1, RPresNum2, RPresNum3; 
  float *density_field = NULL, *HII_field = NULL, *HeII_field = NULL, 
    *HeIII_field = NULL, *Temperature_field = NULL;

  /* create fields */
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  int ivel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }
  int ColourNum = NumberOfBaryonFields;
  if (SphereUseColour)
    FieldType[NumberOfBaryonFields++] = Metallicity; /* fake it with metals */

  if (RadiativeTransfer && (MultiSpecies < 1)) {
    ENZO_FAIL("Grid_PhotonTestInitialize: Radiative Transfer but not MultiSpecies set");
  }

  //   Allocate fields for photo ionization and heating rates
  if (RadiativeTransfer)
    if (MultiSpecies) {
      FieldType[kphHINum    = NumberOfBaryonFields++] = kphHI;
      FieldType[gammaNum    = NumberOfBaryonFields++] = PhotoGamma;
      if (RadiativeTransferHydrogenOnly == FALSE) {
	FieldType[kphHeINum   = NumberOfBaryonFields++] = kphHeI;
	FieldType[kphHeIINum  = NumberOfBaryonFields++] = kphHeII;
      }
      if (MultiSpecies > 1) 
	FieldType[kdissH2INum    = NumberOfBaryonFields++] = kdissH2I;
    } 

  if (RadiationPressure && RadiativeTransfer) {
    FieldType[RPresNum1 = NumberOfBaryonFields++] = RadPressure0;
    FieldType[RPresNum2 = NumberOfBaryonFields++] = RadPressure1;
    FieldType[RPresNum3 = NumberOfBaryonFields++] = RadPressure2;
  }

  NumberOfPhotonPackages = 0;
  PhotonPackages-> NextPackage= NULL;

  /* Return if this doesn't concern us. */

  int ReadData = (ParallelRootGridIO == FALSE ||
		  (ParallelRootGridIO == TRUE && TotalRefinement < 0));

  if (ProcessorNumber != MyProcessorNumber || !ReadData) {
    NumberOfParticles = (SphereUseParticles > 0) ? 1 : 0;
    for (dim = 0; dim < GridRank; dim++)
      NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
  }

  /* Set various units. */

  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.67e-8,
               pi = 3.14159, mh = 1.67e-24, kboltz = 1.381e-16;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, mu = 0.6, mu_data;

  FLOAT a, dadt, ExpansionFactor = 1;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    ExpansionFactor = a/(1.0+InitialRedshift);
    CriticalDensity = 2.78e11*pow(HubbleConstantNow, 2); // in Msolar/Mpc^3
    BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
  } else {
    BoxLength = LengthUnits / 3.086e24;
  }

  /* Compute NFW profile-info. The parameters are SphereCoreRadius (which is
     the "knee radius") and SphereDensity (the overdensity ar the knee
     radius).  The pressure is computed by integrated the hydro-static
     equation and from this comes the temperature and the dm velocity
     dispersion. */

#define NFW_POINTS 500
  float NFWMass[NFW_POINTS], NFWPressure[NFW_POINTS], NFWTemp[NFW_POINTS], x1,
        NFWDensity[NFW_POINTS], NFWSigma[NFW_POINTS], m200;
  FLOAT NFWRadius[NFW_POINTS];
  double dpdr = 0, dpdr_old;
  sphere = 0;
  m200   = 0;
  NFWPressure[0] = 1.0 * kboltz * InitialTemperature / (mu * mh);
  FILE *fptr = fopen("NFWProfile.out", "w");
  for (i = 0; i < NFW_POINTS; i++) {
    NFWRadius[i] = SphereRadius[sphere]*pow(10, -3*(float(i)/NFW_POINTS));
    x1 = NFWRadius[i]/SphereCoreRadius[sphere];
    NFWDensity[i] = SphereDensity[sphere]/(x1*(1.0+x1)*(1.0+x1));
    NFWMass[i] = 4.0*pi*SphereDensity[sphere]*
                      (CriticalDensity/pow(ExpansionFactor, 3)) *
		pow(SphereCoreRadius[sphere]*BoxLength, 3) *
		(log(1.0+x1) - x1/(x1+1.0));  // in Msolar
    dpdr_old = dpdr;
    dpdr = GravConst * NFWMass[i] * SolarMass * 
           NFWDensity[i] / 
           pow(NFWRadius[i]*BoxLength*Mpc, 2);
    if (i > 0)
      NFWPressure[i] = NFWPressure[i-1] -
	0.5*(dpdr+dpdr_old)*(NFWRadius[i]-NFWRadius[i-1])*BoxLength*Mpc;
    NFWTemp[i] = NFWPressure[i]*mu*mh/(kboltz*NFWDensity[i]); // in K
    NFWSigma[i] = sqrt(kboltz * NFWTemp[i] / (mu * mh));  // in cm/s
    float mean_overdensity = 3.0*SphereDensity[sphere] / (x1*x1*x1) *
        (log(1.0+x1) - x1/(x1+1.0));
    fprintf(fptr, "%"ISYM" %"GOUTSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", i, NFWRadius[i], 
	 NFWDensity[i], NFWMass[i], NFWPressure[i], NFWTemp[i], NFWSigma[i],
         mean_overdensity);
    if (mean_overdensity > 200 && m200 == 0)
      m200 = NFWMass[i];
  }
  fprintf(fptr, "#m200 = %"GSYM"\n", m200);
  fclose(fptr);

  /* Loop over the set-up twice, once to count the particles, the second
     time to initialize them. */

  int SetupLoopCount, npart = 0;
  for (SetupLoopCount = 0; SetupLoopCount < 1+min(SphereUseParticles, 1);
       SetupLoopCount++) {

  /* Set densities */

  float BaryonMeanDensity = SphereUseParticles ? 0.174 : 1.0;
  if (SphereUseParticles == 2) BaryonMeanDensity = 0.826;
  float ParticleMeanDensity = 1.0 - BaryonMeanDensity, ParticleCount = 0;

  /* Set particles. */

  if (SphereUseParticles > 0 && SetupLoopCount > 0) {
    /* If particles already exist (coarse particles), then delete. */
    if (NumberOfParticles > 0)
      this->DeleteParticles();
    /* Use count from previous loop to set particle number. */
    NumberOfParticles = npart;
    npart = 0;
    /* Allocate space. */
    this->AllocateNewParticles(NumberOfParticles);
    /* Particle values will be set below. */
  } // end: particle initialization

  /* Set up the baryon field. */
  /* compute size of fields */
  size = 1;
  active_size = 1;
  int ActiveDims[MAX_DIMENSION];
  for (dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
    ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    active_size *= ActiveDims[dim];
  }

  /* allocate fields */
  if (SetupLoopCount == 0)
    {
      for (field = 0; field < NumberOfBaryonFields; field++)
	if (BaryonField[field] == NULL)
	  BaryonField[field] = new float[size];
    }

  /* Initialize radiation fields */

  if (this->InitializeRadiativeTransferFields() == FAIL) {
    ENZO_FAIL("\nError in InitializeRadiativeTransferFields.\n");
  }

  /* Read density field, if given */

  if (DensityFilename != NULL) {
    char *data_filename, *dataset_name;
    hsize_t OutDims[MAX_DIMENSION];
    herr_t h5error;
    hid_t file_id;

    // Parse DensityFilename into filename and dataset name
    char *delim = "/";
    data_filename = strtok(DensityFilename, delim);
    dataset_name = strtok(NULL, delim);

    file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id == -1) ENZO_FAIL("Error opening density field.");
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      OutDims[GridRank-dim-1] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    density_field = new float[active_size];
    this->read_dataset(GridRank, OutDims, dataset_name, file_id,
		       HDF5_REAL, density_field, FALSE, NULL, NULL);
    h5error = H5Fclose(file_id);
    if (h5error == -1) ENZO_FAIL("Error closing density field.");

    /* HII, HeI, HeII, HeIII density fields */

    if (HIIFractionFilename != NULL) {
      data_filename = strtok(HIIFractionFilename, delim);
      dataset_name = strtok(NULL, delim);
      file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file_id == -1) ENZO_FAIL("Error opening HII field.");
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	OutDims[GridRank-dim-1] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
      HII_field = new float[active_size];
      this->read_dataset(GridRank, OutDims, dataset_name, file_id,
			 HDF5_REAL, HII_field, FALSE, NULL, NULL);
      h5error = H5Fclose(file_id);
      if (h5error == -1) ENZO_FAIL("Error closing HII field.");
    }

    if (HeIIFractionFilename != NULL) {
      data_filename = strtok(HeIIFractionFilename, delim);
      dataset_name = strtok(NULL, delim);
      file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file_id == -1) ENZO_FAIL("Error opening HeII field.");
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	OutDims[GridRank-dim-1] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
      HeII_field = new float[active_size];
      this->read_dataset(GridRank, OutDims, dataset_name, file_id,
			 HDF5_REAL, HeII_field, FALSE, NULL, NULL);
      h5error = H5Fclose(file_id);
      if (h5error == -1) ENZO_FAIL("Error closing HeII field.");
    }

    if (HeIIIFractionFilename != NULL) {
      data_filename = strtok(HeIIIFractionFilename, delim);
      dataset_name = strtok(NULL, delim);
      file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file_id == -1) ENZO_FAIL("Error opening HeIII field.");
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	OutDims[GridRank-dim-1] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
      HeIII_field = new float[active_size];
      this->read_dataset(GridRank, OutDims, dataset_name, file_id,
			 HDF5_REAL, HeIII_field, FALSE, NULL, NULL);
      h5error = H5Fclose(file_id);
      if (h5error == -1) ENZO_FAIL("Error closing HeIII field.");
    }

    if (TemperatureFilename != NULL) {
      data_filename = strtok(TemperatureFilename, delim);
      dataset_name = strtok(NULL, delim);
      file_id = H5Fopen(data_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if (file_id == -1) ENZO_FAIL("Error opening temperature field.");
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	OutDims[GridRank-dim-1] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
      Temperature_field = new float[active_size];
      this->read_dataset(GridRank, OutDims, dataset_name, file_id,
			 HDF5_REAL, Temperature_field, FALSE, NULL, NULL);
      h5error = H5Fclose(file_id);
      if (h5error == -1) ENZO_FAIL("Error closing temperature field.");
    }


  } // ENDIF DensityFilename

  /* Loop over the mesh. */
  float density, dens1, Velocity[MAX_DIMENSION],
    temperature, temp1, sigma, sigma1, colour, outer_radius;
  float HII_Fraction, HeII_Fraction, HeIII_Fraction, H2I_Fraction;
  FLOAT r, x, y = 0, z = 0;
  int n = 0;

  float SphereRotationalPeriod[MAX_SPHERES];
  float HydrostaticTemperature[MAX_SPHERES];
  float VelocityKep = 0;
  float VelocitySound[MAX_SPHERES];
  double SphereMass, SphereCoreMass, SphereCoreDens;
  float alpha, beta, theta;
  float Scale_Factor[MAX_SPHERES];

  /* Pre-compute cloud properties before looping over mesh */

  for (sphere = 0; sphere < NumberOfSpheres; sphere++) {
    Scale_Factor[sphere] = SphereCutOff[sphere] / SphereRadius[sphere];
    HydrostaticTemperature[sphere] = (2*M_PI * GravConst * mh) / 
      (3.0*kboltz) * (SphereDensity[sphere] * DensityUnits) * 
      pow(SphereRadius[sphere] * LengthUnits, 2.0);
    if (SphereFracKeplerianRot[sphere] > 0.0) {
//      if (SphereType[sphere] == 7)
//	HydrostaticTemperature[sphere] *= 1.0 - SphereFracKeplerianRot[sphere];
      if (SphereType[sphere] == 5) {
	SphereCoreDens = (SphereDensity[sphere]*DensityUnits) * 
	  pow(SphereCoreRadius[sphere] / SphereRadius[sphere], -2);
	SphereCoreMass = (4*pi/3) * 
	  pow((SphereCoreRadius[sphere]*LengthUnits),3) * (SphereCoreDens);
	SphereMass = (4*pi*SphereDensity[sphere]*DensityUnits * 
		      pow(SphereRadius[sphere]*LengthUnits, 2) *
		      ((SphereRadius[sphere] - SphereCoreRadius[sphere]) * 
		       LengthUnits)) + SphereCoreMass;
	printf("\nSphere Mass (M_sun): %"FSYM"\n", SphereMass/SolarMass);
      }
      else if (SphereType[sphere] == 1 || SphereType[sphere] == 7) {
	SphereMass = double(4*pi/3) *
	  pow((SphereRadius[sphere]*LengthUnits), 3) *
	  double(SphereDensity[sphere]*DensityUnits);
	printf("\nSphere Mass (M_sun): %"FSYM"\n", SphereMass/SolarMass);
      } 
      else if (SphereType[sphere] == 6) {
	VelocitySound[sphere] = sqrt((SphereTemperature[sphere] * Gamma)/mu);
	SphereMass = pow(VelocitySound[sphere],3) / 
	  (sqrt(4*pi*pow((GravConst)/(4*pi),3) * 
		(SphereDensity[sphere]))) * 
	  ph_q(SphereCutOff[sphere]);
	SphereMass = SphereMass*DensityUnits*pow(LengthUnits,3);
	printf("\nSphere Mass (M_sun): %"FSYM"\n", SphereMass/SolarMass);
      }
      
      VelocityKep = sqrt(6.673e-8*SphereMass/(SphereRadius[sphere]*(LengthUnits)));
      SphereRotationalPeriod[sphere] = 2*pi*SphereRadius[sphere]*(LengthUnits)/
	(SphereFracKeplerianRot[sphere]*VelocityKep);
      SphereRotationalPeriod[sphere] = SphereRotationalPeriod[sphere]/(TimeUnits);

      printf("\nKeplerian Rotation Period (s): %"FSYM"\n", SphereRotationalPeriod[sphere] 
	     * SphereFracKeplerianRot[sphere]*TimeUnits);
      printf("\nSphere Rotation Period (s): %"FSYM"\n", SphereRotationalPeriod[sphere]
	     * TimeUnits);
    } else
      SphereRotationalPeriod[sphere] = 0.0;



    // Calculate speed of sound for this sphere
    VelocitySound[sphere] = sqrt((SphereTemperature[sphere] * Gamma)/mu);
    printf("\nVelocitySound (cm s^-1): %"FSYM"\n", VelocitySound[sphere] * 
	   (LengthUnits/TimeUnits));

  } // ENDFOR sphere

  /* Loop over mesh */

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++, n++) {

	/* Compute position */

	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	/* Loop over spheres. */

	if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
	    j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
	    k >= GridStartIndex[2] && k <= GridEndIndex[2]) {
	  cindex = (i-GridStartIndex[0]) + ActiveDims[0] *
	    ((j-GridStartIndex[1]) + (k-GridStartIndex[2])*ActiveDims[1]);
	  if (density_field != NULL)
	    density = max(density_field[cindex], 1e-6);
	  else
	    density = 1.0;
	  if (HII_field != NULL)
	    HII_Fraction = HII_field[cindex];
	  else
	    HII_Fraction = PhotonTestInitialFractionHII;
	  if (HeII_field != NULL)
	    HeII_Fraction = HeII_field[cindex];
	  else
	    HeII_Fraction = PhotonTestInitialFractionHeII;
	  if (HeIII_field != NULL)
	    HeIII_Fraction = HeIII_field[cindex];
	  else
	    HeIII_Fraction = PhotonTestInitialFractionHeIII;
	  if (Temperature_field != NULL)
	    temperature = temp1 = Temperature_field[cindex];
	  else
	    temperature = temp1 = InitialTemperature;
	}

	H2I_Fraction = PhotonTestInitialFractionH2I;
	sigma = sigma1 = 0;
	colour = 1.0e-10;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Velocity[dim] = 0;

	for (sphere = 0; sphere < NumberOfSpheres; sphere++) {

	  /* Find distance from center. */

	  r = sqrt(pow(fabs(x-SpherePosition[sphere][0]), 2) +
		   pow(fabs(y-SpherePosition[sphere][1]), 2) +
		   pow(fabs(z-SpherePosition[sphere][2]), 2) );
	  r = max(r, 0.1*CellWidth[0][0]);

	  outer_radius = (SphereSmoothSurface[sphere] == TRUE) ? 
	    SphereSmoothRadius[sphere]*SphereRadius[sphere] : SphereRadius[sphere];
	  if (r < outer_radius) {

	    /* Compute Cartesian coordinates for rotational properties */

            FLOAT xpos, ypos, zpos, drad;

	    xpos = x-SpherePosition[sphere][0] - (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
	    ypos = y-SpherePosition[sphere][1] - (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
	    zpos = z-SpherePosition[sphere][2] - (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);
	    drad = sqrt(xpos*xpos + ypos*ypos + zpos*zpos);
	    alpha = 2*pi*(FLOAT(rand())/FLOAT(RAND_MAX));
	    beta = (pi/2)*(((2*FLOAT(rand()))/FLOAT(RAND_MAX)) - 1);

	    /* Compute spherical coordinate theta */

	    if (xpos != 0) {
	      if (xpos > 0 && ypos >= 0)
		theta = atan(ypos/xpos);
	      else if (xpos < 0 && ypos >= 0)
		theta = pi + atan(ypos/xpos);
	      else if (xpos < 0 && ypos < 0)
		theta = pi + atan(ypos/xpos);
	      else if (xpos > 0 && ypos < 0)
		theta = 2*pi + atan(ypos/xpos);
             } else if (xpos == 0 && ypos > 0)
	      theta = 3.14159 / 2.0;
	    else if (xpos == 0 && ypos < 0)
	      theta = (3*3.14159) / 2.0;
	    else
	      theta = 0.0;

	    // Find out which shell the cell is in
	    a = ph_Ang(SphereAng1[sphere],SphereAng2[sphere],SphereRadius[sphere],r);

	    /* Start with solid body rotation and then add in a
	       velocity of a fraction of sound speed in a random
	       direction to create turbulence */
	    
	    if (SphereRotationalPeriod[sphere] > 0) {
	      Velocity[0] = -2*pi*(ypos*cos(a) + zpos*sin(a)) / 
		SphereRotationalPeriod[sphere];
	      Velocity[1] = 2*pi*(xpos*cos(a)) / SphereRotationalPeriod[sphere];
	      Velocity[2] = 2*pi*(xpos*sin(a)) / SphereRotationalPeriod[sphere];
	    }
	    Velocity[0] += SphereTurbulence[sphere] * 
	      ph_Maxwellian(VelocitySound[sphere], LengthUnits/TimeUnits, mu, Gamma);
	    Velocity[1] += SphereTurbulence[sphere] * 
	      ph_Maxwellian(VelocitySound[sphere], LengthUnits/TimeUnits, mu, Gamma);
             if (dim == 0 || dim == 3)
	       Velocity[2] += SphereTurbulence[sphere] * 
		 ph_Maxwellian(VelocitySound[sphere], LengthUnits/TimeUnits, mu, Gamma);
	     m = 0;

	    /* 1) Uniform */

	    if (SphereType[sphere] == 1)
	      dens1 = SphereDensity[sphere];

	    /* 2) r^-2 power law */

	    if (SphereType[sphere] == 2)
	      dens1 = SphereDensity[sphere]*pow(r/SphereRadius[sphere], -2);

	    /* 3) NFW profile (use look-up table for temperature and
	          velocity dispersion)*/

	    if (SphereType[sphere] == 3) {
	      x1 = r/SphereCoreRadius[sphere];
	      dens1 = SphereDensity[sphere]/(x1*(1.0+x1)*(1.0+x1));
	      for (m = 1; m < NFW_POINTS; m++)
		if (r > NFWRadius[m]) {
		  temp1 = NFWTemp[m] + (NFWTemp[m-1] - NFWTemp[m])*
		    (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
		  sigma1 = NFWSigma[m] + (NFWSigma[m-1] - NFWSigma[m])*
		    (r - NFWRadius[m])/(NFWRadius[m-1] - NFWRadius[m]);
		  break;
		}
	    }

	    /* 4) Gaussian */
	    if (SphereType[sphere] == 4) {
	      dens1 = SphereDensity[sphere]*
                      PEXP(-0.5*pow(r/SphereCoreRadius[sphere], 2));
	    }

	    /* 5) r^-2 power law with core radius */

	    if (SphereType[sphere] == 5) {
	      if (r < SphereCoreRadius[sphere]) {
		dens1 = SphereDensity[sphere]*pow(SphereCoreRadius[sphere]/
						  SphereRadius[sphere], -2);
                //dens1 = dens1*(1.1 + 0.1*cos(2*theta));
	      } else {
		dens1 = SphereDensity[sphere]*pow(r/SphereRadius[sphere], -2);
                //dens1 = dens1*(1.1 + 0.1*cos(2*theta));
	      }
	    }

	    /* 6) Free expansion blastwave (see FreeExpansionInitialize) */

	    if (SphereType[sphere] == 6) {

	      const float DensitySlope = 9.0;
	      double M_ej, E_ej, r_max, v_max, BlastTime, v_core, normalization,
		speed;

	      M_ej = SphereDensity[sphere] * (4.0*M_PI/3.0) * 
		POW(SphereRadius[sphere]*LengthUnits, 3.0) * DensityUnits;
	      E_ej = M_ej * kboltz * SphereTemperature[sphere] / (mh*mu);
	      r_max = LengthUnits * SphereRadius[sphere];

	      // Temperature parameter is a proxy for total energy
	      // (convert to velocity)
	      v_max = 0.333333 * sqrt(2.0 * SphereTemperature[sphere] / 
				      TemperatureUnits / ((Gamma-1.0)*mu));
	      BlastTime = SphereRadius[sphere] / v_max;
	      v_core = sqrt( (10.0 * E_ej * (DensitySlope-5)) /
			     (3.0 * M_ej * (DensitySlope-3)) );
	      normalization = (10.0 * (DensitySlope-5) * E_ej) / 
		(4.0 * M_PI * DensitySlope) / POW(v_core, 5.0);
	      printf("v_max = %g\n", v_max * VelocityUnits * 1e-5);
	      
	      v_core /= VelocityUnits;
	      speed = r / BlastTime;
	      if (speed <= v_core)
		dens1 = normalization / POW(BlastTime*TimeUnits, 3.0) / DensityUnits;
	      else
		dens1 = normalization / POW(BlastTime*TimeUnits, 3.0) /
		  POW(speed/v_core, DensitySlope) / DensityUnits;
	      dens1 = max(dens1, SphereDensity[sphere]);
	      Velocity[0] = speed * xpos / r;
	      Velocity[1] = speed * ypos / r;
	      Velocity[2] = speed * zpos / r;
	      temp1 = 10.0;

	    } // ENDIF SphereType 6

	    /* 7) Uniform density in hydrostatic equilbrium */

	    if (SphereType[sphere] == 7) {
	      dens1 = SphereDensity[sphere];
	      temp1 = HydrostaticTemperature[sphere] *
		pow(r / SphereRadius[sphere], 2.0);
	      temp1 = max(temp1, 1.0);
	    } // ENDIF type 7

	    /* 10) disk (ok, it's not a sphere, so shoot me) */

	    if (SphereType[sphere] == 10) {

	      FLOAT xpos, ypos, zpos, xpos1, ypos1, zpos1, zheight, drad;
	      FLOAT ScaleHeightz = SphereCoreRadius[sphere]/6.0,
	            ScaleHeightR = SphereCoreRadius[sphere];

	      /* Loop over dims if using Zeus (since vel's face-centered). */

	      for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
		   dim++) {

		/* Compute position. */

		xpos = x-SpherePosition[sphere][0] - 
		  (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
		ypos = y-SpherePosition[sphere][1] -
		  (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
		zpos = z-SpherePosition[sphere][2] -
		  (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);

		/* Compute z and r_perp (SphereVelocity is angular momentum 
		   and must have unit length). */    

		zheight = SphereVelocity[sphere][0]*xpos + 
		          SphereVelocity[sphere][1]*ypos +
		          SphereVelocity[sphere][2]*zpos;
		xpos1 = xpos - zheight*SphereVelocity[sphere][0];
		ypos1 = ypos - zheight*SphereVelocity[sphere][1];
		zpos1 = zpos - zheight*SphereVelocity[sphere][2];
		drad = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);

		/* If we're above the disk, then exit. */

//		if (zheight > max(5.0*ScaleHeightz, 2.0*CellWidth[0][0]))
//		  continue;

		/* Compute density (Kruit & Searle 1982). */

		if (dim == 0)
		  dens1 = SphereDensity[sphere]*PEXP(-drad/ScaleHeightR)/
		    pow(cosh(zheight/max(ScaleHeightz, CellWidth[0][0])), 2);

		if (dens1 < density)
		  break;

		/* Compute velocity magnitude (divided by drad). 
		   This assumes PointSourceGravityPosition and Sphere center 
		   are the same.  This should be fixed to use the disk mass
		   as well, but that's a bit tricky. */

//		float vel = sqrt(PointSourceGravityConstant/drad)/drad;

		float accel = 0;
		if (PointSourceGravity == 1)
		  accel = PointSourceGravityConstant/
		    (pow(drad,3) + pow(PointSourceGravityCoreRadius, 3));
		if (PointSourceGravity == 2) {
		  x1 = drad/PointSourceGravityCoreRadius;
		  accel = PointSourceGravityConstant*(log(1+x1)-x1/(1+x1))/
		           pow(drad, 3);
		}
		
		float vel = sqrt(accel);
		
		/* Compute velocty: L x r_perp. */

		if (dim == 0 || dim == 1)
		  Velocity[0] = vel*(SphereVelocity[sphere][1]*zpos1 -
				     SphereVelocity[sphere][2]*ypos1);
		if (dim == 0 || dim == 2)
		  Velocity[1] = vel*(SphereVelocity[sphere][2]*xpos1 +
				     SphereVelocity[sphere][0]*zpos1);
		if (dim == 0 || dim == 3)
		  Velocity[2] = vel*(SphereVelocity[sphere][0]*ypos1 -
				     SphereVelocity[sphere][1]*xpos1);

	      } // end: loop over dims

	    } // end: disk
	    
	    /* If the density is larger than the background (or the previous
	       sphere), then set the velocity. */

	    if (dens1 > density) {
	      density = dens1;
	      if (SphereType[sphere] != 7 && SphereType[sphere] != 9) {
            if (temp1 == InitialTemperature) {
              if (SphereConstantPressure[sphere] == TRUE) {
                temperature = SphereTemperature[sphere] * 
                  (SphereDensity[sphere] / dens1);
              } else {
                temperature = SphereTemperature[sphere];
              }
            } else {
              temperature = temp1;
            }
          }

	      sigma = sigma1;
	      if (SphereType[sphere] != 6 &&
		  SphereType[sphere] != 10 &&
		  SphereRotationalPeriod[sphere] <= 0)
		for (dim = 0; dim < GridRank; dim++)
		  Velocity[dim] = SphereVelocity[sphere][dim];
	      if (sphere == 0)
		colour = dens1; /* only mark first sphere */
	      HII_Fraction = SphereHII[sphere];
	      HeII_Fraction = SphereHeII[sphere];
	      HeIII_Fraction = SphereHeIII[sphere];
	      H2I_Fraction = SphereH2I[sphere];
	    }

	  } // end: if (r < SphereRadius)

	  if (SphereSmoothSurface[sphere] == TRUE && 
	      r < SphereSmoothRadius[sphere]*SphereRadius[sphere] &&
		  r > SphereRadius[sphere]) {
	    float ramp = 1.0 - 1.0 * tanh((3.0/(SphereSmoothRadius[sphere]-1.0))*
					  (r/SphereRadius[sphere] - 1.0));
	    ramp = max(ramp, 1.0/density);
	    density *= ramp;
	    if (SphereConstantPressure[sphere] == TRUE) {
	      temperature /= ramp;
	    }
	  } // end: if (SmoothSurface)
	  
	} // end: loop over spheres

	/* Set density. */

	BaryonField[0][n] = density*BaryonMeanDensity;

	/* If doing multi-species (HI, etc.), set these. */

	if (MultiSpecies > 0) {
	  BaryonField[HIINum][n] = HII_Fraction *
	    CoolData.HydrogenFractionByMass * BaryonField[0][n];
	  BaryonField[HeIINum][n] = HeII_Fraction *
	    BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	  BaryonField[HeIIINum][n] = HeIII_Fraction *
	    BaryonField[0][n] * 4.0 * (1.0-CoolData.HydrogenFractionByMass);
	  BaryonField[HeINum][n] = 
	    (1.0 - CoolData.HydrogenFractionByMass)*BaryonField[0][n] -
	    BaryonField[HeIINum][n] - BaryonField[HeIIINum][n];
	}
	if (MultiSpecies > 1) {
	  BaryonField[HMNum][n] = PhotonTestInitialFractionHM*
	    BaryonField[HIINum][n]* pow(temperature,float(0.88));
	  BaryonField[H2IINum][n] = PhotonTestInitialFractionH2II*
	    2.0*BaryonField[HIINum][n]* pow(temperature,float(1.8));
	  if (ComovingCoordinates)
	    BaryonField[H2INum][n] = H2I_Fraction *
	      BaryonField[0][n]*CoolData.HydrogenFractionByMass*pow(301.0,5.1)*
	      pow(OmegaMatterNow, float(1.5))/
	      (OmegaMatterNow*BaryonMeanDensity)/
	      HubbleConstantNow*2.0;
	  else
	    BaryonField[H2INum][n] = H2I_Fraction *
	      BaryonField[0][n]*CoolData.HydrogenFractionByMass*pow(301.0,5.1)/
	      (2.0*BaryonMeanDensity);
	}
	
	BaryonField[HINum][n] = 
	  CoolData.HydrogenFractionByMass*BaryonField[0][n]
	  - BaryonField[HIINum][n];
	if (MultiSpecies > 1)
	  BaryonField[HINum][n] -= BaryonField[HMNum][n]
	    + BaryonField[H2IINum][n]
	    + BaryonField[H2INum][n];
	
	BaryonField[DeNum][n] = BaryonField[HIINum][n] + 
	  0.25*BaryonField[HeIINum][n] + 0.5*BaryonField[HeIIINum][n];
	if (MultiSpecies > 1)
	  BaryonField[DeNum][n] += 0.5*BaryonField[H2IINum][n] - 
	    BaryonField[HMNum][n];
	
	/* Set Deuterium species (assumed to be negligible). */
	
	if (MultiSpecies > 2) {
	  BaryonField[DINum][n] = CoolData.DeuteriumToHydrogenRatio*
	    BaryonField[HINum][n];
	  BaryonField[DIINum][n] = CoolData.DeuteriumToHydrogenRatio*
	    BaryonField[HIINum][n];
	  BaryonField[HDINum][n] = CoolData.DeuteriumToHydrogenRatio*
	    BaryonField[H2INum][n];
	}
	
	
	/* If there is a colour field, set it. */
	
	if (SphereUseColour)
	  BaryonField[ColourNum][n] = colour;
	
	/* Set Velocities. */
	
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[ivel+dim][n] = Velocity[dim] + UniformVelocity[dim];
	
	/* Set energy (thermal and then total if necessary). */

	if (MultiSpecies) {
	  mu_data =  
	    0.25*(BaryonField[HeINum][n]  + BaryonField[HeIINum][n] +
		  BaryonField[HeIIINum][n]                        ) +
	    BaryonField[HINum][n]   + BaryonField[HIINum][n]  +
	    BaryonField[DeNum][n];
	  if (MultiSpecies > 1)
	    mu_data += BaryonField[HMNum][i]   +
	      0.5*(BaryonField[H2INum][i]  + BaryonField[H2IINum][i]);
	  mu_data = BaryonField[0][n] / mu_data;
	} else
	  mu_data = mu;

	BaryonField[1][n] = temperature/TemperatureUnits/
	  ((Gamma-1.0)*mu_data);
	
	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n];
	
	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*pow(BaryonField[ivel+dim][n], 2);

	/* Set particles if being used (generate a number of particle
	   proportional to density). */

	if (SphereUseParticles)
	  if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
	      j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
	      k >= GridStartIndex[2] && k <= GridEndIndex[2]  ) {
	    ParticleCount += density/pow(float(RefineBy), GridRank*level);
	    while (ParticleCount > 1) {
	      if (SetupLoopCount > 0) {
		ParticleMass[npart] = ParticleMeanDensity*
		                      pow(float(RefineBy), GridRank*level);
		ParticleNumber[npart] = PhotonTestParticleCount++;
		ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;

		/* Set random position within cell. */

		ParticlePosition[0][npart] = x + 
		        CellWidth[0][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		ParticlePosition[1][npart] = y +
		        CellWidth[1][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		ParticlePosition[2][npart] = z +
		        CellWidth[2][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);

		/* Set bulk velocity. */

		for (dim = 0; dim < GridRank; dim++)
		  ParticleVelocity[dim][npart] = 
		    Velocity[dim]+UniformVelocity[dim];

		/* Add random velocity; */

		if (sigma != 0)
		  for (dim = 0; dim < GridRank; dim++)
		    ParticleVelocity[dim][npart] += 
		      ph_gasdev()*sigma/VelocityUnits;

	      }
	      npart++;
	      ParticleCount -= 1.0;
	    }
	  } // end: if statement

	} // end loop over grid
	
  } // end loop SetupLoopCount

  /* If needed, set a radiation field in the cell where the
     sources resides to flag the cells by optical depth */

  FLOAT pos[3];
  int ci[3];
  RadiationSourceEntry *Source = GlobalRadiationSources->NextSource;

  if (RefineByOpticalDepth) {

    while (Source != NULL) {
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	pos[dim] = Source->Position[dim];
	ci[dim] = int((pos[dim] - GridLeftEdge[dim]) / CellWidth[dim][0]);
      }

      printf("PhotonTest[%"ISYM"]: Left edge = %"FSYM" %"FSYM" %"FSYM"\n", level, GridLeftEdge[0],
	     GridLeftEdge[1], GridLeftEdge[2]);
      printf("PhotonTest[%"ISYM"]: source (%"FSYM" %"FSYM" %"FSYM") in %"ISYM" %"ISYM" %"ISYM"\n", 
	     level, pos[0], pos[1], pos[2], ci[0], ci[1], ci[2]);

      if (pos[0] >= GridLeftEdge[0] && pos[0] <= GridRightEdge[0] &&
	  pos[1] >= GridLeftEdge[1] && pos[1] <= GridRightEdge[1] &&
	  pos[2] >= GridLeftEdge[2] && pos[2] <= GridRightEdge[2]) {

	index = GRIDINDEX(ci[0], ci[1], ci[2]);
	BaryonField[kphHINum][index] = 1;

	printf("PhotonTest[%"ISYM"]: set kphHI in %"ISYM" %"ISYM" %"ISYM" (%"ISYM")\n", 
	       level, ci[0], ci[1], ci[2], index);

      } /* ENDIF inside grid */

      Source = Source->NextSource;

    } /* ENDWHILE sources */

  } /* ENDIF RefineByOpticalDepth */

  if (SphereUseParticles && debug)
    printf("PhotonTestInitialize: DM NumberOfParticles = %"ISYM"\n", 
	   NumberOfParticles);

  delete [] density_field;
  delete [] HII_field;
  delete [] HeII_field;
  delete [] HeIII_field;
  delete [] Temperature_field;
  
  return SUCCESS;
}

/************************************************************************/
/* Routine to return a Gaussian random deviate with zero mean and unite
   variance (adapted from Numerical Recipes). */

static int ph_gasdev_iset = 0;
static float ph_gasdev_gset;

float ph_gasdev()
{
  float v1, v2, r = 0, fac, ph_gasdev_ret;
  if (ph_gasdev_iset == 0) {

    while (r >= 1 || r == 0) {
      v1 = 2.0*float(rand())/(float(RAND_MAX)) - 1.0;
      v2 = 2.0*float(rand())/(float(RAND_MAX)) - 1.0;
      r = v1*v1 + v2*v2;
    }
    fac = sqrt(-2.0*log(r)/r);
    ph_gasdev_gset = v1*fac;
    ph_gasdev_ret  = v2*fac;
    ph_gasdev_iset = 1;
  } else {
    ph_gasdev_ret  = ph_gasdev_gset;
    ph_gasdev_iset = 0;
  }
  return ph_gasdev_ret;
}

/************************************************************************/

double ph_BE(double r)
{
  double factor;
  factor = 4.8089e-04*pow(r,5) - 1.0173e-02*pow(r,4) + 7.7899e-02*pow(r,3) - 
    2.3299e-01*pow(r,2) + 1.4721e-02*r + 1.0008e+00;
  return factor;
}

/************************************************************************/

double ph_q(double r)
{
  double factor;
  factor = 0.0015970*pow(r,5) - 0.0229113*pow(r,4) + 0.0386709*pow(r,3) + 
    0.7350457*pow(r,2) - 0.5490283*r + 0.0872061;
  return factor;
}

/************************************************************************/

double ph_Ang(double a1, double a2, double R, double r)
{
  return ((a2-a1)/R)*r + a1;
}

/************************************************************************/

double ph_Maxwellian(double c_tilda, double vel_unit, double mu, double gamma)
{
   // Set constants
   double mh = 1.67e-24;
   double kboltz = 1.38e-16;
   double pi = 3.14159;

   // Compute temperature in cgs units
   double c = c_tilda*vel_unit;
   double T = (pow(c,2)*mu*mh)/(kboltz*gamma);

   // Compute random Maxwellian velocity
   double mean = 0;
   double stdev = sqrt((kboltz*T)/(mu*mh));
   double u1 = rand();
   u1 = u1/RAND_MAX;
   double u2 = rand();
   u2 = u2/RAND_MAX;
   double x1 = mean + stdev*sqrt(-2*log(u1))*cos(2*pi*u2);
   double x2 = mean + stdev*sqrt(-2*log(u1))*sin(2*pi*u2);
   return (x1/vel_unit);
}
