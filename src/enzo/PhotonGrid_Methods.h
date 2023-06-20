/*-*-C++-*-*/
/*******************************************************************/
/*             Methods for handling Photon Packages                */
/*******************************************************************/

void SetOriginalProcessorNumber(int num) { OriginalProcessorNumber = num; };
int ReturnOriginalProcessorNumber() { return OriginalProcessorNumber; };
void DeleteSubgridMarker() { delete [] SubgridMarker; SubgridMarker = NULL; };

/* Identify radiation pressure fields */

  int IdentifyRadiationPressureFields(int &RPresNum1, int &RPresNum2,
				      int &RPresNum3);

/* Photons: Move PhotonPackages from FromGrids[] to this one */

   int MoveAllPhotonPackages(int NumberOfGrids, grid* FromGrid[]);

/* Photons: Computes photon timestep */

   float ComputePhotonTimestep(void);
   float ComputePhotonTimestepHII(float DensityUnits, float LengthUnits,
				  float VelocityUnits, float aye, 
				  float Ifront_kph);
   float ComputePhotonTimestepTau(float DensityUnits, float LengthUnits,
				  float VelocityUnits, float aye);

/* Photons: return number of PhotonPackages. */

   int ReturnNumberOfPhotonPackages(void) {return NumberOfPhotonPackages;};

/* Photons: return PhotonPackage pointer. */

   PhotonPackageEntry *ReturnPhotonPackagePointer(void) 
   {return PhotonPackages;};

   PhotonPackageEntry *ReturnPausedPackagePointer(void) 
   {return PausedPhotonPackages;};

/* Photons: set number of photons. */

   void SetNumberOfPhotonPackages(int num) {NumberOfPhotonPackages = num;};

/* remove unused photon packages */ 

   int CleanUpMovedPhotonPackages(void);

/* remove all photon packages */ 

   int DeletePhotonPackages(int DeleteHeadPointer=FALSE);

/* sort photon linked lists */

   int PhotonSortLinkedLists(void);

/* Set Subgrid Marker field */

   int SetSubgridMarkerFromSubgrid(grid *Subgrid);
   int SetSubgridMarkerFromParent(grid *Parent, int level);
   int SetSubgridMarkerFromSibling(grid *Sibling, 
				   FLOAT EdgeOffset[MAX_DIMENSION]);
   int SubgridMarkerPostParallel(HierarchyEntry **Grids[], int *NumberOfGrids);
   int SubgridMarkerPostParallelGZ(grid *Parent, HierarchyEntry **Grids[],
				   int *NumberOfGrids);
   int SetSubgridMarkerIsolatedBoundaries(void);
   int CheckSubgridMarker(void);

/* Return Subgrid Marker for a position */

  int ReturnSubgridMarker(int &cindex, FLOAT x, FLOAT y, FLOAT z);

/* Initialize photoionization and heating fields  */
 
  int InitializeRadiativeTransferFields(void);
  int AllocateInterpolatedRadiation(void);

/* Function that handle the dissociation/ionization effects of photons */
int RadiativeTransferLWShielding(PhotonPackageEntry **PP, FLOAT &dP,
				 FLOAT thisDensity, FLOAT ddr, 
				 int cellindex,  float LengthUnits,
				 int kdissH2INum, int TemperatureField,
				 float geo_correction);

int RadiativeTransferLW(PhotonPackageEntry **PP, FLOAT &dP,
			int cellindex, float tau, FLOAT photonrate, 
			float geo_correction, int kdissH2INum);

int RadiativeTransferH2II(PhotonPackageEntry **PP,
			  int cellindex, float tau, FLOAT photonrate, 
			  float geo_correction, int kdissH2IINum);


int RadiativeTransferIR(PhotonPackageEntry **PP, FLOAT &dP,
			int cellindex, float tau, FLOAT photonrate, 
			FLOAT *excessrate, float geo_correction, 
			int kphHMNum, int gammaNum);

int RadiativeTransferIonization(PhotonPackageEntry **PP, FLOAT *dPi, int cellindex, 
				int species, float tau, FLOAT photonrate, 
				FLOAT *excessrate, float geo_correction,
				const int *kphNum, int gammaNum);

int RadiativeTransferXRays(PhotonPackageEntry **PP, FLOAT *dPi, int cellindex, 
			   int species, FLOAT ddr, float tau, FLOAT geo_correction,
			   FLOAT photonrate, FLOAT *excessrate, float *ion_factor2,
			   float heat_factor, const int *kphNum, int gammaNum);

int RadiativeTransferComptonHeating(PhotonPackageEntry **PP, FLOAT *dPi, int cellindex, 
				    float LengthUnits, FLOAT photonrate, 
				    int TemperatureField, FLOAT ddr, double dN, 
				    float geo_correction, int gammaNum);

/* Functions to calculate the H2II cross section */
float LookUpCrossSectionH2II(float hnu, float T);

/* Tools for setting up temperature field for H2I Shielding */

  int InitializeTemperatureFieldForH2Shield(void);
  int FinalizeTemperatureFieldForH2Shield(void);
  int GetTemperatureFieldNumberForH2Shield(void);

/* Tools for setting up temperature field for H2II Cross Section  */

  int InitializeTemperatureFieldForH2IICrossSection(void);
  int FinalizeTemperatureFieldForH2IICrossSection(void);
  int GetTemperatureFieldNumberForH2IICrossSection(void);

/* Tools for setting up temperature field for Compton heating */

  int InitializeTemperatureFieldForComptonHeating(void);
  int FinalizeTemperatureFieldForComptonHeating(void);
  int GetTemperatureFieldNumberForComptonHeating(void);

/* Flag cells to be refined by optical depth */

  int FlagCellsToBeRefinedByOpticalDepth(void);

/* Add acceleration from radiation pressure */

  int AddRadiationPressureAcceleration(void);

/* Initialize ionized sphere around a source */

  int InitializeSource(RadiationSourceEntry *RS);

/* Communicate photon packages when rebuilding hierarchy */

  int CommunicationSendPhotonPackages(grid *ToGrid, int ToProcessor,
				      int ToNumber, int FromNumber, 
				      PhotonPackageEntry **ToPP);

  int CommunicationSendSubgridMarker(grid *ToGrid, int ToProcessor);

/* Transport Photon Packages */ 

int TransportPhotonPackages(int level, int finest_level,
			    ListOfPhotonsToMove **PhotonsToMove, 
			    int GridNum, grid **Grids0, int nGrids0, 
			    grid *ParentGrid, grid *CurrentGrid);

int ElectronFractionEstimate(float dt);
int RadiationPresent(void) { return HasRadiation; }
void SetRadiation(char value) { HasRadiation = value; }

void InitializePhotonPackages(void) {
  if (PhotonPackages == NULL) {
    PhotonPackages = new PhotonPackageEntry;
    PhotonPackages->NextPackage     = NULL;
    PhotonPackages->PreviousPackage = NULL;
  }
  if (FinishedPhotonPackages == NULL) {
    FinishedPhotonPackages = new PhotonPackageEntry;
    FinishedPhotonPackages->NextPackage = NULL;
    FinishedPhotonPackages->PreviousPackage = NULL;
  }    
  if (PausedPhotonPackages == NULL) {
    PausedPhotonPackages = new PhotonPackageEntry;
    PausedPhotonPackages->NextPackage = NULL;
    PausedPhotonPackages->PreviousPackage = NULL;
  }
  return;
}

void ResetPhotonPackagePointer(void) {
  PhotonPackages->NextPackage     = NULL;
  PhotonPackages->PreviousPackage = NULL;
  PhotonPackages->Photons         = 1.;
  PhotonPackages->Type            = 0;          
  PhotonPackages->Energy          = 0.;        
  PhotonPackages->EmissionTimeInterval= 0.;      
  PhotonPackages->EmissionTime    = 0.;  
  PhotonPackages->CurrentTime     = 0.;   
  PhotonPackages->Radius          = 0.;        
  PhotonPackages->ipix            = 0;         
  PhotonPackages->level           = 0;        
  return;
}

int MoveFinishedPhotonsBack(void) {

  PhotonPackageEntry *FPP = FinishedPhotonPackages->NextPackage;

  if (FPP != NULL) {

    // Find the end of the PhotonPackages list
    PhotonPackageEntry *PP = PhotonPackages;
    while (PP->NextPackage != NULL)
      PP = PP->NextPackage;

    FPP->PreviousPackage = PP;
    PP->NextPackage = FPP;

  }

  FinishedPhotonPackages->PreviousPackage = NULL;
  FinishedPhotonPackages->NextPackage = NULL;

  return SUCCESS;
}

float ReturnTotalNumberOfRaySegments(int RaySegNum) {
  float result = 0.0;
  int i,j,k,index;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++)
	result += BaryonField[RaySegNum][index];
    }  // loop over grid
  return result;
}

/************************************************************************
   FUNCTIONS FOR DEBUGGING PHOTON COUNTS
************************************************************************/
#define NO_DEBUG
#ifdef DEBUG
int ErrorCheckSource(void) {
  PhotonPackageEntry *PP;
  for (PP = PhotonPackages->NextPackage; PP; PP = PP->NextPackage) {
    if (PP->CurrentSource != NULL) {
      if ((PP->CurrentSource->LeafID < 0 ||
	   PP->CurrentSource->LeafID > 10000) &&
	  PP->CurrentSource->LeafID != INT_UNDEFINED) {
	printf("Bad merge...\n");
	return FAIL;
      }
    }
  }
  for (PP = FinishedPhotonPackages->NextPackage; PP; PP = PP->NextPackage) {
    if (PP->CurrentSource != NULL) {
      if ((PP->CurrentSource->LeafID < 0 ||
	   PP->CurrentSource->LeafID > 10000) &&
	  PP->CurrentSource->LeafID != INT_UNDEFINED) {
	printf("Bad merge...\n");
	return FAIL;
      }
    }
  }
  return SUCCESS;
}

int ErrorCheckPhotonNumber(int level) {
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  int count = 0, fcount = 0;
  PhotonPackageEntry *PP;
  for (PP = PhotonPackages->NextPackage; PP; PP = PP->NextPackage)
    count++;
  for (PP = FinishedPhotonPackages->NextPackage; PP; PP = PP->NextPackage)
    fcount++;
  if (count+fcount != NumberOfPhotonPackages) {
    printf("level %"ISYM", grid %"ISYM" (%x)\n", level, this->ID, this);
    printf("-> Mismatch between photon count (%"ISYM", %"ISYM") and "
	   "NumberOfPhotonPackages (%"ISYM")\n", 
	   count, fcount, NumberOfPhotonPackages);
    return FAIL;
  }
  return SUCCESS;
}

int ReturnFinishedPhotonCount(void) {
  int result = 0;
  if (MyProcessorNumber != ProcessorNumber)
    return result;
  PhotonPackageEntry *PP = FinishedPhotonPackages->NextPackage;
  while (PP != NULL) {
    result++;
    PP = PP->NextPackage;
  }
  return result;
}

int ReturnRealPhotonCount(void) {
  int result = 0;
  if (MyProcessorNumber != ProcessorNumber)
    return result;
  PhotonPackageEntry *PP = PhotonPackages->NextPackage;
  while (PP != NULL) {
    result++;
    PP = PP->NextPackage;
  }
  PP = FinishedPhotonPackages->NextPackage;
  while (PP != NULL) {
    result++;
    PP = PP->NextPackage;
  }
  return result;
}
#endif /* DEBUG */
/************************************************************************
   END -- DEBUG FUNCTIONS
************************************************************************/

int CountPhotonNumber(void) {

  if (MyProcessorNumber != ProcessorNumber)
    return 0;

  int nphotons = 0;
  PhotonPackageEntry *PP = PhotonPackages->NextPackage;
  while (PP != NULL) {
    nphotons++;
    PP = PP->NextPackage;
  }

  return nphotons;

}

int CountRadiationCells(void) {

  if (MyProcessorNumber != ProcessorNumber)
    return 0;
  
  int i, dim, ncells, size = 1;

  ncells = 0;

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				      kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    return 0;
  }
  
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  for (i = 0; i < size; i++)
    if (BaryonField[kphHINum][i] > 1) ncells++;

  return ncells;
}

float Max_kph(int &ncells) {

  if (MyProcessorNumber != ProcessorNumber)
    return 0;

  int i, dim, size = 1;
  float max_kph = 0;

  ncells = 0;

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				      kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    return 0;
  }
  
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  for (i = 0; i < size; i++) {
    if (BaryonField[kphHINum][i] > 0) ncells++;
    max_kph = max(max_kph, BaryonField[kphHINum][i]);
  }

  return max_kph;

};

float Min_kph(int &ncells) {

  if (MyProcessorNumber != ProcessorNumber)
    return 0;

  int i, dim, size = 1;
  float min_kph = 1e20;

  ncells = 0;

  int kphHINum, gammaNum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, 
				      kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    return 0;
  }
  
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  for (i = 0; i < size; i++) {
    if (BaryonField[kphHINum][i] > 0) {
      ncells++;
      min_kph = min(min_kph, BaryonField[kphHINum][i]);
    }
  }

  return min_kph;

};

float ReturnMaximumkphIfront(void) { 
  if (MyProcessorNumber == ProcessorNumber)
    return MaximumkphIfront;
  else
    return 0;
};

/* Compares min/max radiation field to estimate if there is an
   ionization front present in this grid. */

int DetectIonizationFrontApprox(float TemperatureUnits);

/* Merge any paused photons that have been stopped at the merging
   radius */

int MergePausedPhotonPackages(void);


/* Regrid a paused photon into its new spherical (HEALPix) grid
   point in preparation for merging. */

int RegridPausedPhotonPackage(PhotonPackageEntry** PP, grid* ParentGrid,
			      grid** MoveToGrid, int &DeltaLevel,
			      int &DeleteMe, const float *DomainWidth,
			      const float LightSpeed);

/* Trace a line thorugh the grid */

int TraceRay(int NumberOfSegments,
	     FLOAT r,
	     FLOAT x, FLOAT y, FLOAT z,
	     FLOAT ux, FLOAT uy, FLOAT uz,
	     FLOAT dr[],
	     long cindex[], int ci[], int cj[], int ck[]);

/* Walk Photon Package one by one */ 

int WalkPhotonPackage(PhotonPackageEntry **PP, 
		      grid **MoveToGrid, grid *ParentGrid, grid *CurrentGrid,
		      grid **Grids0, int nGrids0, int &DeleteMe, int &PauseMe, 
		      int &DeltaLevel, float LightCrossingTime,float LightSpeed,
		      int level, float MinimumPhotonFlux);

int FindPhotonNewGrid(int cindex, FLOAT *r, double *u, int *g,
		      PhotonPackageEntry* &PP,
		      grid* &MoveToGrid, int &DeltaLevel,
		      const float *DomainWidth, int &DeleteMe,
		      grid *ParentGrid);

int PhotonPeriodicBoundary(int &cindex, FLOAT *r, int *g, FLOAT *s,
			   PhotonPackageEntry* &PP, grid* &MoveToGrid, 
			   const float *DomainWidth, int &DeleteMe);

/* Create PhotonPackages for a given radiation sources   */

int Shine(RadiationSourceEntry *RadiationSource);

/* PhotonTest: Initialize grid allowing for up to ten sources  */

#define MAX_SPHERES 10
int PhotonTestInitializeGrid(int NumberOfSpheres,
			     float SphereRadius[MAX_SPHERES],
			     float SphereCoreRadius[MAX_SPHERES],
			     float SphereDensity[MAX_SPHERES],
			     float SphereTemperature[MAX_SPHERES],
			     FLOAT SpherePosition[MAX_SPHERES][MAX_DIMENSION],
			     float SphereVelocity[MAX_SPHERES][MAX_DIMENSION],
			     float SphereFracKeplarianRot[MAX_SPHERES],
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
			     char *TemperatureFilename);

int TestRadiatingStarParticleInitializeGrid(float TestStarParticleStarMass, 
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
					 FLOAT TestStarParticleStarPosition[],
           float TestStarParticleDensity,
           float TestStarParticleEnergy,
           float TestStarParticleVelocity[]);

/************************************************************************/

int StarParticlesToRadSources(FLOAT Time, double ConversionFactor);
int ConvertToCellCenteredRadiation(void);

int ReassignSuperSources(void);

int CorrectRadiationIncompleteness(void);
int FinalizeRadiationFields(void);

//***********************************************************************
// Routines for coupling to FLD solver

int DeleteEmissivity(void);
int CreateEmissivityLW(Star *AllStars, FLOAT TimeFLD, float dtFLD);
int AddRadiationImpulse(int field, double Luminosity, double sigma, 
			FLOAT BirthTime, FLOAT* pos);
