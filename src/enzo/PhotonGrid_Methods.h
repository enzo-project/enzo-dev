/*-*-C++-*-*/
/*******************************************************************/
/*             Methods for handling Photon Packages                */
/*******************************************************************/

/* Identify radiation pressure fields */

  int IdentifyRadiationPressureFields(int &RPresNum1, int &RPresNum2,
				      int &RPresNum3);

/* Photons: Move PhotonPackages from FromGrids[] to this one */

   int MoveAllPhotonPackages(int NumberOfGrids, grid* FromGrid[]);

/* Photons: Computes photon timestep */

//   float ComputePhotonTimestep(void);
   float ComputeRT_TimeStep(void);
   float ComputeRT_TimeStep2(float DensityUnits, float LengthUnits);

/* Photons: return number of PhotonPackages. */

   int ReturnNumberOfPhotonPackages(void) {return NumberOfPhotonPackages;};

/* Photons: return PhotonPackage pointer. */

   PhotonPackageEntry *ReturnPhotonPackagePointer(void) {return PhotonPackages;};

/* Photons: set number of photons. */

   void SetNumberOfPhotonPackages(int num) {NumberOfPhotonPackages = num;};

/* remove unused photon packages */ 

   int CleanUpMovedPhotonPackages(void);

/* remove all photon packages */ 

   int DeletePhotonPackages(void);

/* Set Subgrid Marker field */

   int SetSubgridMarkerFromSubgrid(grid *Subgrid, grid *CurrentGrid);

/* Return Subgrid Marker for a position */

  int ReturnSubgridMarker(int &cindex, FLOAT x, FLOAT y, FLOAT z);

/* Initialize photoionization and heating fields  */
 
  int InitializeRadiativeTransferFields(void);
  int AllocateInterpolatedRadiation(void);

/* Flag cells to be refined by optical depth */

  int FlagCellsToBeRefinedByOpticalDepth(void);

/* Add acceleration from radiation pressure */

  int AddRadiationPressureAcceleration(void);

/* Solve cooling/rate equations coupled to the radiative transfer */

  int SolveCoupledRateEquations();

/* Initialize ionized sphere around a source */

  int InitializeSource(RadiationSourceEntry *RS);

/* Communicate photon packages when rebuilding hierarchy */

  int CommunicationSendPhotonPackages(grid *ToGrid, int ToProcessor,
				      int ToNumber, int FromNumber, 
				      PhotonPackageEntry **ToPP);

/* Transport Photon Packages */ 

int TransportPhotonPackages(int level, ListOfPhotonsToMove **PhotonsToMove, 
			    int GridNum, grid **Grids0, int nGrids0, 
			    grid *ParentGrid, grid *CurrentGrid);

int ElectronFractionEstimate(float dt);
int RadiationPresent(void) { return HasRadiation; }

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

/************************************************************************
   UNUSED FUNCTIONS (FOR DEBUGGING PHOTON COUNTS)
************************************************************************/

//#ifdef UNUSED
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
//#endif /* UNUSED */
/************************************************************************
   END -- UNUSED FUNCTIONS
************************************************************************/

int CountPhotonNumber(void) {

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  NumberOfPhotonPackages = 0;
  PhotonPackageEntry *PP = PhotonPackages->NextPackage;
  while (PP != NULL) {
    NumberOfPhotonPackages++;
    PP = PP->NextPackage;
  }

  return SUCCESS;

}

int CountRadiationCells(void) {

  if (MyProcessorNumber != ProcessorNumber)
    return 0;
  
  int i, dim, ncells, size = 1;

  ncells = 0;

  int kphHINum, gammaHINum, kphHeINum, gammaHeINum, kphHeIINum, gammaHeIINum,
    kdissH2INum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
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

  int kphHINum, gammaHINum, kphHeINum, gammaHeINum, kphHeIINum, gammaHeIINum,
    kdissH2INum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    return 0;
  }
  
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  for (i = 0; i < size; i++) {
    if (BaryonField[kphHINum][i] > 0) ncells++;
    max_kph = max(max_kph, BaryonField[kdissH2INum][i]);
  }

  return max_kph;

};

/* Compares min/max radiation field to estimate if there is an
   ionization front present in this grid. */

int DetectIonizationFrontApprox(float TemperatureUnits);

/* Merge any paused photons that have been stopped at the merging
   radius */

int MergePausedPhotonPackages(void);

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
			 grid **Grids0, int nGrids0,
			 int DensNum, int HINum, int HeINum,
			 int HeIINum, int H2INum,
			 int kphHINum, int gammaHINum,
			 int kphHeINum, int gammaHeINum,
			 int kphHeIINum, int gammaHeIINum,
			 int kdissH2INum, int RPresNum1, int RPresNum2, 
			 int RPresNum3, int &DeleteMe, int &PauseMe,
			 int &DeltaLevel, float DensityUnits, 
			 float TemperatureUnits, float VelocityUnits, 
			 float LengthUnits, float TimeUnits);

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
				int RefineByOpticalDepth);

/************************************************************************/

int ConvertToCellCenteredRadiation(void);

int ReassignSuperSources(void);

int CorrectRadiationIncompleteness(void);
int FinalizeRadiationFields(void);
