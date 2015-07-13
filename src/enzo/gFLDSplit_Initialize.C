/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Split Implicit Problem Class 
/  Initialization routine
/
/  written by: Daniel Reynolds
/  date:       July 2009
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the FLD solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"
#include "CosmologyParameters.h"
// #ifdef _OPENMP
// #include <omp.h>
// #endif

// character strings
EXTERN char outfilename[];


// function prototypes
int InitializeRateData(FLOAT Time);
int FreezeRateData(FLOAT Time, HierarchyEntry &TopGrid);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

// Function prototypes
int CosmoIonizationInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid,
			      TopGridData &MetaData, int local);
int RadHydroStreamTestInitialize(FILE *fptr, FILE *Outfptr,
				 HierarchyEntry &TopGrid,
				 TopGridData &MetaData, int local);
int RadHydroRadShockInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local);
int RadHydroPulseTestInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int RadHydroConstTestInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int RHIonizationTestInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid,
			       TopGridData &MetaData, int local);
int RHIonizationSteepInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);
int RHIonizationClumpInitialize(FILE *fptr, FILE *Outfptr,
				HierarchyEntry &TopGrid,
				TopGridData &MetaData, int local);



int gFLDSplit::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData)
{

  if (debug)  printf("Entering gFLDSplit::Initialize routine\n");


  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("FLD Initialize ERROR: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
    ENZO_FAIL("Error in gFLDSplit_Initialize");
  }

  
// #ifdef _OPENMP
//   // output number of OpenMP threads that will be used in this run
//   int nthreads = omp_get_max_threads();
//   printf("FLD Initialize: MPI task %"ISYM" has %"ISYM" available OpenMP threads\n",
//	  MyProcessorNumber,nthreads);
// #endif


#ifndef MPI_INT
  // in case MPI is not included
  int MPI_PROC_NULL = -3;
  int MPI_COMM_WORLD = 0;
#endif

  // set rank of self-gravity problem to 3
  rank = MetaData.TopGridRank;

  // get processor layout from Grid
  for (dim=0; dim<rank; dim++) 
    layout[dim] = ThisGrid->GridData->GetProcessorLayout(dim);
  
  // get processor location in MPI grid
  for (dim=0; dim<rank; dim++) 
    location[dim] = ThisGrid->GridData->GetProcessorLocation(dim);

  // get neighbor information from grid
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      NBors[dim][face] = ThisGrid->GridData->GetProcessorNeighbors(dim,face);

//   if (debug)  printf("  Initialize: setting default parameters\n");

  // set default module parameters
  Nchem  = 1;           // hydrogen only
  Model  = 1;           // case-B HII recombination coefficient
  ESpectrum = 1;        // T=10^5 blackbody spectrum
  HFrac  = 1.0;         // all Hydrogen
  theta  = 1.0;         // backwards euler implicit time discret.
  maxsubcycles = 1.0;   // step ratio between radiation and hydro
  maxchemsub = 1.0;     // step ratio between chemistry and radiation
  dtnorm = 2.0;         // use 2-norm for time step estimation
  dtgrowth = 1.1;       // 10% allowed growth in dt per step
  ErScale = 1.0;        // no radiation equation scaling
  ecScale = 1.0;        // no energy equation scaling
  NiScale = 1.0;        // no chemistry equation scaling
  int autoscale = 1;    // enable automatic variable scaling
  for (dim=0; dim<rank; dim++)     // set default radiation boundaries to 
    for (face=0; face<2; face++)   //   periodic in each direction
      BdryType[dim][face] = 0;

  // set default opacity parameters
  //   opacity computed as C0 * (rho/C1)^C2
  EnergyOpacityC0 = 1.0;
  EnergyOpacityC1 = 1.0;
  EnergyOpacityC2 = 0.0;

  // set default solver parameters
  initial_guess      = 0;         // previous time step
  sol_tolerance      = 1.0e-8;    // HYPRE solver tolerance
  sol_printl         = 1;         // HYPRE print level
  sol_log            = 1;         // HYPRE logging level
  sol_maxit          = 50;        // HYPRE max multigrid iters
  sol_rlxtype        = 1;         // HYPRE relaxation type
  sol_npre           = 1;         // HYPRE num pre-smoothing steps
  sol_npost          = 1;         // HYPRE num post-smoothing steps
  Krylov_method      = 1;         // BiCGStab outer solver

  // set default ionization parameters
  NGammaDot          = 0.0;       // ionization strength
  EtaRadius          = 0.0;       // single cell
  EtaCenter[0]       = 0.0;       // x-location
  EtaCenter[1]       = 0.0;       // y-location
  EtaCenter[2]       = 0.0;       // z-location

  // set default chemistry constants
  hnu0_HI   = 13.6;      // ionization energy of HI   [eV]
  hnu0_HeI  = 24.6;      // ionization energy of HeI  [eV]
  hnu0_HeII = 54.4;      // ionization energy of HeII [eV]

//   if (debug)  printf("  Initialize: checking input file\n");

  ////////////////////////////////
  // if input file present, over-write defaults with module inputs
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int ret;
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;

  // check whether input file is non-null
  if (MetaData.RadHydroParameterFname != NULL) {
    if ((fptr = fopen(MetaData.RadHydroParameterFname, "r")) == NULL)
      fprintf(stderr,"Error opening RadHydro parameter file %s, using defaults\n",
	      MetaData.RadHydroParameterFname);
    else {

      // read until out of lines
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	ret = 0;
	ret += sscanf(line, "RadHydroESpectrum = %"ISYM, &ESpectrum);
	ret += sscanf(line, "RadHydroChemistry = %"ISYM, &Nchem);
	ret += sscanf(line, "RadHydroHFraction = %"FSYM, &HFrac);
	ret += sscanf(line, "RadHydroModel = %"ISYM, &Model);
	ret += sscanf(line, "RadHydroMaxDt = %"FSYM, &maxdt);
	ret += sscanf(line, "RadHydroMinDt = %"FSYM, &mindt);
	ret += sscanf(line, "RadHydroInitDt = %"FSYM, &initdt);
	ret += sscanf(line, "RadHydroMaxSubcycles = %"FSYM, &maxsubcycles);
	ret += sscanf(line, "RadHydroMaxChemSubcycles = %"FSYM, &maxchemsub);
	ret += sscanf(line, "RadHydroDtNorm = %"FSYM, &dtnorm);
	ret += sscanf(line, "RadHydroDtGrowth = %"FSYM, &dtgrowth);
	ret += sscanf(line, "RadHydroDtRadFac = %"FSYM, &dtfac[0]);
	ret += sscanf(line, "RadHydroDtGasFac = %"FSYM, &dtfac[1]);
	ret += sscanf(line, "RadHydroDtChemFac = %"FSYM, &dtfac[2]);
	ret += sscanf(line, "RadiationScaling = %"FSYM, &ErScale);
	ret += sscanf(line, "EnergyCorrectionScaling = %"FSYM, &ecScale);
	ret += sscanf(line, "ChemistryScaling = %"FSYM, &NiScale);
	ret += sscanf(line, "AutomaticScaling = %"ISYM, &autoscale);
	ret += sscanf(line, "RadHydroTheta = %"FSYM, &theta);
	ret += sscanf(line, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM, 
		      BdryType[0], BdryType[0]+1);
	if (rank > 1) {
	  ret += sscanf(line, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM,
			BdryType[1], BdryType[1]+1);
	  if (rank > 2) {
	    ret += sscanf(line, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM,
			  BdryType[2], BdryType[2]+1);
	  }
	}
	ret += sscanf(line, "RadHydroInitialGuess = %"ISYM, &initial_guess);
	ret += sscanf(line, "RadHydroKrylovMethod = %"ISYM, &Krylov_method);
	ret += sscanf(line, "RadHydroSolTolerance = %"FSYM, &sol_tolerance);
	ret += sscanf(line, "RadHydroMaxMGIters = %i", &sol_maxit);
	ret += sscanf(line, "RadHydroMGRelaxType = %i", &sol_rlxtype);
	ret += sscanf(line, "RadHydroMGPreRelax = %i", &sol_npre);
	ret += sscanf(line, "RadHydroMGPostRelax = %i", &sol_npost);
	ret += sscanf(line, "EnergyOpacityC0 = %"FSYM, &EnergyOpacityC0);
	ret += sscanf(line, "EnergyOpacityC1 = %"FSYM, &EnergyOpacityC1);
	ret += sscanf(line, "EnergyOpacityC2 = %"FSYM, &EnergyOpacityC2);
	ret += sscanf(line, "NGammaDot = %"FSYM, &NGammaDot);
	ret += sscanf(line, "EtaRadius = %"FSYM, &EtaRadius);
	ret += sscanf(line, "EtaCenter = %"FSYM" %"FSYM" %"FSYM, 
		      &(EtaCenter[0]), &(EtaCenter[1]), &(EtaCenter[2]));
	
      }  // end loop over file lines
    }  // end successful file open
  }  // end if file name exists
 
  // clean up
  delete[] dummy;
  rewind(fptr);
  fclose(fptr);

//   if (debug)  printf("  Initialize: verifying inputs\n");

  ////////////////////////////////

  // check that these give appropriate values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++)
      /// ADD NEW BOUNDARY CONDITION TYPES HERE!
      if ((BdryType[dim][face] < 0) || (BdryType[dim][face] > 2)) {
	fprintf(stderr,"gFLDSplit_Initialize Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"gFLDSplit_Initialize Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
      BdryType[dim][0] = 0;
      BdryType[dim][1] = 0;
    }


  // ensure that new BdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      BdryVals[dim][face] = NULL;


//   if (debug)  printf("  Initialize: setting up subdomain information\n");

  // set up subdomain information
  //   EdgeVals gives the location of the left/right edge of the
  //      domain (no bdry) -- start with Enzo grid size
  for (dim=0; dim<rank; dim++) {
    EdgeVals[dim][0] = ThisGrid->GridData->GetGridLeftEdge(dim);
    EdgeVals[dim][1] = ThisGrid->GridData->GetGridRightEdge(dim);
  }

  //   LocDims holds the dimensions of the local domain, 
  //   active cells only (no ghost or boundary cells)
  for (dim=0; dim<rank; dim++)
    LocDims[dim] = ThisGrid->GridData->GetGridEndIndex(dim)
      - ThisGrid->GridData->GetGridStartIndex(dim) + 1;

  // Model gives the physical set of equations to use
  if ((Model != 1) && (Model != 10) && (Model != 4)) {
    fprintf(stderr,"gFLDSplit Initialize: illegal Model = %"ISYM"\n",Model);
    ENZO_FAIL("   Model is unimplemented in this module, exiting.");
  }

  // Nchem gives the number of chemical species
  if ((Nchem < 0) || (Nchem > 10)) {
    fprintf(stderr,"gFLDSplit Initialize: illegal Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"   re-setting Nchem to 0\n");
    Nchem = 1;  // default is Hydrogen only
  }

  // RadHydroHFraction must be between 0 and 1
  if ((HFrac < 0.0) || (HFrac > 1.0)) {
    fprintf(stderr,"gFLDSplit Initialize: illegal RadHydroHFraction = %g\n",
	    HFrac);
    fprintf(stderr,"   re-setting to 1.0\n");
    HFrac = 1.0;  // default is all Hydrogen
  }

  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal MaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal MinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal InitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // maxsubcycles gives the maximum desired ratio between hydro time step 
  // size and radiation time step size (dt_rad <= dt_hydro)
  if (maxsubcycles < 1.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal RadHydroMaxSubcycles = %g\n",maxsubcycles);
    fprintf(stderr,"   re-setting to 1.0\n");
    maxsubcycles = 1.0;    // default is to synchronize steps
  }

  // maxchemsub gives the maximum desired ratio between radiation time step 
  // size and chemistry time step size (dt_chem <= dt_rad)
  if (maxchemsub < 1.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal RadHydroMaxChemSubcycles = %g\n",maxchemsub);
    fprintf(stderr,"   re-setting to %g\n",1.0);
    maxchemsub = 1.0;    // default is to synchronize steps
  }

  // if using Enzo's chemistry module, warn if subcycling radiation, and reset chemsub
  if (RadiativeCooling) {
    if (maxsubcycles > 1.0) {
      fprintf(stderr,"\n**************************************************************\n");
      fprintf(stderr," WARNING: radiation subcycling (RadHydroMaxSubcycles = %g > 1.0)\n",
	      maxsubcycles);
      fprintf(stderr,"          may not work properly with Enzo chemistry module!\n");
      fprintf(stderr,"**************************************************************\n\n");
    }
    maxchemsub = 1.0;
  }

  // a, adot give cosmological expansion & rate
  a = 1.0;
  a0 = 1.0;
  adot = 0.0;
  adot0 = 0.0;

  // *Scale give variable scalings for implicit solver
  if (ErScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal RadiationScaling = %g\n",ErScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    ErScale = 1.0;  // default is no scaling
  }
  if (ecScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal EnergyCorrectionScaling = %g\n",ecScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    ecScale = 1.0;  // default is no scaling
  }
  if (NiScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal ChemistryScaling = %g\n",NiScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    NiScale = 1.0;  // default is no scaling
  }
  autoScale = (autoscale != 0);  // set bool based on integer input
  if (debug)
    printf("gFLDSplit::Initialize p%"ISYM": ErScale = %g, ecScale = %g, NiScale = %g, autoScale = %"ISYM"\n",
	   MyProcessorNumber,ErScale,ecScale,NiScale,autoscale);

  // dtfac gives the desired percent change in values per step
  if (dtfac[0] <= 0.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal DtRadFac = %g\n",dtfac[0]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[0] = huge_number;  // default is no limit
  }
  if (dtfac[1] <= 0.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal DtGasFac = %g\n",dtfac[1]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[1] = huge_number;  // default is no limit
  }
  if (dtfac[2] <= 0.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal DtChemFac = %g\n",dtfac[2]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[2] = huge_number;  // default is no limit
  }

  // dtnorm gives the norm for calculating allowed relative change per step
  if (dtnorm < 0.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal DtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 2.0 (2-norm)\n");
    dtnorm = 2.0;  // default is 2-norm
  }

  // dtgrowth gives the maximum growth factor in dt per step
  if (dtgrowth < 1.0 || dtgrowth > 10.0) {
    fprintf(stderr,"gFLDSplit Initialize: illegal RadHydroDtGrowth = %g\n",dtgrowth);
    fprintf(stderr,"   re-setting to 1.1\n");
    dtgrowth = 1.1;
  }

  // theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"gFLDSplit Initialize: illegal theta = %g\n",
	    theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
  }

  // set flags denoting if this processor is on the external boundary
  for (dim=0; dim<rank; dim++) {
    if (layout[dim]==0) {
      OnBdry[dim][0] = OnBdry[dim][1] = true;
    }
    else {
      OnBdry[dim][0] = (location[dim] == 0);
      OnBdry[dim][1] = (location[dim] == layout[dim]-1);
    }
  }
  if (debug){
    printf("gFLDSplit::Initialize p%"ISYM": rank = %"ISYM", Nchem = %"ISYM", HFrac = %g\n",
	   MyProcessorNumber, rank, Nchem, HFrac);
    printf("gFLDSplit::Initialize p%"ISYM": layout = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,layout[0],layout[1],layout[2]);
    printf("gFLDSplit::Initialize p%"ISYM": location = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,location[0],location[1],location[2]);
  }

  //   for non-periodic domain, unset neighbor info.
  for (dim=0; dim<rank; dim++) {
    if ((OnBdry[dim][0]) && (BdryType[dim][0] != 0))
      NBors[dim][0] = MPI_PROC_NULL;
    if ((OnBdry[dim][1]) && (BdryType[dim][1] != 0))
      NBors[dim][1] = MPI_PROC_NULL;
  }
  if (debug) {
//     printf("gFLDSplit::Initialize p%"ISYM": EdgeVals = (%g:%g,%g:%g,%g:%g)\n",
// 	   MyProcessorNumber, EdgeVals[0][0], EdgeVals[0][1], EdgeVals[1][0],
// 	   EdgeVals[1][1], EdgeVals[2][0], EdgeVals[2][1]);
//     printf("gFLDSplit::Initialize p%"ISYM": OnBdry = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,int(OnBdry[0][0]),int(OnBdry[0][1]),int(OnBdry[1][0]),int(OnBdry[1][1]),int(OnBdry[2][0]),int(OnBdry[2][1]));
//     printf("gFLDSplit::Initialize p%"ISYM": BdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,BdryType[0][0],BdryType[0][1],BdryType[1][0],BdryType[1][1],BdryType[2][0],BdryType[2][1]);
//     printf("gFLDSplit::Initialize p%"ISYM": NBors = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,NBors[0][0],NBors[0][1],NBors[1][0],NBors[1][1],NBors[2][0],NBors[2][1]);
  }

  // get the current units values (used to help set the time step size)
  double MassUnits;
  float TempUnits;
  DenUnits=LenUnits=TempUnits=MassUnits=TimeUnits=VelUnits=aUnits=1.0;
  if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	       &TimeUnits, &VelUnits, &MassUnits, MetaData.Time) == FAIL) 
    ENZO_FAIL("Error in GetUnits.");
  a = 1.0; adot = 0.0;
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(MetaData.Time, &a, &adot) == FAIL) 
      ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");
    aUnits = 1.0/(1.0 + InitialRedshift);
  }

  // copy initial units values into "old" units
  a0 = a;
  adot0 = adot;
  ErUnits0 = ErUnits;
  NiUnits0 = NiUnits;
  DenUnits0 = DenUnits;
  LenUnits0 = LenUnits;

  // dt* gives the time step sizes for each piece of physics
  dtrad  = initdt;                        // use the input value (scaled units)
  dtchem = max(dtrad/maxchemsub,mindt);   // use the subcycled input value

  // set a bound on the global initial dt as a factor of the radiation timestep
  dt = initdt*maxsubcycles;

  // set initial time step into TopGrid
  ThisGrid->GridData->SetMaxRadiationDt(dt);
  
  // compute global dimension information
  for (dim=0; dim<rank; dim++)
    GlobDims[dim] = MetaData.TopGridDims[dim];

  // dx gives grid cell size (comoving, normalized units)
  for (dim=0; dim<rank; dim++)
    dx[dim] = (EdgeVals[dim][1]-EdgeVals[dim][0])/LocDims[dim];

  // compute global index information for this subdomain
  float fCellsLeft;
  for (dim=0; dim<rank; dim++) {

    // the global indexing is easy if we're at the left edge
    if (location[dim]==0)  SolvIndices[dim][0]=0;

    // otherwise we compute the number of intervening cells to left edge
    else {

      // get floating point value for number of cells
      fCellsLeft = (EdgeVals[dim][0] - DomainLeftEdge[dim])/dx[dim];

      // round floating point value to closest integer
      SolvIndices[dim][0] =  (long) (fCellsLeft >= 0.0) ?
	(trunc(fCellsLeft+0.5)) : (trunc(fCellsLeft-0.5));
    }

    // add on local size to obtain right edge indices
    SolvIndices[dim][1] = SolvIndices[dim][0] + LocDims[dim] - 1;
  }

  // store local array sizes (active + ghost)
  for (dim=0; dim<rank; dim++)
    ArrDims[dim] = LocDims[dim] + 2*NumberOfGhostZones;

//   if (debug) {
//     printf("gFLDSplit::Initialize p%"ISYM": SolvIndices = (%i:%i,%i:%i,%i:%i)\n",
// 	   MyProcessorNumber, SolvIndices[0][0], SolvIndices[0][1], SolvIndices[1][0], 
// 	   SolvIndices[1][1], SolvIndices[2][0], SolvIndices[2][1]);
//     printf("gFLDSplit::Initialize p%"ISYM": SolvOff = (%"ISYM",%"ISYM",%"ISYM")\n",
// 	   MyProcessorNumber, SolvOff[0], SolvOff[1], SolvOff[2]);
//     printf("gFLDSplit::Initialize p%"ISYM": LocDims = (%"ISYM",%"ISYM",%"ISYM")\n",
// 	   MyProcessorNumber, LocDims[0], LocDims[1], LocDims[2]);
//     printf("gFLDSplit::Initialize p%"ISYM": ArrDims = (%"ISYM",%"ISYM",%"ISYM")\n",
// 	   MyProcessorNumber, ArrDims[0], ArrDims[1], ArrDims[2]);
//   }

//   if (debug)  printf("  Initialize: setting up EnzoVectors\n");

  // set up vector container for previous time step (empty data)
  int xghosts = NumberOfGhostZones, yghosts=0, zghosts=0;
  if (rank > 1) {
    yghosts = NumberOfGhostZones;
    if (rank > 2) {
      zghosts = NumberOfGhostZones;
    }
  }
  int empty=1;
  U0 = new EnzoVector(LocDims[0], LocDims[1], LocDims[2], 
		      xghosts, xghosts, yghosts, yghosts, zghosts, zghosts, 
		      2+Nchem, NBors[0][0], NBors[0][1], NBors[1][0], 
		      NBors[1][1], NBors[2][0], NBors[2][1], empty);
  GhDims[0][0] = xghosts;
  GhDims[0][1] = xghosts;
  GhDims[1][0] = yghosts;
  GhDims[1][1] = yghosts;
  GhDims[2][0] = zghosts;
  GhDims[2][1] = zghosts;

//   if (debug)
//     printf("gFLDSplit::Initialize p%"ISYM": GhDims = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
// 	   MyProcessorNumber, GhDims[0][0], GhDims[0][1], GhDims[1][0], 
// 	   GhDims[1][1], GhDims[2][0], GhDims[2][1]);

  // set up vectors for temporary storage and Jacobian components
  sol  = U0->clone();
  extsrc = U0->clone();
  FluidEnergyCorrection = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  OpacityE = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  Temperature = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  Temperature0 = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];


//   if (debug)  printf("  Initialize: setting up CoolData object\n");

  // ensure that CoolData object has been set up, and reset Hydrogen fraction
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) 
      ENZO_FAIL("Error in InitializeRateData.");
  CoolData.HydrogenFractionByMass = HFrac;

  // if performing chemistry in this module, un-scale rates for use 
  // within RadHydro solver (handles its own units) 
  if (RadiativeCooling == 0) {
    float mp = 1.67262171e-24;   // Mass of a proton [g]
    float tbase1 = TimeUnits;
    float xbase1 = LenUnits/(a*aUnits);
    float dbase1 = DenUnits*a*a*a*aUnits*aUnits*aUnits;
    float kunit  = (aUnits*aUnits*aUnits*mp) / (dbase1*tbase1);
    float kunit_3bdy  = kunit * (aUnits*aUnits*aUnits*mp) / dbase1;
    float coolunit = (aUnits*aUnits*aUnits*aUnits*aUnits*xbase1*xbase1*mp*mp) 
                     / (tbase1*tbase1*tbase1*dbase1);
    for (i=0; i<CoolData.NumberOfTemperatureBins; i++) {
      RateData.k1[i]      *= kunit;
      RateData.k2[i]      *= kunit;
      RateData.k3[i]      *= kunit;
      RateData.k4[i]      *= kunit;
      RateData.k5[i]      *= kunit;
      RateData.k6[i]      *= kunit;
      RateData.k7[i]      *= kunit;
      RateData.k8[i]      *= kunit;
      RateData.k9[i]      *= kunit;
      RateData.k10[i]     *= kunit;
      RateData.k11[i]     *= kunit;
      RateData.k12[i]     *= kunit;
      RateData.k13[i]     *= kunit;
      RateData.k13dd[i]   *= kunit;
      RateData.k14[i]     *= kunit;
      RateData.k15[i]     *= kunit;
      RateData.k16[i]     *= kunit;
      RateData.k17[i]     *= kunit;
      RateData.k18[i]     *= kunit;
      RateData.k19[i]     *= kunit;
      RateData.k20[i]     *= kunit;
      RateData.k21[i]     *= kunit;
      RateData.k22[i]     *= kunit_3bdy;
      RateData.k50[i]     *= kunit;
      RateData.k51[i]     *= kunit;
      RateData.k52[i]     *= kunit;
      RateData.k53[i]     *= kunit;
      RateData.k54[i]     *= kunit;
      RateData.k55[i]     *= kunit;
      RateData.k56[i]     *= kunit;
      CoolData.ceHI[i]    *= coolunit;
      CoolData.ceHeI[i]   *= coolunit;
      CoolData.ceHeII[i]  *= coolunit;
      CoolData.ciHI[i]    *= coolunit;
      CoolData.ciHeI[i]   *= coolunit;
      CoolData.ciHeIS[i]  *= coolunit;
      CoolData.ciHeII[i]  *= coolunit;
      CoolData.reHII[i]   *= coolunit;
      CoolData.reHeII1[i] *= coolunit;
      CoolData.reHeII2[i] *= coolunit;
      CoolData.reHeIII[i] *= coolunit;
      CoolData.brem[i]    *= coolunit;
    }
    CoolData.comp   *= coolunit;
    CoolData.piHI   *= coolunit;
    CoolData.piHeI  *= coolunit;
    CoolData.piHeII *= coolunit;
  }

//   if (debug)  printf("  Initialize: computing radiation spectrum integrals\n");

  // compute Radiation Energy spectrum integrals
  if (this->ComputeRadiationIntegrals() == FAIL) 
    ENZO_FAIL("gFLDSplit::Initialize Error in computing radiation spectrum integrals");

//   if (debug)  printf("  Initialize: initializing HYPRE data structures\n");

#ifdef USE_HYPRE

#ifdef USE_MPI
  float stime = MPI_Wtime();
#else
  float stime = 0.0;
#endif
  // initialize HYPRE stuff
  //    initialize the diagnostic information
  totIters = 0;

//   if (debug)  printf("     HYPRE_StructGridCreate\n");

  //    set up the grid
  //       create the grid object
  HYPRE_StructGridCreate(MPI_COMM_WORLD, rank, &grid);

//   if (debug)  printf("     HYPRE_StructGridSetExtents\n");

  //       set my grid extents as if we have one part with multiple boxes.
  //       Have each processor describe it's own global extents
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_StructGridSetExtents(grid, ilower, iupper);

//   if (debug)  printf("     HYPRE_StructGridSetPeriodic\n");

  //       set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (BdryType[0][0] == 0)  periodicity[0] = GlobDims[0];
  if (BdryType[1][0] == 0)  periodicity[1] = GlobDims[1];
  if (BdryType[2][0] == 0)  periodicity[2] = GlobDims[2];
  HYPRE_StructGridSetPeriodic(grid, periodicity);
  
//   if (debug)  printf("     HYPRE_StructGridAssemble\n");

  //       assemble the grid
  HYPRE_StructGridAssemble(grid);

//   if (debug)  printf("     HYPRE_StructStencilCreate\n");

  //   set up the stencil
  if (rank == 1) 
    stSize = 3;
  else if (rank == 2)
    stSize = 5;
  else 
    stSize = 7;
  HYPRE_StructStencilCreate(rank, stSize, &stencil);

//   if (debug)  printf("     HYPRE_StructStencilSetElement\n");

  //      set stencil entries
  Eint32 offset[3];
  Eint32 stentry=0;
  //         dependency to x2 left
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = -1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x1 left
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = -1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x0 left
  offset[0] = -1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to self
  offset[0] = 0;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to x0 right
  offset[0] = 1;  offset[1] = 0;  offset[2] = 0;
  HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  //         dependency to x1 right
  if (rank >= 2) {
    offset[0] = 0;  offset[1] = 1;  offset[2] = 0;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }
  //         dependency to x2 right
  if (rank == 3) {
    offset[0] = 0;  offset[1] = 0;  offset[2] = 1;
    HYPRE_StructStencilSetElement(stencil, stentry++, offset);
  }

//   if (debug)  printf("     HYPRE_Struct{Matrix/Vector}*\n");

  //   allocate temporary arrays
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  matentries = new Eflt64[stSize*Nx*Ny*Nz];
  rhsentries = new Eflt64[Nx*Ny*Nz];
  HYPREbuff = new Eflt64[Nx];
  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &P);
  HYPRE_StructMatrixInitialize(P);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_StructVectorInitialize(rhsvec);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solvec);
  HYPRE_StructVectorInitialize(solvec);

#ifdef USE_MPI
  float ftime = MPI_Wtime();
#else
  float ftime = 0.0;
#endif
  HYPREtime += ftime-stime;

#else  // ifdef USE_HYPRE

  ENZO_FAIL("gFLDSplit_Initialize ERROR: module requires USE_HYPRE to be set!");
  
#endif

//   if (debug)  printf("  Initialize: checking MG solver parameters\n");

  //   check HYPRE solver parameters
  if (Krylov_method < 0 || Krylov_method > 2) {
    fprintf(stderr,"Illegal RadHydroKrylovMethod = %"ISYM". Setting to 1 (BiCGStab\n",
	    Krylov_method);
    Krylov_method = 1;
  }
  if (sol_maxit < 0) {
    fprintf(stderr,"Illegal RadHydroMaxMGIters = %i. Setting to 20\n",
	    sol_maxit);
    sol_maxit = 20;
  }
  if ((sol_rlxtype<0) || (sol_rlxtype>3)) {
    fprintf(stderr,"Illegal RadHydroMGRelaxType = %i. Setting to 1\n",
	    sol_rlxtype);
    sol_rlxtype = 1;
  }
  if (sol_npre < 1) {
    fprintf(stderr,"Illegal RadHydroMGPreRelax = %i. Setting to 1\n",
	    sol_npre);
    sol_npre = 1;
  }
  if (sol_npost < 1) {
    fprintf(stderr,"Illegal RadHydroMGPostRelax = %i. Setting to 1\n",
	    sol_npost);
    sol_npost = 1;
  }
  if ((sol_tolerance < 1.0e-15) || (sol_tolerance > 1.0)) {
    fprintf(stderr,"Illegal RadHydroSolTolerance = %g. Setting to 1e-4\n",
	    sol_tolerance);
    sol_tolerance = 1.0e-4;
  }


//   if (debug)  printf("  Initialize: calling local problem initializers\n");

  ////////////////////////////////
  // set up the boundary conditions on the radiation field, 
  // depending on the ProblemType
  float ZERO = 0.0;
  float ONE  = 1.0;
  float SMALL = 1.0e-6;
  fptr = NULL;

  // set boundary conditions based on problem type
  // (default to homogeneous Dirichlet)
  switch (ProblemType) {
    
  // ODE test problem, set BCs based on input.  
  // 0 implies periodic, otherwise set to homogeneous Dirichlet
  case 400:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroConstTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroConstTestInitialize.");
    
    if (BdryType[0][0] != 0) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    if ( MetaData.TopGridRank >= 2 ) {
      if (BdryType[1][0] != 0) {
	if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x1 left radiation BCs.");
	if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x1 right radiation BCs.");
      }
    }
    if ( MetaData.TopGridRank == 3 ) {
      if (BdryType[2][0] != 0) {
	if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
	if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
    }
    break;
    
  // Streaming test problem: set Dirichlet BC to value of 1.0, 
  // or Neumann BC to value of 0.0; leave Periodic BC alone
  case 401:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroStreamTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroStreamTestInitialize.");
    
    //   x0, left
    if (BdryType[0][0] == 1) {
      if (this->SetupBoundary(0,0,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    else if (BdryType[0][0] == 2) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    
    //   x0, right
    if (BdryType[0][1] == 1) {
      if (this->SetupBoundary(0,1,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    else if (BdryType[0][1] == 2) {
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    
    //   x1, left
    if (BdryType[1][0] == 1) {
      if (this->SetupBoundary(1,0,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    else if (BdryType[1][0] == 2) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    
    //   x1, right
    if (BdryType[1][1] == 1) {
      if (this->SetupBoundary(1,1,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    else if (BdryType[1][1] == 2) {
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    
    //   x2, left
    if ( MetaData.TopGridRank == 3 ) {
      if (BdryType[2][0] == 1) {
	if (this->SetupBoundary(2,0,1,&ONE) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
      }
      else if (BdryType[2][0] == 2) {
	if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
      }
      
      //   x2, right
      if (BdryType[2][1] == 1) {
	if (this->SetupBoundary(2,1,1,&ONE) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
      else if (BdryType[2][1] == 2) {
	if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
    }
    break;
    
  // Pulse test problem: set Dirichlet BC to value of 1.0, 
  // or Neumann BC to value of 0.0; leave Periodic BC alone
  case 402:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroPulseTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroPulseTestInitialize.");
    
    //   x0, left
    if (BdryType[0][0] == 1) {
      if (this->SetupBoundary(0,0,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    else if (BdryType[0][0] == 2) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    
    //   x0, right
    if (BdryType[0][1] == 1) {
      if (this->SetupBoundary(0,1,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    else if (BdryType[0][1] == 2) {
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    
    //   x1, left
    if (BdryType[1][0] == 1) {
      if (this->SetupBoundary(1,0,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    else if (BdryType[1][0] == 2) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    
    //   x1, right
    if (BdryType[1][1] == 1) {
      if (this->SetupBoundary(1,1,1,&ONE) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    else if (BdryType[1][1] == 2) {
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    
    //   x2, left
    if ( MetaData.TopGridRank == 3 ) {
      if (BdryType[2][0] == 1) {
	if (this->SetupBoundary(2,0,1,&ONE) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
      }
      else if (BdryType[2][0] == 2) {
	if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 left radiation BCs.");
      }
      
      //   x2, right
      if (BdryType[2][1] == 1) {
	if (this->SetupBoundary(2,1,1,&ONE) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
      else if (BdryType[2][1] == 2) {
	if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	  ENZO_FAIL("Error setting x2 right radiation BCs.");
      }
    }
    break;
    
  // Astrophysical and Lowrie & Edwards Radiating shock test problems: 
  // set Neumann value to 0.0; leave Periodic BC alone
  case 404:
  case 405:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroRadShockInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroRadShockInitialize.");
    //   x0, left
    if (BdryType[0][0] == 2) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    }
    
    //   x0, right
    if (BdryType[0][1] == 2) {
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    
    //   x1, left
    if (BdryType[1][0] == 2) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    }
    
    //   x1, right
    if (BdryType[1][1] == 2) {
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    }
    
    //   x2, left
    if (BdryType[2][0] == 2) {
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 left radiation BCs.");
    }
    
    //   x1, right
    if (BdryType[2][1] == 2) {
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 right radiation BCs.");
    }
    
    break;
    
    
  // Ionization tests 0,1,7,8: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 410:
  case 411:
  case 417:
  case 418:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationTestInitialize.");
    
    //   x0, left
    if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 left radiation BCs.");
    //   x0, right
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 right radiation BCs.");
    //   x1, left
    if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 left radiation BCs.");
    //   x1, right
    if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 right radiation BCs.");
    //   x2, left
    if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 left radiation BCs.");
    //   x2, right
    if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
    
    
  // Ionization test 2: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all non-periodic faces.
  case 412:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationClumpInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationSteepInitialize.");
    
    if (BdryType[0][0] != 0)
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    if (BdryType[0][1] != 0)
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    if (BdryType[1][0] != 0)
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    if (BdryType[1][1] != 0)
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    if (BdryType[2][0] != 0)
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 left radiation BCs.");
    if (BdryType[2][1] != 0)
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
    
    
  // Ionization test 13: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 413:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationSteepInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RHIonizationSteepInitialize.");
    
    //   x0, left
    if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 left radiation BCs.");
    //   x0, right
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 right radiation BCs.");
    //   x1, left
    if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 left radiation BCs.");
    //   x1, right
    if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 right radiation BCs.");
    //   x2, left
    if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 left radiation BCs.");
    //   x2, right
    if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
    
    
    
  // Ionization test 14: periodic boundary conditions on all faces (store no data).
  case 414:
    // first call local problem initializer (to allocate/setup local data)
    if (CosmoIonizationInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in CosmoIonizationInitialize.");
    
    break;
    


  // Ionization test 15: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 415:
    // first call local problem initializer (to allocate/setup local data)
    if (CosmoIonizationInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in CosmoIonizationInitialize.");
    
    //   x0, left
    if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 left radiation BCs.");
    //   x0, right
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x0 right radiation BCs.");
    //   x1, left
    if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 left radiation BCs.");
    //   x1, right
    if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x1 right radiation BCs.");
    //   x2, left
    if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 left radiation BCs.");
    //   x2, right
    if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
      ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
    
    
  // Temperature jump test 16: periodic BCs on all faces
  case 416:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroConstTestInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroConstTestInitialize.");
    break;
    

    
  // Insert new problem intializers here...


  default:
    // set BCs based on inputs, for non periodic set to 0-valued
    if (BdryType[0][0] != 0)
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 left radiation BCs.");
    if (BdryType[0][1] != 0)
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    if (BdryType[1][0] != 0)
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 left radiation BCs.");
    if (BdryType[1][1] != 0)
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x1 right radiation BCs.");
    if (BdryType[2][0] != 0)
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 left radiation BCs.");
    if (BdryType[2][1] != 0)
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x2 right radiation BCs.");
    break;
  }
  ////////////////////////////////

  // if using an isothermal "model", freeze rate data, now that ICs exist
  if (Model == 4) 
    if (FreezeRateData(MetaData.Time, TopGrid) == FAIL) 
      ENZO_FAIL("Error in FreezeRateData.");


  if (debug)  printf("  Initialize: outputting parameters to log file\n");

  // output RadHydro solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      fprintf(stderr,"Error opening parameter output file %s!!\n", 
	      outfilename);
      ENZO_FAIL("Error in gFLDSplit_Initialize");
    }
    else {
      // write parameters to log file and close
      this->WriteParameters(outfptr);
      fclose(outfptr);
    }
  }

  return SUCCESS;
}
#endif
