/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Implicit Problem Class problem 
/  initialization routine
/
/  written by: Daniel Reynolds
/  date:       September, 2006
/  modified1:  June 12, 2007 by John Hayes; added code to scan for MarshakParms
/              input data.
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the FLD solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"
#include "CosmologyParameters.h"

// character strings
EXTERN char outfilename[];


// function prototypes
int InitializeRateData(FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

// Problem initializer prototypes
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
int RadHydroGreyMarshakWaveInitialize(FILE *fptr, FILE *Outfptr,
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



int gFLDProblem::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData)
{

  if (debug)  printf("Entering gFLDProblem::Initialize routine\n");


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
    ENZO_FAIL("Error in gFLDProblem_Initialize");
  }

#ifndef MPI_INT
  // in case MPI is not included
  int MPI_PROC_NULL = -3;
  int MPI_COMM_WORLD = 0;
#endif

  // set rank of self-gravity problem to 3
//   rank = 3;
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
  HFrac  = 1.0;         // all Hydrogen
  ESpectrum = 1;        // T=10^5 blackbody spectrum
  theta  = 1.0;         // backwards euler implicit time discret.
  dtnorm = 2.0;         // use 2-norm for time step estimation
  LimType = 4;          // ZEUS limiter
  ErScale = 1.0;        // no radiation equation scaling
  ecScale = 1.0;        // no energy equation scaling
  NiScale = 1.0;        // no chemistry equation scaling
  for (dim=0; dim<rank; dim++)     // set default radiation boundaries to 
    for (face=0; face<2; face++)   //   periodic in each direction
      BdryType[dim][face] = 0;

  // set default opacity parameters
  //   opacities computed as C0 * (rho/C1)^C2 & (T/C3)^C4
  PlanckOpacityC0 = 1.0;
  PlanckOpacityC1 = 1.0;
  PlanckOpacityC2 = 0.0;
  PlanckOpacityC3 = 1.0;
  PlanckOpacityC4 = 0.0;
  EnergyOpacityC0 = 1.0;
  EnergyOpacityC1 = 1.0;
  EnergyOpacityC2 = 0.0;
  EnergyOpacityC3 = 1.0;
  EnergyOpacityC4 = 0.0;

  // set default solver parameters
  approx_jac         = 0;         // analytical jacobian
  initial_guess      = 0;         // previous time step
  AnalyticChem       = 1;         // use analytical QSS approach
  newt_linesearch    = 1;         // use a linesearch in the Newton solver
  newt_maxit         = 20;        // 20 Newton iterations
  newt_norm          = 0;         // standard RMS norm
  newt_INconst       = 1.0e-8;    // inexact Newton forcing constant
  newt_tol           = 1.0e-6;    // default nonlinear tolerance
  newt_MinLinesearch = 1.0e-12;   // minimum linesearch step length
  sol_relch          = 0;         // HYPRE relative change stopping crit.
  sol_printl         = 1;         // HYPRE print level
  sol_log            = 1;         // HYPRE logging level
  sol_zeroguess      = 1;         // HYPRE uses a zero initial guess
  sol_maxit          = 50;        // HYPRE max multigrid iters
  sol_rlxtype        = 1;         // HYPRE relaxation type
  sol_npre           = 1;         // HYPRE num pre-smoothing steps
  sol_npost          = 1;         // HYPRE num post-smoothing steps

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
  int numMarshakParms = 1;

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
	ret += sscanf(line, "RadHydroDtNorm = %"FSYM, &dtnorm);
	ret += sscanf(line, "RadHydroDtRadFac = %"FSYM, &dtfac[0]);
	ret += sscanf(line, "RadHydroDtGasFac = %"FSYM, &dtfac[1]);
	ret += sscanf(line, "RadHydroDtChemFac = %"FSYM, &dtfac[2]);
	ret += sscanf(line, "RadiationScaling = %"FSYM, &ErScale);
	ret += sscanf(line, "EnergyCorrectionScaling = %"FSYM, &ecScale);
	ret += sscanf(line, "ChemistryScaling = %"FSYM, &NiScale);
	ret += sscanf(line, "RadHydroTheta = %"FSYM, &theta);
	ret += sscanf(line, "RadHydroLimiterType = %"ISYM, &LimType);
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
	ret += sscanf(line, "RadHydroAprxJacobian = %"ISYM, &approx_jac);
	ret += sscanf(line, "RadHydroInitialGuess = %"ISYM, &initial_guess);
	ret += sscanf(line, "RadHydroAnalyticChem = %"ISYM, &AnalyticChem);
	ret += sscanf(line, "RadHydroNewtLinesearch = %"ISYM, &newt_linesearch);
	ret += sscanf(line, "RadHydroNewtIters = %"ISYM, &newt_maxit);
	ret += sscanf(line, "RadHydroNewtNorm = %"ISYM, &newt_norm);
	ret += sscanf(line, "RadHydroINConst = %"FSYM, &newt_INconst);
	ret += sscanf(line, "RadHydroNewtTolerance = %"FSYM, &newt_tol);
	ret += sscanf(line, "RadHydroMinLinesearch = %"FSYM,
		      &newt_MinLinesearch);
	ret += sscanf(line, "RadHydroMaxMGIters = %i", &sol_maxit);
	ret += sscanf(line, "RadHydroMGRelaxType = %i", &sol_rlxtype);
	ret += sscanf(line, "RadHydroMGPreRelax = %i", &sol_npre);
	ret += sscanf(line, "RadHydroMGPostRelax = %i", &sol_npost);
	ret += sscanf(line, "PlanckOpacityC0 = %"FSYM, &PlanckOpacityC0);
	ret += sscanf(line, "PlanckOpacityC1 = %"FSYM, &PlanckOpacityC1);
	ret += sscanf(line, "PlanckOpacityC2 = %"FSYM, &PlanckOpacityC2);
	ret += sscanf(line, "PlanckOpacityC3 = %"FSYM, &PlanckOpacityC3);
	ret += sscanf(line, "PlanckOpacityC4 = %"FSYM, &PlanckOpacityC4);
	ret += sscanf(line, "EnergyOpacityC0 = %"FSYM, &EnergyOpacityC0);
	ret += sscanf(line, "EnergyOpacityC1 = %"FSYM, &EnergyOpacityC1);
	ret += sscanf(line, "EnergyOpacityC2 = %"FSYM, &EnergyOpacityC2);
	ret += sscanf(line, "EnergyOpacityC3 = %"FSYM, &EnergyOpacityC3);
	ret += sscanf(line, "EnergyOpacityC4 = %"FSYM, &EnergyOpacityC4);
	
      }  // end loop over file lines

      // if doing a Marshak-type problem (20 <= Model < 30), input 
      // additional Marshak parameters from the input file
      if ( Model >= 20 && Model <= 29 ) {
        MarshakParms = new float[numMarshakParms];
	rewind(fptr);
        while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	  ret = 0;
	  ret += sscanf(line, "SuOlsonGreyEps = %"FSYM, &MarshakParms[0]);
        }  // end loop over file lines
//         if (debug) printf("gFLDProblem_Initialize: SuOlsonGreyEps = %g\n",MarshakParms[0]);
      }  // end Model IF statement

      // if doing an ionization problem (ProblemTypes 410-415),  
      // input additional parameters 
      rewind(fptr);
      if ((ProblemType >= 410) && (ProblemType <= 415)) {
	IonizationParms[0] = 0.0;  // set defaults
	IonizationParms[1] = 0.0;
	IonizationParms[2] = 0.0;
	IonizationParms[3] = 0.0;
	IonizationParms[4] = 0.0;
        while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	  ret = 0;
	  ret += sscanf(line, "NGammaDot = %"FSYM, &IonizationParms[0]);
	  ret += sscanf(line, "EtaRadius = %"FSYM, &IonizationParms[1]);
	  ret += sscanf(line, "EtaCenter = %"FSYM" %"FSYM" %"FSYM, 
			&IonizationParms[2], &IonizationParms[3], &IonizationParms[4]);
        }  // end loop over file lines
        if (debug) {
	  printf("gFLDProblem_Initialize: NGammaDot = %g\n",IonizationParms[0]);
	  printf("gFLDProblem_Initialize: EtaRadius = %g\n",IonizationParms[1]);
	  printf("gFLDProblem_Initialize: EtaCenter = %g %g %g\n",
		 IonizationParms[2],IonizationParms[3],IonizationParms[4]);
	}
      }  // end ProblemType IF statement

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
	fprintf(stderr,"gFLDProblem_Initialize Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"gFLDProblem_Initialize Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
      BdryType[dim][0] = 0;
      BdryType[dim][1] = 0;
    }


  // ensure that new EBdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      EBdryVals[dim][face] = NULL;


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

  // Nchem gives the number of chemical species
  if ((Nchem < 0) || (Nchem > 10)) {
    fprintf(stderr,"gFLDProblem Initialize: illegal Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"   re-setting Nchem to 0\n");
    Nchem = 0;  // default is no chemistry
  }

  // RadHydroHFraction must be between 0 and 1
  if ((HFrac < 0.0) || (HFrac > 1.0)) {
    fprintf(stderr,"gFLDProblem Initialize: illegal RadHydroHFraction = %g\n",
	    HFrac);
    fprintf(stderr,"   re-setting to 1.0\n");
    HFrac = 1.0;  // default is all Hydrogen
  }

  // LimType gives the limiter formula to use (see header)
  if ((LimType < 0) || (LimType > 4)) {
    fprintf(stderr,"gFLDProblem Initialize: illegal LimType = %"ISYM"\n",LimType);
    fprintf(stderr,"   re-setting LimType to 4 (ZEUS)\n");
    LimType = 0;  // default is ZEUS limiter
  }

  // check that AnalyticChem is enabled for this Model
  if ((AnalyticChem) && ((Model != 1) && (Model != 4) || (Nchem == 0))) {
    AnalyticChem = 0;  // disable AnalyticChem for this problem
  }

  // AnalyticChem requires an approximate Jacobian
  if ((AnalyticChem) && (approx_jac==0)) {
    approx_jac = 1;  // turn on approximate Jacobian
  }

  // maxdt gives the maximum radiation time step size
  if (maxdt <= 0.0) {
    fprintf(stderr,"gFLDProblem Initialize: illegal MaxDt = %g\n",maxdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    maxdt = huge_number;  // default is no limit
  }

  // mindt gives the minimum radiation time step size
  if (mindt < 0.0) {
    fprintf(stderr,"gFLDProblem Initialize: illegal MinDt = %g\n",mindt);
    fprintf(stderr,"   re-setting to %g\n",0.0);
    mindt = 0.0;  // default is 0.0
  }

  // initdt gives the initial time step size
  if (initdt <= 0.0) {
    fprintf(stderr,"gFLDProblem Initialize: illegal InitDt = %g\n",initdt);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    initdt = huge_number;  // default is no limit
  }

  // dt gives the time step size (initialize to zero)
  dt = 0.0;

  // a, adot give cosmological expansion & rate
  a = 1.0;
  adot = 0.0;

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
  if (debug)
    printf("gFLDProblem::Initialize p%"ISYM": ErScale = %g, ecScale = %g, NiScale = %g\n",
	   MyProcessorNumber,ErScale,ecScale,NiScale);

  // dtfac gives the desired percent change in values per step
  if (dtfac[0] <= 0.0) {
    fprintf(stderr,"gFLDProblem Initialize: illegal DtRadFac = %g\n",dtfac[0]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[0] = huge_number;  // default is no limit
  }
  if (dtfac[1] <= 0.0) {
    fprintf(stderr,"gFLDProblem Initialize: illegal DtGasFac = %g\n",dtfac[1]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[1] = huge_number;  // default is no limit
  }
  if (dtfac[2] <= 0.0) {
    fprintf(stderr,"gFLDProblem Initialize: illegal DtChemFac = %g\n",dtfac[2]);
    fprintf(stderr,"   re-setting to %g\n",huge_number);
    dtfac[2] = huge_number;  // default is no limit
  }

  // dtnorm gives the norm for calculating allowed relative change per step
  if (dtnorm < 0.0) {
    fprintf(stderr,"gFLDProblem Initialize: illegal DtNorm = %g\n",dtnorm);
    fprintf(stderr,"   re-setting to 2.0 (2-norm)\n");
    dtnorm = 2.0;  // default is 2-norm
  }

  // Theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"gFLDProblem Initialize: illegal theta = %g\n",
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
    printf("gFLDProblem::Initialize p%"ISYM": rank = %"ISYM", Nchem = %"ISYM", HFrac = %g\n",
	   MyProcessorNumber, rank, Nchem, HFrac);
    printf("gFLDProblem::Initialize p%"ISYM": layout = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,layout[0],layout[1],layout[2]);
    printf("gFLDProblem::Initialize p%"ISYM": location = (%"ISYM",%"ISYM",%"ISYM")\n",MyProcessorNumber,location[0],location[1],location[2]);
  }

  //   for non-periodic domain, unset neighbor info.
  for (dim=0; dim<rank; dim++) {
    if ((OnBdry[dim][0]) && (BdryType[dim][0] != 0))
      NBors[dim][0] = MPI_PROC_NULL;
    if ((OnBdry[dim][1]) && (BdryType[dim][1] != 0))
      NBors[dim][1] = MPI_PROC_NULL;
  }
  if (debug) {
    printf("gFLDProblem::Initialize p%"ISYM": EdgeVals = (%g:%g,%g:%g,%g:%g)\n",
	   MyProcessorNumber, EdgeVals[0][0], EdgeVals[0][1], EdgeVals[1][0],
	   EdgeVals[1][1], EdgeVals[2][0], EdgeVals[2][1]);
    printf("gFLDProblem::Initialize p%"ISYM": OnBdry = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,int(OnBdry[0][0]),int(OnBdry[0][1]),int(OnBdry[1][0]),int(OnBdry[1][1]),int(OnBdry[2][0]),int(OnBdry[2][1]));
    printf("gFLDProblem::Initialize p%"ISYM": BdryType = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,BdryType[0][0],BdryType[0][1],BdryType[1][0],BdryType[1][1],BdryType[2][0],BdryType[2][1]);
    printf("gFLDProblem::Initialize p%"ISYM": NBors = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,NBors[0][0],NBors[0][1],NBors[1][0],NBors[1][1],NBors[2][0],NBors[2][1]);
  }

  // set initial time step into TopGrid
  ThisGrid->GridData->SetMaxRadiationDt(initdt);
  
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

  // adjust SolvIndices, SolvOff for Dirichlet boundary zones
  for (dim=0; dim<rank; dim++) {
    SolvOff[dim] = 0;
    if (OnBdry[dim][0]  && (BdryType[dim][0] == 1)) {
      SolvIndices[dim][0] -= 1;
      SolvOff[dim] = 1;
    }
    if (OnBdry[dim][1]  && (BdryType[dim][1] == 1))
      SolvIndices[dim][1] += 1;
  }

  // store local array sizes (active + ghost)
  for (dim=0; dim<rank; dim++)
    ArrDims[dim] = LocDims[dim] + 2*NumberOfGhostZones;

  if (debug) {
    printf("gFLDProblem::Initialize p%"ISYM": SolvIndices = (%i:%i,%i:%i,%i:%i)\n",
	   MyProcessorNumber, SolvIndices[0][0], SolvIndices[0][1], SolvIndices[1][0], 
	   SolvIndices[1][1], SolvIndices[2][0], SolvIndices[2][1]);
    printf("gFLDProblem::Initialize p%"ISYM": SolvOff = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber, SolvOff[0], SolvOff[1], SolvOff[2]);
    printf("gFLDProblem::Initialize p%"ISYM": LocDims = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber, LocDims[0], LocDims[1], LocDims[2]);
    printf("gFLDProblem::Initialize p%"ISYM": ArrDims = (%"ISYM",%"ISYM",%"ISYM")\n",
	   MyProcessorNumber, ArrDims[0], ArrDims[1], ArrDims[2]);
  }

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

  if (debug)
    printf("gFLDProblem::Initialize p%"ISYM": GhDims = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	   MyProcessorNumber, GhDims[0][0], GhDims[0][1], GhDims[1][0], 
	   GhDims[1][1], GhDims[2][0], GhDims[2][1]);

  // set up vectors for temporary storage and Jacobian components
  sol  = U0->clone();
  tmp1 = U0->clone();
  tmp2 = U0->clone();
  tmp3 = U0->clone();
  rhs  = U0->clone();
  rhs0 = U0->clone();
  extsrc = U0->clone();
  FluidEnergyCorrection = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  Temp = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  OpacityP = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  OpacityE = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
//   Eta = new float[ArrDims[0]*ArrDims[1]*ArrDims[2]];
  L = (EnzoVector **) new EnzoVector*[2+Nchem];
  for (i=0; i<(2+Nchem); i++)  L[i] = U0->clone();

  // set up InexactNewton solver
  INSolve = new InexactNewtonSolver(sol);


//   if (debug)  printf("  Initialize: setting up CoolData object\n");

  // ensure that CoolData object has been set up
  if (CoolData.ceHI == NULL) 
    if (InitializeRateData(MetaData.Time) == FAIL) 
      ENZO_FAIL("Error in InitializeRateData.");
  // un-scale rates for use within RadHydro solver (handles its own units)
  {
    DenUnits = LenUnits = TempUnits = TimeUnits = VelUnits = MassUnits = aUnits = 1.0;
    if (GetUnits(&DenUnits, &LenUnits, &TempUnits, 
	         &TimeUnits, &VelUnits, &MassUnits, MetaData.Time) == FAIL) 
      ENZO_FAIL("Error in GetUnits.");
    float mp = 1.67262171e-24;   // Mass of a proton [g]
    a = 1.0; adot = 0.0;
    if (ComovingCoordinates) {
      if (CosmologyComputeExpansionFactor(MetaData.Time, &a, &adot) == FAIL) 
        ENZO_FAIL("Error in CosmologyComputeExpansionFactor.");
      aUnits = 1.0/(1.0 + InitialRedshift);
    }
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
    ENZO_FAIL("gFLDProblem::Initialize Error in computing radiation spectrum integrals");

//   if (debug)  printf("  Initialize: initializing HYPRE data structures\n");

#ifdef USE_HYPRE
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
  Ptmpvec = new Eflt64[stSize*Nx*Ny*Nz];
  HYPREbuff = new Eflt64[Nx];
  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &P);
  HYPRE_StructMatrixInitialize(P);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_StructVectorInitialize(rhsvec);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solvec);
  HYPRE_StructVectorInitialize(solvec);

#else  // ifdef USE_HYPRE

  ENZO_FAIL("gFLDProblem_Initialize ERROR: module requires USE_HYPRE to be set!");
  
#endif

//   if (debug)  printf("  Initialize: checking MG solver parameters\n");

  //   check MG solver parameters
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

//   if (debug)  printf("  Initialize: checking Newton solver parameters\n");

  //   check Newton solver parameters
  if (newt_maxit < 1) {
    fprintf(stderr,"Illegal RadHydroNewtIters = %"ISYM". Setting to 20\n",
	    newt_maxit);
    newt_maxit = 20;
  }
  if ((newt_norm < 0) || (newt_norm > 5)) {
    fprintf(stderr,"Illegal RadHydroNewtNorm = %"ISYM". Setting to 0\n",
	    newt_norm);
    newt_norm = 0;
  }
  if ((newt_tol < 1.0e-15) || (newt_tol > 1.0)) {
    fprintf(stderr,"Illegal RadHydroNewtTolerance = %g. Setting to 1e-4\n",
	    newt_tol);
    newt_tol = 1.0e-4;
  }
  if ((newt_INconst <= 0.0e0) || (newt_INconst >= 1.)) {
    fprintf(stderr,"Illegal RadHydroINConst = %g. Setting to 1.0e-8\n",
	    newt_INconst);
    newt_INconst = 1.0e-8;
  }
  if ((newt_MinLinesearch < 1.0e-15) || (newt_MinLinesearch > 1.0e-3)) {
    fprintf(stderr,"Illegal RadHydroMinLinesearch = %g. Setting to 1e-12\n",
	    newt_MinLinesearch);
    newt_MinLinesearch = 1.0e-12;
  }
  // if linesearch enabled (default), force INconst to be small
  if ((newt_linesearch) && (newt_INconst > 1e-4)) {
    if (debug) 
      printf("Warning: when using linesearch, RadHydroINConst should be small.  Resettingto 1.0e-4\n");
    newt_INconst = 1.0e-4;
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
    if (RadHydroConstTestInitialize(fptr, fptr, TopGrid, 
				    MetaData, 1) == FAIL) 
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
    if (RadHydroStreamTestInitialize(fptr, fptr, TopGrid, 
				     MetaData, 1) == FAIL) 
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
    if (RadHydroPulseTestInitialize(fptr, fptr, TopGrid, 
				    MetaData, 1) == FAIL) 
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
    
  // Grey Marshak test problem: set Dirichlet BC to value of 1.0, 
  // or Neumann BC to value of 0.0; leave Periodic BC alone
  case 403:
    // first call local problem initializer (to allocate/setup local data)
    if (RadHydroGreyMarshakWaveInitialize(fptr, fptr, TopGrid, 
					  MetaData, 1) == FAIL) 
      ENZO_FAIL("Error in RadHydroGreyMarshakWaveInitialize.");
    
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
      //if (this->SetupBoundary(0,1,1,&ONE) == FAIL) 
      if (this->SetupBoundary(0,1,1,&SMALL) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
    }
    else if (BdryType[0][1] == 2) {
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL) 
	ENZO_FAIL("Error setting x0 right radiation BCs.");
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
    
    
  // Ionization tests 0 and 1: set zero-gradient (homogeneous Neumann)
  // boundary conditions on all faces.
  case 410:
  case 411:
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
  // boundary conditions on all faces.
  case 412:
    // first call local problem initializer (to allocate/setup local data)
    if (RHIonizationClumpInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
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
    // set BC on all faces to homogeneous (if not periodic)
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

//   if (debug)  printf("  Initialize: outputting parameters to log file\n");

  // output RadHydro solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      fprintf(stderr,"Error opening parameter output file %s!!\n", 
	      outfilename);
      return FAIL;
    }
    else {
      fprintf(outfptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
      fprintf(outfptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
      fprintf(outfptr, "RadHydroModel = %"ISYM"\n", Model);
      fprintf(outfptr, "RadHydroMaxDt = %g\n", maxdt);
      fprintf(outfptr, "RadHydroMinDt = %g\n", mindt);
      fprintf(outfptr, "RadHydroInitDt = %g\n", initdt);
      fprintf(outfptr, "RadHydroDtNorm = %"FSYM"\n", dtnorm);
      fprintf(outfptr, "RadHydroDtRadFac = %g\n", dtfac[0]);
      fprintf(outfptr, "RadHydroDtGasFac = %g\n", dtfac[1]);
      fprintf(outfptr, "RadHydroDtChemFac = %g\n", dtfac[2]);
      fprintf(outfptr, "RadiationScaling = %g\n", ErScale);
      fprintf(outfptr, "EnergyCorrectionScaling = %g\n", ecScale);
      fprintf(outfptr, "ChemistryScaling = %g\n", NiScale);
      fprintf(outfptr, "RadHydroTheta = %g\n", theta);
      fprintf(outfptr, "RadHydroLimiterType = %"ISYM"\n", LimType);
      fprintf(outfptr, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	      BdryType[0][0], BdryType[0][1]);
      if (rank > 1) {
	fprintf(outfptr, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
		BdryType[1][0], BdryType[1][1]);
	if (rank > 2) {
	  fprintf(outfptr, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
		  BdryType[2][0], BdryType[2][1]);
	}
      }
      fprintf(outfptr, "RadHydroAprxJacobian = %"ISYM"\n", approx_jac);    
      fprintf(outfptr, "RadHydroInitialGuess = %"ISYM"\n", initial_guess);    
      fprintf(outfptr, "RadHydroAnalyticChem = %"ISYM"\n", AnalyticChem);
      fprintf(outfptr, "RadHydroNewtLinesearch = %"ISYM"\n", newt_linesearch);
      fprintf(outfptr, "RadHydroNewtIters = %"ISYM"\n", newt_maxit);    
      fprintf(outfptr, "RadHydroNewtNorm = %"ISYM"\n", newt_norm);    
      fprintf(outfptr, "RadHydroINConst = %g\n", newt_INconst);    
      fprintf(outfptr, "RadHydroNewtTolerance = %g\n", newt_tol);    
      fprintf(outfptr, "RadHydroMinLinesearch = %g\n", 
	      newt_MinLinesearch);    
      fprintf(outfptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
      fprintf(outfptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
      fprintf(outfptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
      fprintf(outfptr, "RadHydroMGPostRelax = %i\n", sol_npost);    
      fprintf(outfptr, "PlanckOpacityC0 = %g\n", PlanckOpacityC0);
      fprintf(outfptr, "PlanckOpacityC1 = %g\n", PlanckOpacityC1);
      fprintf(outfptr, "PlanckOpacityC2 = %g\n", PlanckOpacityC2);
      fprintf(outfptr, "PlanckOpacityC3 = %g\n", PlanckOpacityC3);
      fprintf(outfptr, "PlanckOpacityC4 = %g\n", PlanckOpacityC4);
      fprintf(outfptr, "EnergyOpacityC0 = %g\n", EnergyOpacityC0);
      fprintf(outfptr, "EnergyOpacityC1 = %g\n", EnergyOpacityC1);
      fprintf(outfptr, "EnergyOpacityC2 = %g\n", EnergyOpacityC2);
      fprintf(outfptr, "EnergyOpacityC3 = %g\n", EnergyOpacityC3);
      fprintf(outfptr, "EnergyOpacityC4 = %g\n", EnergyOpacityC4);
      
      // close parameter file
      fclose(outfptr);
    }
  }

  return SUCCESS;
}
#endif
