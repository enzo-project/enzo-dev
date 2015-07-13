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
/  Free-streaming Radiation Implicit Problem Class
/  Problem initialization routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified:   
/
/  PURPOSE: Allocates all necessary internal memory for problem 
/           definition and associated linear solver.  This begins the
/           interface between Enzo and the FS solver module, so any
/           and all grid/index transformations must be performed and 
/           stored here.
/
/           This routine should only be called once (at the beginning 
/           of the simulation), so it is where necessary memory should 
/           be allocated (as opposed to the constructor, that is 
/           called prior to knowledge of the local grid sizes).
/
************************************************************************/
#ifdef TRANSFER
#include "FSProb.h"
#include "CosmologyParameters.h"

// character strings
EXTERN char outfilename[];



// Problem initializer prototypes
int FSMultiSourceInitialize(FILE *fptr, FILE *Outfptr,
			    HierarchyEntry &TopGrid,
			    TopGridData &MetaData, int local);



// Short routine to reset emissivity source magnitude
int FSProb::SetEmissivity(float NGammaDot_new)
{
  NGammaDot = NGammaDot_new;
  return SUCCESS;
}


// Main initializer routine
int FSProb::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData)
{

  //  if (debug)  printf("Entering FSProb::Initialize routine\n");


  // find the grid corresponding to this process from the Hierarcy
  HierarchyEntry *ThisGrid = &TopGrid;
  int i, dim, face, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("FSProb Initialize ERROR: p%"ISYM" could not locate his grid\n",
	   MyProcessorNumber);
    return FAIL;
  }

  // set rank of problem 
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


  // set default module parameters
  theta   = 1.0;        // backwards euler implicit time discret.
  LimType = 4;          // Zeus limiter
  EScale  = 1.0;        // no radiation equation scaling
  kappa0  = 1.0e-31;    // negligible opacity
  kappa_h2on = 0;       // no spatially dependent (use background) opacity
  for (dim=0; dim<rank; dim++)       // set default radiation boundaries to 
    for (face=0; face<2; face++)     // periodic in each direction
      BdryType[dim][face] = 0;

  // set default solver parameters
  initial_guess      = 0;         // previous time step
  sol_tolerance      = 1e-5;      // solver tolerance
  sol_printl         = 0;         // HYPRE print level
  sol_log            = 0;         // HYPRE logging level
  sol_maxit          = 50;        // HYPRE max multigrid iters
  sol_rlxtype        = 2;         // HYPRE relaxation type
  sol_npre           = 5;         // HYPRE num pre-smoothing steps
  sol_npost          = 5;         // HYPRE num post-smoothing steps

  // set default ionization parameters
  NGammaDot          = 0.0;       // ionization strength
  EtaRadius          = 0.0;       // single cell
  EtaCenter[0]       = 0.0;       // x-location
  EtaCenter[1]       = 0.0;       // y-location
  EtaCenter[2]       = 0.0;       // z-location

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
	ret += sscanf(line, "FSRadiationScaling = %"FSYM, &EScale);
	ret += sscanf(line, "FSRadiationTheta = %"FSYM, &theta);
	ret += sscanf(line, "FSRadiationOpacity = %"FSYM, &kappa0);
	ret += sscanf(line, "FSRadiationH2OpacityOn = %"ISYM, &kappa_h2on);
	ret += sscanf(line, "FSRadiationNGammaDot = %lf", &NGammaDot);
	ret += sscanf(line, "FSRadiationEtaRadius = %"FSYM, &EtaRadius);
	ret += sscanf(line, "FSRadiationEtaCenter = %"FSYM" %"FSYM" %"FSYM, 
		      &(EtaCenter[0]), &(EtaCenter[1]), &(EtaCenter[2]));
	ret += sscanf(line, "FSRadiationLimiterType = %"ISYM, &LimType);
	ret += sscanf(line, "FSRadiationBoundaryX0Faces = %"ISYM" %"ISYM, 
		      BdryType[0], BdryType[0]+1);
	if (rank > 1) {
	  ret += sscanf(line, "FSRadiationBoundaryX1Faces = %"ISYM" %"ISYM,
			BdryType[1], BdryType[1]+1);
	  if (rank > 2) {
	    ret += sscanf(line, "FSRadiationBoundaryX2Faces = %"ISYM" %"ISYM,
			  BdryType[2], BdryType[2]+1);
	  }
	}
	ret += sscanf(line, "FSRadiationMaxDt = %"FSYM, &maxdt);
	ret += sscanf(line, "FSRadiationInitialGuess = %"ISYM, &initial_guess);
	ret += sscanf(line, "FSRadiationTolerance = %g", &sol_tolerance);
	ret += sscanf(line, "FSRadiationMaxMGIters = %i", &sol_maxit);
	ret += sscanf(line, "FSRadiationMGRelaxType = %i", &sol_rlxtype);
	ret += sscanf(line, "FSRadiationMGPreRelax = %i", &sol_npre);
	ret += sscanf(line, "FSRadiationMGPostRelax = %i", &sol_npost);
	
      }  // end loop over file lines

    }  // end successful file open
  }  // end if file name exists
 
  // clean up
  delete[] dummy;
  rewind(fptr);
  fclose(fptr);

  ////////////////////////////////

  // set maximum time step into TopGrid (if used)
  if (maxdt > 0.0)
    ThisGrid->GridData->SetMaxRadiationDt(maxdt);
  
  // check that these give appropriate values, otherwise set dim to periodic
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++)
      /// ADD NEW BOUNDARY CONDITION TYPES HERE!
      if ((BdryType[dim][face] < 0) || (BdryType[dim][face] > 2)) {
	fprintf(stderr,"FSProb_Initialize Warning: re-setting BC to periodic, dim %"ISYM", face %"ISYM"\n",dim,face);
	BdryType[dim][face] = 0;
      }

  // check that periodic faces match
  for (dim=0; dim<rank; dim++) 
    if ((BdryType[dim][0]*BdryType[dim][1] == 0) && 
	(BdryType[dim][0]+BdryType[dim][1] != 0)) {
      fprintf(stderr,"FSProb_Initialize Warning: non-matching periodic BCs, dim %"ISYM"\n",dim);
      BdryType[dim][0] = 0;
      BdryType[dim][1] = 0;
    }

  // ensure that new EBdryVals array pointers are set to NULL
  for (dim=0; dim<rank; dim++) 
    for (face=0; face<2; face++) 
      BdryVals[dim][face] = NULL;


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

  // dt gives the time step size (initialize to zero)
  dt = 0.0;

  // dt_suggest gives the desired time step size (initialize to large)
  dt_suggest = huge_number;

  // a, adot give cosmological expansion & rate
  a = 1.0;
  adot = 0.0;

  // EScale gives variable scalings for implicit solver
  if (EScale <= 0.0) {
    fprintf(stderr,"Initialize: illegal FSRadiationScaling = %g\n",EScale);
    fprintf(stderr,"   re-setting to 1.0\n");
    EScale = 1.0;  // default is no scaling
  }
  if (debug)
    printf("FSProb::Initialize: EScale = %g\n", EScale);

  // Theta gives the implicit time-stepping method (1->BE, 0.5->CN, 0->FE)
  if ((theta < 0.0) || (theta > 1.0)) {
    fprintf(stderr,"FSProb Initialize: illegal FSRadiationTheta = %g\n",
	    theta);
    fprintf(stderr,"   re-setting theta to 1.0 (Backwards Euler)\n");
    theta = 1.0;  // default is backwards Euler
  }

  // kappa gives the background opacity
  if (kappa0 < 0.0) {
    fprintf(stderr,"FSProb Initialize: illegal FSRadiationOpacity = %g < 0\n",
	    kappa0);
    fprintf(stderr,"   re-setting to 1e-31\n");
    kappa0 = 1.0e-31;
  }

  // unless this is used for LW radiation, disable spatially-dependent opacity
  if ((RadiativeTransferFLD > 1) && kappa_h2on) {
    fprintf(stderr,"FSProb_Initialize Warning: kappa_h2on disabled for FLD-only problem\n");
    kappa_h2on = 0;
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

  //   for non-periodic domain, unset neighbor info.
#ifndef MPI_INT
  int MPI_PROC_NULL = -3;
  int MPI_COMM_WORLD = 0;
#endif
  for (dim=0; dim<rank; dim++) {
    if ((OnBdry[dim][0]) && (BdryType[dim][0] != 0))
      NBors[dim][0] = MPI_PROC_NULL;
    if ((OnBdry[dim][1]) && (BdryType[dim][1] != 0))
      NBors[dim][1] = MPI_PROC_NULL;
  }
  
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
		      1, NBors[0][0], NBors[0][1], NBors[1][0], 
		      NBors[1][1], NBors[2][0], NBors[2][1], empty);
  GhDims[0][0] = xghosts;
  GhDims[0][1] = xghosts;
  GhDims[1][0] = yghosts;
  GhDims[1][1] = yghosts;
  GhDims[2][0] = zghosts;
  GhDims[2][1] = zghosts;

  // set up vectors for temporary storage
  sol = U0->clone();        // linear system solution
  extsrc = U0->clone();     // emissivity sources
  if (kappa_h2on == 1)      // opacity (don't allocate if unused)
    kappa = U0->clone();
  else
    kappa = new EnzoVector(LocDims[0], LocDims[1], LocDims[2], xghosts, 
			   xghosts, yghosts, yghosts, zghosts, zghosts, 
			   1, NBors[0][0], NBors[0][1], NBors[1][0], 
			   NBors[1][1], NBors[2][0], NBors[2][1], empty);

  // initialize HYPRE stuff
#ifdef USE_HYPRE
  //    initialize the diagnostic information
  totIters = 0;

  //    set up the grid
  //       create the grid object
  HYPRE_StructGridCreate(MPI_COMM_WORLD, rank, &grid);

  //       set my grid extents as if we have one part with multiple boxes.
  //       Have each processor describe it's own global extents
  Eint32 ilower[3] = {SolvIndices[0][0], SolvIndices[1][0], SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1], SolvIndices[1][1], SolvIndices[2][1]};
  HYPRE_StructGridSetExtents(grid, ilower, iupper);

  //       set grid periodicity
  Eint32 periodicity[3] = {0, 0, 0};
  if (BdryType[0][0] == 0)  periodicity[0] = GlobDims[0];
  if (BdryType[1][0] == 0)  periodicity[1] = GlobDims[1];
  if (BdryType[2][0] == 0)  periodicity[2] = GlobDims[2];
  HYPRE_StructGridSetPeriodic(grid, periodicity);
  
  //       assemble the grid
  HYPRE_StructGridAssemble(grid);

  //   set up the stencil
  if (rank == 1) 
    stSize = 3;
  else if (rank == 2)
    stSize = 5;
  else 
    stSize = 7;
  HYPRE_StructStencilCreate(rank, stSize, &stencil);

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

  //   allocate temporary arrays
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  matentries = new Eflt64[stSize*Nx*Ny*Nz];
  rhsentries = new Eflt64[Nx*Ny*Nz];
  HYPREbuff = new Eflt64[Nx];
  HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &J);
  HYPRE_StructMatrixInitialize(J);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhsvec);
  HYPRE_StructVectorInitialize(rhsvec);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solvec);
  HYPRE_StructVectorInitialize(solvec);

#else  // ifdef USE_HYPRE

  ENZO_FAIL("FSProb usage requires HYPRE to be enabled!");

#endif


  //   check MG solver parameters
  if (sol_maxit < 0) {
    fprintf(stderr,"Illegal FSRadiationMaxMGIters = %i. Setting to 20\n",
	    sol_maxit);
    sol_maxit = 20;
  }
  if ((sol_rlxtype<0) || (sol_rlxtype>3)) {
    fprintf(stderr,"Illegal FSRadiationMGRelaxType = %i. Setting to 1\n",
	    sol_rlxtype);
    sol_rlxtype = 1;
  }
  if (sol_npre < 0) {
    fprintf(stderr,"Illegal FSRadiationMGPreRelax = %i. Setting to 1\n",
	    sol_npre);
    sol_npre = 1;
  }
  if (sol_npost < 0) {
    fprintf(stderr,"Illegal FSRadiationMGPostRelax = %i. Setting to 1\n",
	    sol_npost);
    sol_npost = 1;
  }
  if (sol_tolerance < 1.0e-15) {
    fprintf(stderr,"Illegal FSRadiationTolerance = %g. Setting to 1e-4\n",
	    sol_tolerance);
    sol_tolerance = 1.0e-4;
  }




  ////////////////////////////////
  // set up any problem-specific local data initializers here, 
  // depending on the ProblemType
  float ZERO = 0.0;
  float ONE = 1.0;
  fptr = NULL;
  switch (ProblemType) {
    
  // FSMultiSource Test
  case 450:

    // first call local problem initializer (to allocate/setup local data)
    if (FSMultiSourceInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("FSProb Initialize: Error in FSMultiSourceInitialize");

    // set BCs based on input, 0 implies periodic, otherwise set to zero-valued
    if (BdryType[0][0] != 0) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x0 left)");
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x0 right)");
    }	
    if (BdryType[1][0] != 0) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x1 left)");
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x1 right)");
    }	
    if (BdryType[2][0] != 0) {
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x2 left)");
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x2 right)");
    }	

    break;

  // FSRadWave Test
  case 451:

    ONE = 1.0e-15;
    // first call local problem initializer (to allocate/setup local data)
    // [call this other problem init because it sets a homogeneous field]
    if (FSMultiSourceInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("FSProb Initialize: Error in FSMultiSourceInitialize");

    // set x-left BC to Dirichlet, x-right to Neumann, leave others as periodic
    if (this->SetupBoundary(0,0,1,&ONE) == FAIL)
      ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x0 left)");
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL)
      ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x0 right)");

    break;

  // FSRadPoint Test
  case 452:

    // first call local problem initializer (to allocate/setup local data)
    // [call this other problem init because it sets a homogeneous field]
    if (FSMultiSourceInitialize(fptr, fptr, TopGrid, MetaData, 1) == FAIL) 
      ENZO_FAIL("FSProb Initialize: Error in FSMultiSourceInitialize");

    // set all boundaries to homogeneous Neumann
    if (this->SetupBoundary(0,0,1,&ZERO) == FAIL)
      ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x0 left)");
    if (this->SetupBoundary(0,1,1,&ZERO) == FAIL)
      ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x0 right)");
    if (this->SetupBoundary(1,0,1,&ZERO) == FAIL)
      ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x1 left)");
    if (this->SetupBoundary(1,1,1,&ZERO) == FAIL)
      ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x1 right)");
    if (this->SetupBoundary(2,0,1,&ZERO) == FAIL)
      ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x2 left)");
    if (this->SetupBoundary(2,1,1,&ZERO) == FAIL)
      ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x2 right)");

    break;

  // running a multi-frequency problem, don't initialize data but do set up BCs
  case 460:
  case 462:

    // set BCs based on input, 0 implies periodic, otherwise set to zero-valued
    if (BdryType[0][0] != 0) {
      if (this->SetupBoundary(0,0,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x0 left)");
      if (this->SetupBoundary(0,1,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x0 right)");
    }	
    if (BdryType[1][0] != 0) {
      if (this->SetupBoundary(1,0,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x1 left)");
      if (this->SetupBoundary(1,1,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x1 right)");
    }	
    if (BdryType[2][0] != 0) {
      if (this->SetupBoundary(2,0,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x2 left)");
      if (this->SetupBoundary(2,1,1,&ZERO) == FAIL)
	ENZO_FAIL("FSProb Initialize: Error in SetupBoundary (x2 right)");
    }	

    break;

  default:

    // by default do not call any local problem initializer -- assumes
    // the problem has been set up elsewhere, and that BCs are periodic
    break;

  }
  ////////////////////////////////


  // output FS problem solver parameters to output log file 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    FILE *outfptr;
    if ((outfptr = fopen(outfilename, "a")) == NULL) {
      fprintf(stderr,"Error opening parameter output file %s!!\n", 
	      outfilename);
      return FAIL;
    }
    else {
      fprintf(outfptr, "FSRadiationScaling = %g\n", EScale);
      fprintf(outfptr, "FSRadiationTheta = %g\n", theta);
      fprintf(outfptr, "FSRadiationOpacity = %g\n", kappa0);
      fprintf(outfptr, "FSRadiationH2OpacityOn = %"ISYM"\n", kappa_h2on);
      fprintf(outfptr, "FSRadiationNGammaDot = %g\n", NGammaDot);
      fprintf(outfptr, "FSRadiationEtaRadius = %g\n", EtaRadius);
      fprintf(outfptr, "FSRadiationEtaCenter = %g %g %g\n", 
	      EtaCenter[0], EtaCenter[1], EtaCenter[2]);
      fprintf(outfptr, "FSRadiationLimiterType = %"ISYM"\n", LimType);
      fprintf(outfptr, "FSRadiationBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	      BdryType[0][0], BdryType[0][1]);
      if (rank > 1) {
	fprintf(outfptr, "FSRadiationBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
		BdryType[1][0], BdryType[1][1]);
	if (rank > 2) {
	  fprintf(outfptr, "FSRadiationBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
		  BdryType[2][0], BdryType[2][1]);
	}
      }
      fprintf(outfptr, "FSRadiationMaxDt = %g\n", maxdt);
      fprintf(outfptr, "FSRadiationInitialGuess = %"ISYM"\n", initial_guess);
      fprintf(outfptr, "FSRadiationTolerance = %g\n", sol_tolerance);    
      fprintf(outfptr, "FSRadiationMaxMGIters = %i\n", sol_maxit);    
      fprintf(outfptr, "FSRadiationMGRelaxType = %i\n", sol_rlxtype);    
      fprintf(outfptr, "FSRadiationMGPreRelax = %i\n", sol_npre);    
      fprintf(outfptr, "FSRadiationMGPostRelax = %i\n", sol_npost);    
      
      // close parameter file
      fclose(outfptr);
    }
  }

  //  if (debug)  printf("Leaving FSProb::Initialize routine\n");

  return SUCCESS;
}
#endif
