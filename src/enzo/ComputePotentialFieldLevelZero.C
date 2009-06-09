/***********************************************************************
/
/  COMPUTE THE POTENTIAL FIELD
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:  Daniel R. Reynolds
/  date:       February, 2006
/  modified2:  Robert Harkness (to eliminate Ncpu^2 scaling)
/  date:       December, 2007
/  modified3:  Robert Harkness
/  date:       February, 2008
/
/  PURPOSE:
/
************************************************************************/
/***********************************************************************
/
/  COMPUTE THE GRAVITATIONAL POTENTIAL FIELD
/
/  This file extends Greg Bryan's Original code that uses an FFT-based
/  algorithm to solve for the gravitational potential under periodic 
/  boundary conditions on the root (level 0) grid.
/
/  Additional functionality has been added to allow for isolating 
/  (Dirichlet) boundary conditions on the root grid.  This solve calls
/  the MGMPI library for solution of the poisson equation.
/
/  NOTE: both approaches compute and store relevant solver information 
/  during the first call to the routine.  Neither of these routines 
/  'clean up' after themselves upon exit of the program, i.e. memory 
/  is allocated but never freed, requiring that the compiler take care 
/  of the remaining data upon program completion.  A future version of 
/  this gravity solver module may include a C++ class for the solver, 
/  which is initialized at the same point as the Enzo grids, stores 
/  all necessary information internally and privately, and is cleared 
/  prior to exit of the overall Enzo program.
/
/  written by: Greg Bryan
/  date:       January, 1998
/
/  modified1:  Daniel R. Reynolds
/  date:       August, 2005
/
************************************************************************/

 
#ifdef USE_MPI
#include <mpi.h>
#endif /* USE_MPI */
#include <string.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "communication.h" 

/* Function prototypes */

int CommunicationParallelFFT(region *InRegion, int NumberOfInRegions,
			     region **OutRegion, int *NumberOfOutRegions,
			     int DomainDim[], int Rank,
			     int direction, int TransposeOnCompletion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

#ifdef ISOLATED_GRAVITY
int PrepareIsolatedGreensFunction(region *GreensFunction, int proc, 
				  int DomainDim[], TopGridData *MetaData);
#endif /* ISOLATED_GRAVITY */

#ifdef FAST_SIB
int ComputePotentialFieldLevelZeroPer(TopGridData *MetaData,
				      SiblingGridList SiblingList[],
				      HierarchyEntry *Grids[],
				      int NumberOfGrids);
#else
int ComputePotentialFieldLevelZeroPer(TopGridData *MetaData, 
				      HierarchyEntry *Grids[], 
				      int NumberOfGrids);
#endif




#ifdef FAST_SIB
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   SiblingGridList SiblingList[],
				   HierarchyEntry *Grids[],
				   int NumberOfGrids)
#else
int ComputePotentialFieldLevelZero(TopGridData *MetaData,
				   HierarchyEntry *Grids[], 
				   int NumberOfGrids)
#endif
{

  /* call the periodic solver */

#ifdef FAST_SIB
  if (ComputePotentialFieldLevelZeroPer(MetaData, SiblingList, Grids, NumberOfGrids) == FAIL) {
    fprintf(stderr, "Error in ComputePotentialFieldLevelZeroPer.\n");
    return FAIL;
  }   
#else
  if (ComputePotentialFieldLevelZeroPer(MetaData, Grids, NumberOfGrids) == FAIL) {
    fprintf(stderr, "Error in ComputePotentialFieldLevelZeroPer.\n");
    return FAIL;
  }
#endif
    
  return SUCCESS;
}




/******************************************************************/
/*  ComputePotentialFieldLevelZeroPer performs a root-grid        */
/*  potential field solver using periodic boundary conditions,    */
/*  via an FFT-based solution strategy.  This solver just calls   */
/*  the pre-existing code that Greg Bryan wrote.                  */
/******************************************************************/

#ifdef FAST_SIB
int ComputePotentialFieldLevelZeroPer(TopGridData *MetaData,
				      SiblingGridList SiblingList[],
				      HierarchyEntry *Grids[],
				      int NumberOfGrids)
#else
int ComputePotentialFieldLevelZeroPer(TopGridData *MetaData,
				      HierarchyEntry *Grids[], 
				      int NumberOfGrids)
#endif
{

  /* Static declarations (for Green's function). */
 
  static int FirstCall = TRUE, NumberOfGreensRegions;
  static region *GreensRegion;
 
  /* Declarations. */
 
  region *OutRegion = NULL;
  int NumberOfOutRegions, DomainDim[MAX_DIMENSION];
  int i, j, n, grid1, grid2, dim, TransposeOnCompletion;
 
  /* Allocate space for grid info. */
 
  int NumberOfRegions = NumberOfGrids;
  region *InitialRegion = new region[NumberOfRegions];
 
  /* Compute adot/a at time = t+1/2dt (time-centered). */
 
  FLOAT a = 1, dadt, MidTime = Grids[0]->GridData->ReturnTime() +
                           0.5*Grids[0]->GridData->ReturnTimeStep();
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionFactor(MidTime, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
    }
 
  /* If we are doing an isolated case then we do not have to tranpose the field
     back because the Greens function is computed in real space and then transformed,
     so the multiplication can be done on two un-transformed fields. This saves time. */

  if (MetaData->GravityBoundary == TopGridPeriodic)
    TransposeOnCompletion = TRUE;
  else
    TransposeOnCompletion = FALSE;
 
  /* ------------------------------------------------------------------- */
  /* If this is the first time this routine has been called, then generate
     the Green's function. */
 
  if (FirstCall) {
 
    if (MetaData->GravityBoundary == TopGridPeriodic) {
 
      /* Periodic -- Prepare in k-space. */
 
      NumberOfGreensRegions = NumberOfGrids;
      GreensRegion = new region[NumberOfGreensRegions];
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
	if (Grids[grid1]->GridData->PreparePeriodicGreensFunction(
					     &(GreensRegion[grid1])) == FAIL) {
	  fprintf(stderr, "Error in grid->PreparePeriodicGreensFunction.\n");
	  return FAIL;
	}
 
    } else {
 
#ifdef ISOLATED_GRAVITY
 
      region *TempRegion = new region[NumberOfProcessors];

      /* Generate Greens function in real space. */

      int proc;
      for (proc = 0; proc < NumberOfProcessors; proc++)
	if (PrepareIsolatedGreensFunction(&TempRegion[proc], proc, DomainDim,
					  MetaData) 
	    == FAIL) {
	  fprintf(stderr, "Error in PrepareIsolatedGreensFunction.\n");
	  return FAIL;
	}

      /* Forward FFT Greens function. */

      //      TransposeOnCompletion = FALSE;  // for isolated case we can skip transpose back
      if (CommunicationParallelFFT(TempRegion, NumberOfProcessors,
				   &GreensRegion, &NumberOfGreensRegions,
				   DomainDim, MetaData->TopGridRank,
				   FFT_FORWARD, TransposeOnCompletion) == FAIL) {
	fprintf(stderr, "Error in CommunicationParallelFFT.\n");
	return FAIL;
      }

      /* Clean up. */
      
      if (GreensRegion != TempRegion)
	delete [] TempRegion;

#else  /* ISOLATED_GRAVITY */

      fprintf(stderr, "You do not have isolated BC's turned on.  Recompile.\n");
      return FAIL;

#endif /* ISOLATED_GRAVITY */
 
    } // end: if (Periodic)
 
    FirstCall = FALSE;
 
  } // end: if (FirstCall)
 
  /* ------------------------------------------------------------------- */
  /* Generate FFT regions for density field. */
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (Grids[grid1]->GridData->PrepareFFT(&InitialRegion[grid1],
					  GRAVITATING_MASS_FIELD, DomainDim)
	== FAIL) {
      fprintf(stderr, "Error in grid->PrepareFFT.\n");
      return FAIL;
    }
 
#ifdef ISOLATED_GRAVITY
  /* If doing isolated BC's then double the domain size. */

  if (MetaData->GravityBoundary == TopGridIsolated) {
    for (dim = 0; dim < MetaData->TopGridRank; dim++)
      DomainDim[dim] *= 2;
    DomainDim[0] -= 2; /* correct for real-to-complex extra 2 */
  }
#endif /* ISOLATED_GRAVITY */

  /* Forward FFT density field. */
 
  if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
			       &OutRegion, &NumberOfOutRegions,
			       DomainDim, MetaData->TopGridRank,
#ifdef ISOLATED_GRAVITY
			       FFT_FORWARD, TransposeOnCompletion) == FAIL) {
#else  /* ISOLATED_GRAVITY */
			       FFT_FORWARD, TRUE) == FAIL) {
#endif /* ISOLATED_GRAVITY */
    fprintf(stderr, "Error in CommunicationParallelFFT.\n");
    return FAIL;
  }
 
  /* Quick error check. */
 
  if (NumberOfOutRegions != NumberOfGreensRegions) {
    fprintf(stderr, "OutRegion(%"ISYM") != GreensRegion(%"ISYM")\n", NumberOfOutRegions,
	    NumberOfGreensRegions);
    return FAIL;
  }
 
  /* Compute coefficient for Greens function. */
 
  float coef = GravitationalConstant/a;

  //  for (int dim = 0; dim < MetaData->TopGridRank; dim++)
  //    coef *= (DomainRightEdge[dim] - DomainLeftEdge[dim])/float(DomainDim[dim]);
			
  /* Multiply density by Green's function to get potential. */
 
  for (i = 0; i < NumberOfGreensRegions; i++)
    if (OutRegion[i].Data != NULL) {
      int size = OutRegion[i].RegionDim[0]*OutRegion[i].RegionDim[1]*
	         OutRegion[i].RegionDim[2];

#ifdef ISOLATED_GRAVITY 
      /* In periodic case, the Greens function is purely real. */

      if (MetaData->GravityBoundary == TopGridPeriodic) {
	for (n = 0, j = 0; j < size; j += 2, n++) {
	  OutRegion[i].Data[j  ] *= coef*GreensRegion[i].Data[n];
	  OutRegion[i].Data[j+1] *= coef*GreensRegion[i].Data[n];
	}
      }

      /* In the isolated case, we have compute the Greens function by transforming
	 a real function so we must do a proper complex multiplication (this is a
	 convolution with the real function defined in PrepareIsolatedGreensFunction). */

      if (MetaData->GravityBoundary == TopGridIsolated) {
	float real_part, imag_part;
	for (j = 0; j < size; j += 2) {
	  real_part = OutRegion[i].Data[j  ]*GreensRegion[i].Data[j  ] -
	              OutRegion[i].Data[j+1]*GreensRegion[i].Data[j+1];
	  imag_part = OutRegion[i].Data[j+1]*GreensRegion[i].Data[j  ] +
	              OutRegion[i].Data[j  ]*GreensRegion[i].Data[j+1];
	  OutRegion[i].Data[j  ] = real_part;
	  OutRegion[i].Data[j+1] = imag_part;
	}
      }

#else /* ISOLATED_GRAVITY */
      for (n = 0, j = 0; j < size; j += 2, n++) {
	OutRegion[i].Data[j  ] *= coef*GreensRegion[i].Data[n];
	OutRegion[i].Data[j+1] *= coef*GreensRegion[i].Data[n];
      }
#endif /* ISOLATED_GRAVITY */
    }
 
  /* Inverse FFT potential field. */
 
  if (CommunicationParallelFFT(InitialRegion, NumberOfRegions,
			       &OutRegion, &NumberOfOutRegions,
			       DomainDim, MetaData->TopGridRank,
#ifdef ISOLATED_GRAVITY
			       FFT_INVERSE, TransposeOnCompletion) == FAIL) {
#else  /* ISOLATED_GRAVITY */
			       FFT_INVERSE, TRUE) == FAIL) {
#endif /* ISOLATED_GRAVITY */
    fprintf(stderr, "Error in CommunicationParallelFFT.\n");
    return FAIL;
  }
 
  /* Copy Potential in active region into while grid. */
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (Grids[grid1]->GridData->FinishFFT(&InitialRegion[grid1], POTENTIAL_FIELD,
			       DomainDim) == FAIL) {
      fprintf(stderr, "Error in grid->FinishFFT.\n");
      return FAIL;
    }
 
  /* Update boundary regions of potential
     (first set BCTempL/R which are fluid BC's because that's the format
      that CheckForOverlap takes). */
 
  boundary_type BCTempLeft[MAX_DIMENSION], BCTempRight[MAX_DIMENSION];

  if (Grids[0]->GridData->ReturnGravityBoundaryType() == TopGridPeriodic) {
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      BCTempLeft[dim] = BCTempRight[dim] = periodic;
  } else {

#ifdef ISOLATED_GRAVITY
    for (int dim = 0; dim < MAX_DIMENSION; dim++)
      BCTempLeft[dim] = BCTempRight[dim] = reflecting; // doesn't matter as long as not periodic
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->SetIsolatedPotentialBoundary();
#else  /* ISOLATED_GRAVITY */
    fprintf(stderr, "recompile with isolated boundary conditions turned on!\n");
    return FAIL;
#endif /* ISOLATED_GRAVITY */
  }

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* Post the receives */

  CommunicationDirection = COMMUNICATION_POST_RECEIVE;
  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;

#ifdef FAST_SIB
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
      if (Grids[grid1]->GridData->CheckForOverlap(
				      SiblingList[grid1].GridList[grid2],
				      BCTempLeft,
				      BCTempRight,
				      &grid::CopyPotentialField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	return FAIL;
      }
#else
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				      BCTempLeft, BCTempRight,
     	                              &grid::CopyPotentialField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	return FAIL;
      }
#endif

#ifdef FORCE_MSG_PROGRESS 
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* Send the data */

  CommunicationDirection = COMMUNICATION_SEND;

#ifdef FAST_SIB
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    for (grid2 = 0; grid2 < SiblingList[grid1].NumberOfSiblings; grid2++)
      if (Grids[grid1]->GridData->CheckForOverlap(
				      SiblingList[grid1].GridList[grid2],
				      BCTempLeft,
				      BCTempRight,
				      &grid::CopyPotentialField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	return FAIL;
     }
#else
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    for (grid2 = 0; grid2 < NumberOfGrids; grid2++)
      if (Grids[grid1]->GridData->CheckForOverlap(Grids[grid2]->GridData,
				      BCTempLeft, BCTempRight,
     	                              &grid::CopyPotentialField) == FAIL) {
	fprintf(stderr, "Error in grid->CopyPotentialField.\n");
	return FAIL;
      }
#endif 

#ifdef FORCE_MSG_PROGRESS
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* Process the receives */

  if (CommunicationReceiveHandler() == FAIL)
    return FAIL;

  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
 
  /* Clean up. */
 
  delete [] InitialRegion;
  if (OutRegion != InitialRegion)
    delete [] OutRegion;
 
  if (CopyGravPotential)
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    {
      // fprintf(stderr, "Call CP from ComputePotentialFieldLevelZero\n");
      Grids[grid1]->GridData->CopyPotentialToBaryonField();
    }
 
  return SUCCESS;
}
