/***********************************************************************
/
/  COMMUNICATION ROUTINE: LOAD BALANCE ROOT GRIDS
/
/  written by: John Wise
/  date:       July, 2009
/  modified1: Michael Kuhlen
/             Rewrote ASCII input and added HDF5 input.
/  date:      October, 2010
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "communication.h"
#include "CommunicationUtilities.h"

#define LOAD_BALANCE_RATIO 1.05

int Enzo_Dims_create(int nnodes, int ndims, int *dims);
int DetermineNumberOfNodes(void);
int LoadBalanceSimulatedAnnealing(int NumberOfGrids, int NumberOfNodes, 
				  int* &ProcessorNumbers, double* NumberOfCells, 
				  double* NumberOfSubcells);
int LoadBalanceHilbertCurveRootGrids(FLOAT *GridCenters[], int *CellCount,
				     int NumberOfGrids, int* &RootProcessors);

int InitialLoadBalanceRootGrids(FILE *fptr, hid_t Hfile_id, int TopGridRank,
				int TopGridDim, int &NumberOfRootGrids,
				int* &RootProcessors)
{

  /* Determine number of nodes and thus cores/node */

  
  int NumberOfNodes;
  NumberOfNodes = DetermineNumberOfNodes();

  if (NumberOfProcessors == 1 || LoadBalancing == 0)
    return SUCCESS;

  /* If this is a zoom-in calculation, we shouldn't load balance the
     root grids by nodes but the finest static grid. */

//  if (StaticRefineRegionLevel[0] != INT_UNDEFINED)
//    return SUCCESS;

  /* Declarations */

  bool FinishedLevelZero;
  char line[MAX_LINE_LENGTH];
  char group_name[MAX_LINE_LENGTH];
  int i, j, dim, GridID, Rank, ThisLevel, dummy, GridDims[MAX_DIMENSION];
  int ThisTask, *GridMap;
  int *Workload;
  int size, grid_num, RootGridID, GridPosition[MAX_DIMENSION];
  int Layout[MAX_DIMENSION], LayoutTemp[MAX_DIMENSION];
  double *NumberOfCells, *NumberOfSubgridCells;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT *GridCenters[MAX_DIMENSION];
  float ThisCellWidth;

  int NumberOfDaughterGrids;
  int DaughterGridDims[MAX_DIMENSION];
  hid_t group_id, dset_id, attr_id;
  hid_t daughter_group_id;
  herr_t h5_status;
  herr_t h5_error = -1;

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Layout[dim] = 0;
    GridPosition[dim] = 0;
  }

  // Root processor does all of the work, and broadcasts the new
  // processor numbers

  if (MyProcessorNumber == ROOT_PROCESSOR) {

  // Count number of level-0 grids (assume that all level-0 grids at
  // the beginning)

    //    if (HierarchyFileInputFormat == 0) {
    if (HierarchyFileInputFormat % 2 == 0) {
      group_id = H5Gopen(Hfile_id, "/Level0");
      attr_id = H5Aopen_name(group_id, "NumberOfGrids");
      h5_status = H5Aread(attr_id, HDF5_INT, &NumberOfRootGrids);
      h5_status = H5Aclose(attr_id);
      h5_status = H5Gclose(group_id);
    } // ENDIF HDF5 input

    if (HierarchyFileInputFormat == 1) {
    
      NumberOfRootGrids = 0;
      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
	
	if (sscanf(line, "Grid = %"ISYM, &dummy) > 0)
	  NumberOfRootGrids++;
	
	// Reached the last root grid
	if (strstr(line, "NextGridNextLevel") != NULL)
	  break;

      } // ENDWHILE lines
      rewind(fptr);
    } // ENDIF ASCII input

  } // ENDIF root processor

  // If we're doing normal load balancing, synchronize and exit.
  if (LoadBalancing == 1) {
#ifdef USE_MPI
    MPI_Bcast(&NumberOfRootGrids, 1, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);
#endif /* USE_MPI */
    return SUCCESS;
  } // ENDIF LoadBalancing == 1

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    RootProcessors = new int[NumberOfRootGrids];
    
    if (LoadBalancing == 2 || LoadBalancing == 3) {
      // Compute (processor number)->(grid number) map
      Enzo_Dims_create(NumberOfRootGrids, TopGridRank, LayoutTemp);
      for (dim = 0; dim < TopGridRank; dim++)
	Layout[TopGridRank-1-dim] = LayoutTemp[dim];
      for (dim = TopGridRank; dim < MAX_DIMENSION; dim++)
	Layout[TopGridRank-1-dim] = 0;
      
      GridMap = new int[NumberOfRootGrids];
      NumberOfCells = new double[NumberOfRootGrids];
      NumberOfSubgridCells = new double[NumberOfRootGrids];
      for (i = 0; i < NumberOfRootGrids; i++) {
	NumberOfCells[i] = 0;
	NumberOfSubgridCells[i] = 0;
      }
    }
    
    if (LoadBalancing == 4) {
      Workload = new int[NumberOfRootGrids];
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	GridCenters[dim] = new FLOAT[NumberOfRootGrids];

      for (i = 0; i < NumberOfRootGrids; i++)
	Workload[i] = 0;
    }

    // Read the root grids from the hierarchy file
    fgets(line, MAX_LINE_LENGTH, fptr);
    for(i=0;i<NumberOfRootGrids;i++) {
      
      //    if (HierarchyFileInputFormat == 0) {
      if (HierarchyFileInputFormat % 2 == 0) {
	GridID = i+1;
	
	sprintf(group_name,"/Level0/Grid%"GROUP_TAG_FORMAT""ISYM, GridID);
	group_id = H5Gopen(Hfile_id, group_name);
	
	attr_id = H5Aopen_name(group_id, "Task");
	h5_status = H5Aread(attr_id, HDF5_INT, &ThisTask);
	h5_status = H5Aclose(attr_id);
	
	attr_id = H5Aopen_name(group_id, "GridRank");
	h5_status = H5Aread(attr_id, HDF5_INT, &Rank);
	h5_status = H5Aclose(attr_id);

	dset_id = H5Dopen(group_id, "GridDimension");
	h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) GridDims);
	h5_status = H5Dclose(dset_id);

	dset_id = H5Dopen(group_id, "GridLeftEdge");
	h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) LeftEdge);
	h5_status = H5Dclose(dset_id);

	dset_id = H5Dopen(group_id, "GridRightEdge");
	h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) RightEdge);
	h5_status = H5Dclose(dset_id);

	// if LoadBalancing == 2 or 3, get number of cells in daughter
	// grids
	if (LoadBalancing == 2 || LoadBalancing == 3) {
	  attr_id = H5Aopen_name(group_id, "NumberOfDaughterGrids");
	  h5_status = H5Aread(attr_id, HDF5_INT, &NumberOfDaughterGrids);
	  h5_status = H5Aclose(attr_id);
	  
	  for (j=0;j<NumberOfDaughterGrids;j++) {
	    sprintf(group_name,"DaughterGrids/DaughterGrid%"ISYM, j);
	    daughter_group_id = H5Gopen(group_id, group_name);

	    dset_id = H5Dopen(daughter_group_id, "GridDimensions");
	    h5_status = H5Dread(dset_id, HDF5_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) DaughterGridDims);
	    h5_status = H5Dclose(dset_id);

	    for (dim = 0, size = 1; dim < Rank; dim++)
	      size *= GridDims[dim];

	    NumberOfSubgridCells[i] += size;
	  }
	  h5_status = H5Gclose(daughter_group_id);
	}

	h5_status = H5Gclose(group_id);

      } // ENDIF HDF5 input

      if (HierarchyFileInputFormat == 1) {
	
	Rank = 1;
	line[0] = 0;
	while (line[0] != '\n') {
	  fgets(line, MAX_LINE_LENGTH, fptr);
	  if (feof(fptr)) break;
	  
	  sscanf(line, "Grid = %"ISYM, &GridID);
	  sscanf(line, "Task = %"ISYM, &ThisTask);
	  sscanf(line, "GridRank = %"ISYM, &Rank);
	  switch (Rank) {
	  case 1:
	    sscanf(line, "GridDimension = %"ISYM, GridDims);
	    sscanf(line, "GridLeftEdge = %"PSYM, LeftEdge);
	    sscanf(line, "GridRightEdge = %"PSYM, RightEdge);
	    break;
	  case 2:
	    sscanf(line, "GridDimension = %"ISYM" %"ISYM, GridDims, GridDims+1);
	    sscanf(line, "GridLeftEdge = %"PSYM" %"PSYM, LeftEdge, LeftEdge+1);
	    sscanf(line, "GridRightEdge = %"PSYM" %"PSYM, RightEdge, RightEdge+1);
	    break;
	  case 3:
	    sscanf(line, "GridDimension = %"ISYM" %"ISYM" %"ISYM,
		   GridDims, GridDims+1, GridDims+2);
	    sscanf(line, "GridLeftEdge = %"PSYM" %"PSYM" %"PSYM,
		   LeftEdge, LeftEdge+1, LeftEdge+2);
	    sscanf(line, "GridRightEdge = %"PSYM" %"PSYM" %"PSYM,
		   RightEdge, RightEdge+1, RightEdge+2);
	    break;
	  default:
	    ENZO_FAIL("Bad grid rank!");
	  } // ENDSWITCH Rank
	  
	} // ENDWHILE line[0] != '\n' (i.e. end of this grid's hierarchy information)
      } // ENDIF ASCII input
      
      // Compute level and other derived quantities
      for (dim = 0, size = 1; dim < Rank; dim++) {
	GridPosition[dim] =
	  int(Layout[dim] * (LeftEdge[dim] - DomainLeftEdge[dim]) /
	      (DomainRightEdge[dim] - DomainLeftEdge[dim]));
	
	if (LoadBalancing == 4) {
	  GridCenters[dim][GridID-1] =
	    (LeftEdge[dim] - DomainLeftEdge[dim]) /
	    (DomainRightEdge[dim] - DomainLeftEdge[dim]);
	}
	
	size *= GridDims[dim];
      }
      
      RootProcessors[GridID-1] = ThisTask;
      
      if (LoadBalancing == 2 || LoadBalancing == 3) {
	
	grid_num = GridPosition[0] +
	  Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
	
	GridMap[grid_num] = GridID-1;
	NumberOfCells[GridID-1] += size;

      } // ENDIF LoadBalancing == 2 or 3
      
      if (LoadBalancing == 4) {
	Workload[GridID-1] = size;	
	
	if (GridID == NumberOfRootGrids)
	  break;
      }  // ENDIF LoadBalancing == 4
      
    } // loop over root grids
   

    // If LoadBalancing == 2 or 3 and ASCII input, then we need to
    // continue reading the hierarchy file until we get all the
    // level=1 grids in order to get NumberOfSubgridCells
    // (Yuck!! -- All hail HDF5 hierarchy file!)

    if ( (HierarchyFileInputFormat == 1) &&
	 (LoadBalancing == 2 || LoadBalancing == 3) ) {

      while(-1) {
	// check for (line==NULL) in case there's no level>0 grids
	if (feof(fptr)) break;
	
	Rank = 1;
	line[0] = 0;
	while (line[0] != '\n') {
	  fgets(line, MAX_LINE_LENGTH, fptr);
	  if (feof(fptr)) break;
	  
	  sscanf(line, "Grid = %"ISYM, &GridID);
	  sscanf(line, "Task = %"ISYM, &ThisTask);
	  sscanf(line, "GridRank = %"ISYM, &Rank);
	  switch (Rank) {
	  case 1:
	    sscanf(line, "GridDimension = %"ISYM, GridDims);
	    sscanf(line, "GridLeftEdge = %"PSYM, LeftEdge);
	    sscanf(line, "GridRightEdge = %"PSYM, RightEdge);
	    break;
	  case 2:
	    sscanf(line, "GridDimension = %"ISYM" %"ISYM, GridDims, GridDims+1);
	    sscanf(line, "GridLeftEdge = %"PSYM" %"PSYM, LeftEdge, LeftEdge+1);
	    sscanf(line, "GridRightEdge = %"PSYM" %"PSYM, RightEdge, RightEdge+1);
	    break;
	  case 3:
	    sscanf(line, "GridDimension = %"ISYM" %"ISYM" %"ISYM,
		   GridDims, GridDims+1, GridDims+2);
	    sscanf(line, "GridLeftEdge = %"PSYM" %"PSYM" %"PSYM,
		   LeftEdge, LeftEdge+1, LeftEdge+2);
	    sscanf(line, "GridRightEdge = %"PSYM" %"PSYM" %"PSYM,
		   RightEdge, RightEdge+1, RightEdge+2);
	    break;
	  default:
	    ENZO_FAIL("Bad grid rank!");
	  } // ENDSWITCH Rank
	  
	} // ENDWHILE line[0] != '\n' (i.e. end of this grid's hierarchy information)

	for (dim = 0, size = 1; dim < Rank; dim++) {
	  GridPosition[dim] =
	    int(Layout[dim] * (LeftEdge[dim] - DomainLeftEdge[dim]) /
		(DomainRightEdge[dim] - DomainLeftEdge[dim]));
	  
	  size *= GridDims[dim];
	}
      
	grid_num = GridPosition[0] +
	  Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
	
	ThisCellWidth = (RightEdge[0] - LeftEdge[0]) /
	  (GridDims[0] - 2*NumberOfGhostZones);
	ThisLevel = nint(-logf(TopGridDim * ThisCellWidth) / M_LN2);
	
	if (ThisLevel == 1) {	  
	  RootGridID = GridMap[grid_num];
	  NumberOfSubgridCells[RootGridID] += size;
	}
      } // ENDWHILE(-1), i.e. loop over remaining grids

    } // ENDIF (ASCII input AND LoadBalancing == 2 or 3)
    
    
    // reset file pointer for ASCII input
    if (HierarchyFileInputFormat == 1)
      rewind(fptr);

    
    if (LoadBalancing == 2 || LoadBalancing == 3) {
      /* Now that we have the number of subgrid cells in each topgrid, we
	 can distribute the topgrid tiles so each compute node have a
	 similar total number of subgrid cells.  Use simulated annealing
	 to load balance. */

      LoadBalanceSimulatedAnnealing(NumberOfRootGrids, NumberOfNodes, 
				    RootProcessors, NumberOfCells,
				    NumberOfSubgridCells);
      delete [] GridMap;
      delete [] NumberOfCells;
      delete [] NumberOfSubgridCells;
    }      
    
    if (LoadBalancing == 4) {

      LoadBalanceHilbertCurveRootGrids(GridCenters, Workload, 
				       NumberOfRootGrids, RootProcessors);
      
      delete [] Workload;
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	delete [] GridCenters[dim];
    }  // ENDIF LoadBalancing == 4
    
  } // ENDIF ROOT_PROCESSOR

  
#ifdef USE_MPI
  MPI_Bcast(&NumberOfRootGrids, 1, IntDataType, ROOT_PROCESSOR, MPI_COMM_WORLD);
  if (NumberOfRootGrids > 1) {
    if (MyProcessorNumber != ROOT_PROCESSOR)
      RootProcessors = new int[NumberOfRootGrids];
    MPI_Bcast(RootProcessors, NumberOfRootGrids, IntDataType, ROOT_PROCESSOR, 
	      MPI_COMM_WORLD);
  } else {
    delete [] RootProcessors;
    RootProcessors = NULL;
  }
#endif

  return SUCCESS;
}
