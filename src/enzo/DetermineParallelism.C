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
/  This routine examines the top grid layout to determine the MPI 
/  processor parallelism information:  overall processor layout, 
/  individual processor location, neighbor information.
/
/  written by: Daniel R. Reynolds
/  date:       April, 2006
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdio.h>
#include <map>
#include <string>
#include <math.h>

/* Original includes for Enzo */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

int DetermineParallelism(HierarchyEntry *TopGrid, TopGridData &MetaData)
{

  /* Determine grid corresponding to this process from the Hierarchy */ 
  HierarchyEntry *ThisGrid = TopGrid;
  int i, j, foundgrid=0;
  for (i=0; i<=MAX_NUMBER_OF_SUBGRIDS; i++) {
    if (MyProcessorNumber != ThisGrid->GridData->ReturnProcessorNumber()) 
      ThisGrid = ThisGrid->NextGridThisLevel;
    else {foundgrid=1; break;}
  }
  if (foundgrid == 0) {
    printf("Error: proc %"ISYM" could not locate his grid\n",MyProcessorNumber);
    return FAIL;
  }

  /* If MPI not used, just set single-processor information and return */
#ifndef USE_MPI
  ThisGrid->GridData->SetProcessorLayout(0,1);
  ThisGrid->GridData->SetProcessorLayout(1,1);
  ThisGrid->GridData->SetProcessorLayout(2,1);
  ThisGrid->GridData->SetProcessorLocation(0,0);
  ThisGrid->GridData->SetProcessorLocation(1,0);
  ThisGrid->GridData->SetProcessorLocation(2,0);
  ThisGrid->GridData->SetProcessorNeighbors(0,0,0);
  ThisGrid->GridData->SetProcessorNeighbors(0,1,0);
  ThisGrid->GridData->SetProcessorNeighbors(1,0,0);
  ThisGrid->GridData->SetProcessorNeighbors(1,1,0);
  ThisGrid->GridData->SetProcessorNeighbors(2,0,0);
  ThisGrid->GridData->SetProcessorNeighbors(2,1,0);
  
  return SUCCESS;
#else

  /* Set up temporary variables/arrays */
  FLOAT MyLeftEdge[3]  = {ThisGrid->GridData->GetGridLeftEdge(0),
			  ThisGrid->GridData->GetGridLeftEdge(1),
			  ThisGrid->GridData->GetGridLeftEdge(2)};
  FLOAT MyRightEdge[3] = {ThisGrid->GridData->GetGridRightEdge(0),
			  ThisGrid->GridData->GetGridRightEdge(1),
			  ThisGrid->GridData->GetGridRightEdge(2)};
  FLOAT tol = 1.0e-12;

  /* Create file exchange buffers for x2 communication */ 
  FLOAT *X0ProcMesh = new FLOAT[MAX_NUMBER_OF_TASKS];
  FLOAT *X1ProcMesh = new FLOAT[MAX_NUMBER_OF_TASKS];
  FLOAT *X2ProcMesh = new FLOAT[MAX_NUMBER_OF_TASKS];

  /* Set my processor edges into the mesh */
  X0ProcMesh[0] = MyLeftEdge[0];  
  X1ProcMesh[0] = MyLeftEdge[1];
  X2ProcMesh[0] = MyLeftEdge[2];
  int x0procs=1, x1procs=1, x2procs=1;

  /* Set periodicity-related faces */
  FLOAT MyNeighborFace[3][2];
  if (fabs(MyLeftEdge[0] - DomainLeftEdge[0]) < tol) 
    MyNeighborFace[0][0] = DomainRightEdge[0];
  else
    MyNeighborFace[0][0] = MyLeftEdge[0];
  if (fabs(MyRightEdge[0] - DomainRightEdge[0]) < tol) 
    MyNeighborFace[0][1] = DomainLeftEdge[0];
  else
    MyNeighborFace[0][1] = MyRightEdge[0];
  if (fabs(MyLeftEdge[1] - DomainLeftEdge[1]) < tol) 
    MyNeighborFace[1][0] = DomainRightEdge[1];
  else
    MyNeighborFace[1][0] = MyLeftEdge[1];
  if (fabs(MyRightEdge[1] - DomainRightEdge[1]) < tol) 
    MyNeighborFace[1][1] = DomainLeftEdge[1];
  else
    MyNeighborFace[1][1] = MyRightEdge[1];
  if (fabs(MyLeftEdge[2] - DomainLeftEdge[2]) < tol) 
    MyNeighborFace[2][0] = DomainRightEdge[2];
  else
    MyNeighborFace[2][0] = MyLeftEdge[2];
  if (fabs(MyRightEdge[2] - DomainRightEdge[2]) < tol) 
    MyNeighborFace[2][1] = DomainLeftEdge[2];
  else
    MyNeighborFace[2][1] = MyRightEdge[2];

  /* Iterate over grids, storing information and finding neighbors */ 
  MPI_Arg NProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
  FLOAT x0Left, x0Right, x1Left, x1Right, x2Left, x2Right;
  bool NotEvenOnce;
  HierarchyEntry *ThatGrid = TopGrid;
  int ThatProc;
  for (i=0; i<NProcs; i++) {

    /* get processor number for that grid */
    ThatProc = ThatGrid->GridData->ReturnProcessorNumber();
    
    /* get bounds from that grid */
    x0Left  = ThatGrid->GridData->GetGridLeftEdge(0);
    x0Right = ThatGrid->GridData->GetGridRightEdge(0);
    x1Left  = ThatGrid->GridData->GetGridLeftEdge(1);
    x1Right = ThatGrid->GridData->GetGridRightEdge(1);
    x2Left  = ThatGrid->GridData->GetGridLeftEdge(2);
    x2Right = ThatGrid->GridData->GetGridRightEdge(2);


    /* check for left x0 neighbor */
    if ((fabs(MyNeighborFace[0][0]-x0Right) < tol) 
	&& (fabs(MyLeftEdge[1]-x1Left) < tol)
	&& (fabs(MyLeftEdge[2]-x2Left) < tol) ) {
      ThisGrid->GridData->SetProcessorNeighbors(0, 0, ThatProc);
    }
      
    /* check for right x0 neighbor */
    if ((fabs(MyNeighborFace[0][1]-x0Left) < tol) 
	&& (fabs(MyLeftEdge[1]-x1Left) < tol)
	&& (fabs(MyLeftEdge[2]-x2Left) < tol) ) {
      ThisGrid->GridData->SetProcessorNeighbors(0, 1, ThatProc);
    }
      
    /* check for left x1 neighbor */
    if ((fabs(MyNeighborFace[1][0]-x1Right) < tol) 
	&& (fabs(MyLeftEdge[0]-x0Left) < tol)
	&& (fabs(MyLeftEdge[2]-x2Left) < tol) ) {
      ThisGrid->GridData->SetProcessorNeighbors(1, 0, ThatProc);
    }
      
    /* check for right x1 neighbor */
    if ((fabs(MyNeighborFace[1][1]-x1Left) < tol) 
	&& (fabs(MyLeftEdge[0]-x0Left) < tol)
	&& (fabs(MyLeftEdge[2]-x2Left) < tol) ) {
      ThisGrid->GridData->SetProcessorNeighbors(1, 1, ThatProc);
    }
      
    /* check for left x2 neighbor */
    if ((fabs(MyNeighborFace[2][0]-x2Right) < tol) 
	&& (fabs(MyLeftEdge[0]-x0Left) < tol)
	&& (fabs(MyLeftEdge[1]-x1Left) < tol) ) {
      ThisGrid->GridData->SetProcessorNeighbors(2, 0, ThatProc);
    }
      
    /* check for right x2 neighbor */
    if ((fabs(MyNeighborFace[2][1]-x2Left) < tol) 
	&& (fabs(MyLeftEdge[0]-x0Left) < tol)
	&& (fabs(MyLeftEdge[1]-x1Left) < tol) ) {
      ThisGrid->GridData->SetProcessorNeighbors(2, 1, ThatProc);
    }

    /* insert grid left edges if they have not yet been added */
    /*    x0 direction */
    NotEvenOnce = true;
    for (j=0; j<x0procs; j++)
      if (fabs(x0Left-X0ProcMesh[j]) < tol) {
	NotEvenOnce = false;
	break;
      }
    if (NotEvenOnce) {
      X0ProcMesh[x0procs] = x0Left;
      x0procs++;
    }

    /*    x1 direction */
    NotEvenOnce = true;
    for (j=0; j<x1procs; j++)
      if (fabs(x1Left-X1ProcMesh[j]) < tol) {
	NotEvenOnce = false;
	break;
      }
    if (NotEvenOnce) {
      X1ProcMesh[x1procs] = x1Left;
      x1procs++;
    }

    /*    x2 direction */
    NotEvenOnce = true;
    for (j=0; j<x2procs; j++)
      if (fabs(x2Left-X2ProcMesh[j]) < tol) {
	NotEvenOnce = false;
	break;
      }
    if (NotEvenOnce) {
      X2ProcMesh[x2procs] = x2Left;
      x2procs++;
    }

    /* move on to next grid */
    ThatGrid = ThatGrid->NextGridThisLevel;    
  }

  /* Store the resulting processor layout */
  ThisGrid->GridData->SetProcessorLayout(0, x0procs);
  ThisGrid->GridData->SetProcessorLayout(1, x1procs);
  ThisGrid->GridData->SetProcessorLayout(2, x2procs);

  /* Determine this processor's location in the overall proc grid */
  /*    x0 direction */
  int location=0;
  for (j=0; j<x0procs; j++) 
    if ((MyLeftEdge[0]-X0ProcMesh[j]) > tol) location++;
  ThisGrid->GridData->SetProcessorLocation(0, location);

  /*    x1 direction */
  location=0;
  for (j=0; j<x1procs; j++) 
    if ((MyLeftEdge[1]-X1ProcMesh[j]) > tol) location++;
  ThisGrid->GridData->SetProcessorLocation(1, location);

  /*    x2 direction */
  location=0;
  for (j=0; j<x2procs; j++) 
    if ((MyLeftEdge[2]-X2ProcMesh[j]) > tol) location++;
  ThisGrid->GridData->SetProcessorLocation(2, location);


  /* Clean up */
  delete [] X0ProcMesh;
  delete [] X1ProcMesh;
  delete [] X2ProcMesh;


  /* Return */
  return SUCCESS;

#endif
}

/******************************************************************/
