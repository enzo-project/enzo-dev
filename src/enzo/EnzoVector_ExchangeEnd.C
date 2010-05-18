/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Daniel R. Reynolds
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  This routine fills in ghost cells for EnzoVectors grids based on 
/  neighbor information.
/
/  Note: we perform all x0 communication first, then all x1 
/  communication, then all x2 communication.  This ensures that values 
/  at edges and corners of the box are communicated appropriately.
/
/  written by: Daniel R. Reynolds
/  date:       May, 2006
/
************************************************************************/
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "EnzoVector.h"

//  Vector Boundary Communication Routine
//  (may be used for parallelism, or even for single-proc. periodic BCs)
int EnzoVector::exchange_end()
{
#ifdef USE_MPI

  // some local variables
  int i, j, k, l, idx;
  MPI_Arg myrank;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int x2len = Nx2 + Ng2l + Ng2r;

  // Get MPI processor rank
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // allocate MPI status object
  MPI_Status status;

  float *mydata=NULL;

  /////////////////////////
  // wait for receives to complete

  // wait for x0L data, and update ghost cells
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x0l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0L wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<Ng0l; i++)
	    mydata[(k*x1len + j)*x0len + i] = RecvBufX0l[idx++];
    }
  }
  // wait for x0R data, and update ghost cells
  if (Nbors[0][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x0r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0R wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=x0len-Ng0r; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = RecvBufX0r[idx++];
    }
  }      
  // wait for x1L data, and update ghost cells
  if (Nbors[1][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x1l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1L wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<x2len; k++)
	for (j=0; j<Ng1l; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = RecvBufX1l[idx++];
    }
  }
  // wait for x1R data, and update ghost cells
  if (Nbors[1][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x1r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1R wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<x2len; k++)
	for (j=x1len-Ng1r; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = RecvBufX1r[idx++];
    }
  }      
  // wait for x2L data, and update ghost cells
  if (Nbors[2][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x2l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2L wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<Ng2l; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = RecvBufX2l[idx++];
    }
  }
  // wait for x2R data, and update ghost cells
  if (Nbors[2][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x2r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2R wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=x2len-Ng2r; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = RecvBufX2r[idx++];
    }
  }


      
  /////////////////////////
  // wait for deliveries to complete

  // wait to ensure proper delivery of x0R data
  if (Nbors[0][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x0l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x0L data
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x0r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0R wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x1R data
  if (Nbors[1][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x1l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x1L data
  if (Nbors[1][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x1r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1R wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x2R data
  if (Nbors[2][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x2l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x2L data
  if (Nbors[2][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x2r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2R wait error\n",myrank);
      return FAIL;
    }
  }
#endif  // end if USE_MPI

  return SUCCESS;
}

/******************************************************************/



//  Vector Boundary Communication Routine (single component)
//  (may be used for parallelism, or even for single-proc. periodic BCs)
int EnzoVector::exchange_end_component(int ivar)
{
#ifdef USE_MPI

  // some local variables
  int i, j, k, l, idx;
  MPI_Arg myrank;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int x2len = Nx2 + Ng2l + Ng2r;

  // Get MPI processor rank
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // allocate MPI status object
  MPI_Status status;

  float *mydata = data[ivar];

  /////////////////////////
  // wait for receives to complete

  // wait for x0L data, and update ghost cells
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x0l_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0L wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (k=0; k<x2len; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<Ng0l; i++)
	  mydata[(k*x1len + j)*x0len + i] = RecvBufX0l_comp[idx++];
  }
  // wait for x0R data, and update ghost cells
  if (Nbors[0][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x0r_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0R wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (k=0; k<x2len; k++)
      for (j=0; j<x1len; j++)
	for (i=x0len-Ng0r; i<x0len; i++)
	  mydata[(k*x1len + j)*x0len + i] = RecvBufX0r_comp[idx++];
  }      
  // wait for x1L data, and update ghost cells
  if (Nbors[1][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x1l_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1L wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (k=0; k<x2len; k++)
      for (j=0; j<Ng1l; j++)
	for (i=0; i<x0len; i++)
	  mydata[(k*x1len + j)*x0len + i] = RecvBufX1l_comp[idx++];
  }
  // wait for x1R data, and update ghost cells
  if (Nbors[1][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x1r_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1R wait error\n",myrank);
      return FAIL;
    }
    idx=0;
    for (k=0; k<x2len; k++)
      for (j=x1len-Ng1r; j<x1len; j++)
	for (i=0; i<x0len; i++)
	  mydata[(k*x1len + j)*x0len + i] = RecvBufX1r_comp[idx++];
  }      
  // wait for x2L data, and update ghost cells
  if (Nbors[2][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x2l_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2L wait error\n",myrank);
      return FAIL;
    }
    idx=0; 
    for (k=0; k<Ng2l; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<x0len; i++)
	  mydata[(k*x1len + j)*x0len + i] = RecvBufX2l_comp[idx++];
  }
  // wait for x2R data, and update ghost cells
  if (Nbors[2][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_recv_x2r_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2R wait error\n",myrank);
      return FAIL;
    }
    idx=0; 
    for (k=x2len-Ng2r; k<x2len; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<x0len; i++)
	  mydata[(k*x1len + j)*x0len + i] = RecvBufX2r_comp[idx++];
  }


      
  /////////////////////////
  // wait for deliveries to complete

  // wait to ensure proper delivery of x0R data
  if (Nbors[0][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x0l_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x0L data
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x0r_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0R wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x1R data
  if (Nbors[1][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x1l_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x1L data
  if (Nbors[1][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x1r_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1R wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x2R data
  if (Nbors[2][1] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x2l_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x2L data
  if (Nbors[2][0] != MPI_PROC_NULL) {
    if (MPI_Wait((MPI_Request *) id_send_x2r_comp, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2R wait error\n",myrank);
      return FAIL;
    }
  }
#endif  // end if USE_MPI

  return SUCCESS;
}

/******************************************************************/
