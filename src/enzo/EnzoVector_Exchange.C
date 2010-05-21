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
int EnzoVector::exchange()
{
#ifdef USE_MPI

  //  printf("Entering EnzoVector::exchange\n");

  // some local variables
  int i, j, k, l, idx;
  MPI_Arg myrank;
  MPI_Arg one=1;
  MPI_Datatype FDataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Datatype IDataType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int x2len = Nx2 + Ng2l + Ng2r;

  // Get MPI processor rank
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // set exchange tags
  MPI_Arg msg_xch_x0l = 10000;
  MPI_Arg msg_xch_x0r = 10001;
  MPI_Arg msg_xch_x1l = 10002;
  MPI_Arg msg_xch_x1r = 10003;
  MPI_Arg msg_xch_x2l = 10004;
  MPI_Arg msg_xch_x2r = 10005;

  // allocate MPI status object
  MPI_Status status;

  // allocate request IDs
  MPI_Request recv_x0l, send_x0l;
  MPI_Request recv_x0r, send_x0r;
  MPI_Request recv_x1l, send_x1l;
  MPI_Request recv_x1r, send_x1r;
  MPI_Request recv_x2l, send_x2l;
  MPI_Request recv_x2r, send_x2r;


  /////////////////////////////////////////////////////////
  // If first time called for this vector, set up information for exchanges
  if (NborGhosts[0][0] == 0) {

    // Determine send buffer sizes (initialize to 0)
    NborGhosts[0][0] = 0;  NborGhosts[0][1] = 0;
    NborGhosts[1][0] = 0;  NborGhosts[1][1] = 0;
    NborGhosts[2][0] = 0;  NborGhosts[2][1] = 0;

//     //   open receive value for x0L proc ghosts
//     if (Nbors[0][0] != MPI_PROC_NULL)
//       if (MPI_Irecv(&NborGhosts[0][0], one, IDataType, 
// 		    MPI_Arg(Nbors[0][0]), msg_xch_x0l, 
// 		    MPI_COMM_WORLD, &recv_x0l) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x0L receive error\n",myrank);
// 	return FAIL;
//       }
//     //   open receive value for x0R proc ghosts
//     if (Nbors[0][1] != MPI_PROC_NULL)
//       if (MPI_Irecv(&NborGhosts[0][1], one, IDataType, 
// 		    MPI_Arg(Nbors[0][1]), msg_xch_x0r, 
// 		    MPI_COMM_WORLD, &recv_x0r) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x0R receive error\n",myrank);
// 	return FAIL;
//       }
//     //   open receive value for x1L proc ghosts
//     if (Nbors[1][0] != MPI_PROC_NULL)
//       if (MPI_Irecv(&NborGhosts[1][0], one, IDataType, 
// 		    MPI_Arg(Nbors[1][0]), msg_xch_x1l, 
// 		    MPI_COMM_WORLD, &recv_x1l) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x1L receive error\n",myrank);
// 	return FAIL;
//       }
//     //   open receive value for x1R proc ghosts
//     if (Nbors[1][1] != MPI_PROC_NULL)
//       if (MPI_Irecv(&NborGhosts[1][1], one, IDataType, 
// 		    MPI_Arg(Nbors[1][1]), msg_xch_x1r, 
// 		    MPI_COMM_WORLD, &recv_x1r) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x1R receive error\n",myrank);
// 	return FAIL;
//       }
//     //   open receive value for x2L proc ghosts
//     if (Nbors[2][0] != MPI_PROC_NULL)
//       if (MPI_Irecv(&NborGhosts[2][0], one, IDataType, 
// 		    MPI_Arg(Nbors[2][0]), msg_xch_x2l, 
// 		    MPI_COMM_WORLD, &recv_x2l) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x2L receive error\n",myrank);
// 	return FAIL;
//       }
//     //   open receive value for x2R proc ghosts
//     if (Nbors[2][1] != MPI_PROC_NULL)
//       if (MPI_Irecv(&NborGhosts[2][1], one, IDataType, 
// 		    MPI_Arg(Nbors[2][1]), msg_xch_x2r, 
// 		    MPI_COMM_WORLD, &recv_x2r) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x2R receive error\n",myrank);
// 	return FAIL;
//       }


//     //   fill and send value for x0R proc ghosts
//     if (Nbors[0][1] != MPI_PROC_NULL)
//       if (MPI_Isend(&Ng0r, one, IDataType, 
// 		    MPI_Arg(Nbors[0][1]), msg_xch_x0l, 
// 		    MPI_COMM_WORLD, &send_x0l) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x0R send error\n",myrank);
// 	return FAIL;
//       }
//     //   fill and send value for x0L proc ghosts
//     if (Nbors[0][0] != MPI_PROC_NULL)
//       if (MPI_Isend(&Ng0l, one, IDataType, 
// 		    MPI_Arg(Nbors[0][0]), msg_xch_x0r, 
// 		    MPI_COMM_WORLD, &send_x0r) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x0L send error\n",myrank);
// 	return FAIL;
//       }
//     //   fill and send value for x1R proc ghosts
//     if (Nbors[1][1] != MPI_PROC_NULL)
//       if (MPI_Isend(&Ng1r, one, IDataType, 
// 		    MPI_Arg(Nbors[1][1]), msg_xch_x1l, 
// 		    MPI_COMM_WORLD, &send_x1l) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x1R send error\n",myrank);
// 	return FAIL;
//       }
//     //   fill and send value for x1L proc ghosts
//     if (Nbors[1][0] != MPI_PROC_NULL)
//       if (MPI_Isend(&Ng1l, one, IDataType, 
// 		    MPI_Arg(Nbors[1][0]), msg_xch_x1r, 
// 		    MPI_COMM_WORLD, &send_x1r) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x1L send error\n",myrank);
// 	return FAIL;
//       }
//     //   fill and send value for x2R proc ghosts
//     if (Nbors[2][1] != MPI_PROC_NULL)
//       if (MPI_Isend(&Ng2r, one, IDataType, 
// 		    MPI_Arg(Nbors[2][1]), msg_xch_x2l, 
// 		    MPI_COMM_WORLD, &send_x2l) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x2R send error\n",myrank);
// 	return FAIL;
//       }
//     //   fill and send value for x2L proc ghosts
//     if (Nbors[2][0] != MPI_PROC_NULL)
//       if (MPI_Isend(&Ng2l, one, IDataType, 
// 		    MPI_Arg(Nbors[2][0]), msg_xch_x2r, 
// 		    MPI_COMM_WORLD, &send_x2r) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x2L send error\n",myrank);
// 	return FAIL;
//       }


//     //   wait for x0L proc ghosts
//     if (Nbors[0][0] != MPI_PROC_NULL)
//       if (MPI_Wait(&recv_x0l, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x0L wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for x0R proc ghosts
//     if (Nbors[0][1] != MPI_PROC_NULL)
//       if (MPI_Wait(&recv_x0r, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x0R wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for x1L proc ghosts
//     if (Nbors[1][0] != MPI_PROC_NULL)
//       if (MPI_Wait(&recv_x1l, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x1L wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for x1R proc ghosts
//     if (Nbors[1][1] != MPI_PROC_NULL)
//       if (MPI_Wait(&recv_x1r, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x1R wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for x2L proc ghosts
//     if (Nbors[2][0] != MPI_PROC_NULL)
//       if (MPI_Wait(&recv_x2l, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x2L wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for x2R proc ghosts
//     if (Nbors[2][1] != MPI_PROC_NULL)
//       if (MPI_Wait(&recv_x2r, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x2R wait error\n",myrank);
// 	return FAIL;
//       }


//     //   wait for proper delivery of x0R proc ghosts
//     if (Nbors[0][1] != MPI_PROC_NULL)
//       if (MPI_Wait(&send_x0l, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x0R wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for proper delivery of x0L proc ghosts
//     if (Nbors[0][0] != MPI_PROC_NULL)
//       if (MPI_Wait(&send_x0r, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x0L wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for proper delivery of x1R proc ghosts
//     if (Nbors[1][1] != MPI_PROC_NULL)
//       if (MPI_Wait(&send_x1l, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x1R wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for proper delivery of x1L proc ghosts
//     if (Nbors[1][0] != MPI_PROC_NULL)
//       if (MPI_Wait(&send_x1r, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x1L wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for proper delivery of x2R proc ghosts
//     if (Nbors[2][1] != MPI_PROC_NULL)
//       if (MPI_Wait(&send_x2l, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x2R wait error\n",myrank);
// 	return FAIL;
//       }
//     //   wait for proper delivery of x2L proc ghosts
//     if (Nbors[2][0] != MPI_PROC_NULL)
//       if (MPI_Wait(&send_x2r, &status) != 0) {
// 	fprintf(stderr,"EnzoVector_Exchange p%i: x2L wait error\n",myrank);
// 	return FAIL;
//       }
//     MPI_Barrier(MPI_COMM_WORLD);
    NborGhosts[0][0] = Ng0r;  NborGhosts[0][1] = Ng0l;
    NborGhosts[1][0] = Ng1r;  NborGhosts[1][1] = Ng1l;
    NborGhosts[2][0] = Ng2r;  NborGhosts[2][1] = Ng2l;

    /////////////////////////////////////////////////////////
    // set buffer sizes, allocate send/receive buffers
    x0Lsendsize = x1len*x2len*NborGhosts[0][0]*Nspecies;
    x0Rsendsize = x1len*x2len*NborGhosts[0][1]*Nspecies;
    x0Lrecvsize = x1len*x2len*Ng0l*Nspecies;
    x0Rrecvsize = x1len*x2len*Ng0r*Nspecies;
    
    x1Lsendsize = x0len*x2len*NborGhosts[1][0]*Nspecies;
    x1Rsendsize = x0len*x2len*NborGhosts[1][1]*Nspecies;
    x1Lrecvsize = x0len*x2len*Ng1l*Nspecies;
    x1Rrecvsize = x0len*x2len*Ng1r*Nspecies;
    
    x2Lsendsize = x0len*x1len*NborGhosts[2][0]*Nspecies;
    x2Rsendsize = x0len*x1len*NborGhosts[2][1]*Nspecies;
    x2Lrecvsize = x0len*x1len*Ng2l*Nspecies;
    x2Rrecvsize = x0len*x1len*Ng2r*Nspecies;
    
    SendBufX0l = new float[x0Lsendsize];
    SendBufX0r = new float[x0Rsendsize];

    SendBufX1l = new float[x1Lsendsize];
    SendBufX1r = new float[x1Rsendsize];

    SendBufX2l = new float[x2Lsendsize];
    SendBufX2r = new float[x2Rsendsize];

    RecvBufX0l = new float[x0Lrecvsize];
    RecvBufX0r = new float[x0Rrecvsize];

    RecvBufX1l = new float[x1Lrecvsize];
    RecvBufX1r = new float[x1Rrecvsize];

    RecvBufX2l = new float[x2Lrecvsize];
    RecvBufX2r = new float[x2Rrecvsize];

  }  // end if NborGhosts[0][0] == 0  (first call to Exchange)


  float *mydata=NULL;

  /////////////////////////
  // open receive buffer communication channels 

  // open receive buffer for x0L boundary
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX0l, x0Lrecvsize, FDataType, 
		  MPI_Arg(Nbors[0][0]), msg_xch_x0l, 
		  MPI_COMM_WORLD, &recv_x0l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0L receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x0R boundary
  if (Nbors[0][1] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX0r, x0Rrecvsize, FDataType, 
		  MPI_Arg(Nbors[0][1]), msg_xch_x0r, 
		  MPI_COMM_WORLD, &recv_x0r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0R receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x1L boundary
  if (Nbors[1][0] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX1l, x1Lrecvsize, FDataType, 
		  MPI_Arg(Nbors[1][0]), msg_xch_x1l, 
		  MPI_COMM_WORLD, &recv_x1l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1L receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x1R boundary
  if (Nbors[1][1] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX1r, x1Rrecvsize, FDataType, 
		  MPI_Arg(Nbors[1][1]), msg_xch_x1r, 
		  MPI_COMM_WORLD, &recv_x1r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1R receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x2L boundary
  if (Nbors[2][0] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX2l, x2Lrecvsize, FDataType, 
		  MPI_Arg(Nbors[2][0]), msg_xch_x2l, 
		  MPI_COMM_WORLD, &recv_x2l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2L receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x2R boundary
  if (Nbors[2][1] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX2r, x2Rrecvsize, FDataType, 
		  MPI_Arg(Nbors[2][1]), msg_xch_x2r, 
		  MPI_COMM_WORLD, &recv_x2r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2R receive error\n",myrank);
      return FAIL;
    }
  }
  

  /////////////////////////
  // fill and send boundary data
  
  // fill and send buffer for x0R boundary
  if (Nbors[0][1] != MPI_PROC_NULL) {
    idx = 0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=x0len-Ng0r-NborGhosts[0][1]; i<x0len-Ng0r; i++)
	    SendBufX0l[idx++] = mydata[(k*x1len + j)*x0len + i];
    }
    if (MPI_Isend(SendBufX0l, x0Rsendsize, FDataType, 
		  MPI_Arg(Nbors[0][1]), msg_xch_x0l, 
		  MPI_COMM_WORLD, &send_x0l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0R send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x0L boundary
  if (Nbors[0][0] != MPI_PROC_NULL) {
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=Ng0l; i<Ng0l+NborGhosts[0][0]; i++)
	    SendBufX0r[idx++] = mydata[(k*x1len + j)*x0len + i];
    }
    if (MPI_Isend(SendBufX0r, x0Lsendsize, FDataType, 
		  MPI_Arg(Nbors[0][0]), msg_xch_x0r, 
		  MPI_COMM_WORLD, &send_x0r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0L send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x1R boundary
  if (Nbors[1][1] != MPI_PROC_NULL) {
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<x2len; k++)
	for (j=x1len-Ng1r-NborGhosts[1][1]; j<x1len-Ng1r; j++)
	  for (i=0; i<x0len; i++)
	    SendBufX1l[idx++] = mydata[(k*x1len + j)*x0len + i];
    }
    if (MPI_Isend(SendBufX1l, x1Rsendsize, FDataType, 
		  MPI_Arg(Nbors[1][1]), msg_xch_x1l, 
		  MPI_COMM_WORLD, &send_x1l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1R send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x1L boundary
  if (Nbors[1][0] != MPI_PROC_NULL) {
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=0; k<x2len; k++)
	for (j=Ng1l; j<Ng1l+NborGhosts[1][0]; j++)
	  for (i=0; i<x0len; i++)
	    SendBufX1r[idx++] = mydata[(k*x1len + j)*x0len + i];
    }
    if (MPI_Isend(SendBufX1r, x1Lsendsize, FDataType, 
		  MPI_Arg(Nbors[1][0]), msg_xch_x1r, 
		  MPI_COMM_WORLD, &send_x1r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1L send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x2R boundary
  if (Nbors[2][1] != MPI_PROC_NULL) {
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=x2len-Ng2r-NborGhosts[2][1]; k<x2len-Ng2r; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    SendBufX2l[idx++] = mydata[(k*x1len + j)*x0len + i];
    }
    if (MPI_Isend(SendBufX2l, x2Rsendsize, FDataType, 
		  MPI_Arg(Nbors[2][1]), msg_xch_x2l, 
		  MPI_COMM_WORLD, &send_x2l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2R send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x2L boundary
  if (Nbors[2][0] != MPI_PROC_NULL) {
    idx=0;
    for (l=0; l<Nspecies; l++) {
      mydata = data[l];
      for (k=Ng2l; k<Ng2l+NborGhosts[2][0]; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    SendBufX2r[idx++] = mydata[(k*x1len + j)*x0len + i];
    }
    if (MPI_Isend(SendBufX2r, x2Lsendsize, FDataType, 
		  MPI_Arg(Nbors[2][0]), msg_xch_x2r, 
		  MPI_COMM_WORLD, &send_x2r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2L send error\n",myrank);
      return FAIL;
    }
  }
  


  /////////////////////////
  // wait for receives to complete

  // wait for x0L data, and update ghost cells
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Wait(&recv_x0l, &status) != 0) {
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
    if (MPI_Wait(&recv_x0r, &status) != 0) {
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
    if (MPI_Wait(&recv_x1l, &status) != 0) {
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
    if (MPI_Wait(&recv_x1r, &status) != 0) {
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
    if (MPI_Wait(&recv_x2l, &status) != 0) {
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
    if (MPI_Wait(&recv_x2r, &status) != 0) {
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
    if (MPI_Wait(&send_x0l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x0L data
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x0r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x0R wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x1R data
  if (Nbors[1][1] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x1l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x1L data
  if (Nbors[1][0] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x1r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x1R wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x2R data
  if (Nbors[2][1] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x2l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x2L data
  if (Nbors[2][0] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x2r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange p%i: x2R wait error\n",myrank);
      return FAIL;
    }
  }

  //  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("Exiting EnzoVector::exchange\n");

  return SUCCESS;

#else
  // If MPI not used, just copy one boundary to ghosts of other boundary

  // some local variables
  int i, j, k;
#ifndef MPI_INT
  int MPI_PROC_NULL = -3;
#endif
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int x2len = Nx2 + Ng2l + Ng2r;

  // first check that right and left procs are the same 
  // (should be either MPI_PROC_NULL or 0, but definitely the same)
  if (Nbors[0][0] != Nbors[0][1]) {
    fprintf(stderr,"EnzoVector::exchange error in x0 neighbors");
    return FAIL;
  }
  if (Nbors[1][0] != Nbors[1][1]) {
    fprintf(stderr,"EnzoVector::exchange error in x1 neighbors");
    return FAIL;
  }
  if (Nbors[2][0] != Nbors[2][1]) {
    fprintf(stderr,"EnzoVector::exchange error in x2 neighbors");
    return FAIL;
  }


  // iterate over all of the species, filling ghosts with opposite face data

  float *mydata;
  for (int idat=0; idat<Nspecies; idat++) {

    // extract this variable from the vector array
    mydata = data[idat];
    
    // Update x0L boundaries
    if (Nbors[0][0] != MPI_PROC_NULL) 
      for (k=0; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<Ng0l; i++)
	    mydata[(k*x1len + j)*x0len + i] = 
	      mydata[(k*x1len + j)*x0len + x0len-Ng0r-Ng0l + i];
    
    // Update x0R boundaries
    if (Nbors[0][1] != MPI_PROC_NULL) 
      for (k=0; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<Ng0r; i++)
	    mydata[(k*x1len + j)*x0len + x0len-Ng0r + i] = 
	      mydata[(k*x1len + j)*x0len + Ng0l + i];
    
    // Update x1L boundaries
    if (Nbors[1][0] != MPI_PROC_NULL) 
      for (k=0; k<x2len; k++)
	for (j=0; j<Ng1l; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = 
	      mydata[(k*x1len + x1len-Ng0r-Ng0r + j)*x0len + i];
    
    // Update x1R boundaries
    if (Nbors[1][1] != MPI_PROC_NULL) 
      for (k=0; k<x2len; k++)
	for (j=0; j<Ng1r; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + x1len-Ng1r + j)*x0len + i] = 
	      mydata[(k*x1len + Ng1l + j)*x0len + i];
    
    // Update x2L boundaries
    if (Nbors[2][0] != MPI_PROC_NULL) 
      for (k=0; k<Ng2l; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = 
	      mydata[((k + x2len-Ng2r-Ng2l)*x1len + j)*x0len + i];

    // Update x2R boundaries
    if (Nbors[2][1] != MPI_PROC_NULL) 
      for (k=0; k<Ng2r; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    mydata[((k+x2len-Ng2r)*x1len + j)*x0len + i] = 
	      mydata[((k+Ng2l)*x1len + j)*x0len + i];
    
  }  // end for idat

  return SUCCESS;

#endif  // end if USE_MPI
}

/******************************************************************/



//  Vector Boundary Communication Routine (one component only)
//  (may be used for parallelism, or even for single-proc. periodic BCs)
int EnzoVector::exchange_component(int ivar)
{
#ifdef USE_MPI

  //  printf("Entering EnzoVector::exchange_component\n");

  // some local variables
  int i, j, k, l, idx;
  MPI_Arg myrank;
  MPI_Arg one=1;
  MPI_Datatype FDataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Datatype IDataType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int x2len = Nx2 + Ng2l + Ng2r;

  // Get MPI processor rank
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // set exchange tags
  MPI_Arg msg_xch_x0l = 10000 + (ivar+1)*10;
  MPI_Arg msg_xch_x0r = 10001 + (ivar+1)*10;
  MPI_Arg msg_xch_x1l = 10002 + (ivar+1)*10;
  MPI_Arg msg_xch_x1r = 10003 + (ivar+1)*10;
  MPI_Arg msg_xch_x2l = 10004 + (ivar+1)*10;
  MPI_Arg msg_xch_x2r = 10005 + (ivar+1)*10;

  // allocate MPI status object
  MPI_Status status;

  // allocate request IDs
  MPI_Request recv_x0l, send_x0l;
  MPI_Request recv_x0r, send_x0r;
  MPI_Request recv_x1l, send_x1l;
  MPI_Request recv_x1r, send_x1r;
  MPI_Request recv_x2l, send_x2l;
  MPI_Request recv_x2r, send_x2r;

  // set send buffer sizes
  NborGhosts[0][0] = Ng0r;  NborGhosts[0][1] = Ng0l;
  NborGhosts[1][0] = Ng1r;  NborGhosts[1][1] = Ng1l;
  NborGhosts[2][0] = Ng2r;  NborGhosts[2][1] = Ng2l;
  
  x0Lsendsize_comp = x1len*x2len*NborGhosts[0][0];
  x0Rsendsize_comp = x1len*x2len*NborGhosts[0][1];
  x0Lrecvsize_comp = x1len*x2len*Ng0l;
  x0Rrecvsize_comp = x1len*x2len*Ng0r;
  
  x1Lsendsize_comp = x0len*x2len*NborGhosts[1][0];
  x1Rsendsize_comp = x0len*x2len*NborGhosts[1][1];
  x1Lrecvsize_comp = x0len*x2len*Ng1l;
  x1Rrecvsize_comp = x0len*x2len*Ng1r;
  
  x2Lsendsize_comp = x0len*x1len*NborGhosts[2][0];
  x2Rsendsize_comp = x0len*x1len*NborGhosts[2][1];
  x2Lrecvsize_comp = x0len*x1len*Ng2l;
  x2Rrecvsize_comp = x0len*x1len*Ng2r;
    
  /////////////////////////////////////////////////////////
  // If first time called for this vector, set up information for exchanges
  if (SendBufX0l_comp == NULL) {

    /////////////////////////////////////////////////////////
    // allocate send/receive buffers
    SendBufX0l_comp = new float[x0Lsendsize_comp];
    SendBufX0r_comp = new float[x0Rsendsize_comp];

    SendBufX1l_comp = new float[x1Lsendsize_comp];
    SendBufX1r_comp = new float[x1Rsendsize_comp];

    SendBufX2l_comp = new float[x2Lsendsize_comp];
    SendBufX2r_comp = new float[x2Rsendsize_comp];

    RecvBufX0l_comp = new float[x0Lrecvsize_comp];
    RecvBufX0r_comp = new float[x0Rrecvsize_comp];

    RecvBufX1l_comp = new float[x1Lrecvsize_comp];
    RecvBufX1r_comp = new float[x1Rrecvsize_comp];

    RecvBufX2l_comp = new float[x2Lrecvsize_comp];
    RecvBufX2r_comp = new float[x2Rrecvsize_comp];

  }  // end if NborGhosts[0][0] == 0  (first call to Exchange)


  float *mydata = data[ivar];

  /////////////////////////
  // open receive buffer communication channels 

  // open receive buffer for x0L boundary
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX0l_comp, x0Lrecvsize_comp, FDataType, 
		  MPI_Arg(Nbors[0][0]), msg_xch_x0l, 
		  MPI_COMM_WORLD, &recv_x0l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0L receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x0R boundary
  if (Nbors[0][1] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX0r_comp, x0Rrecvsize_comp, FDataType, 
		  MPI_Arg(Nbors[0][1]), msg_xch_x0r, 
		  MPI_COMM_WORLD, &recv_x0r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0R receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x1L boundary
  if (Nbors[1][0] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX1l_comp, x1Lrecvsize_comp, FDataType, 
		  MPI_Arg(Nbors[1][0]), msg_xch_x1l, 
		  MPI_COMM_WORLD, &recv_x1l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1L receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x1R boundary
  if (Nbors[1][1] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX1r_comp, x1Rrecvsize_comp, FDataType, 
		  MPI_Arg(Nbors[1][1]), msg_xch_x1r, 
		  MPI_COMM_WORLD, &recv_x1r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1R receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x2L boundary
  if (Nbors[2][0] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX2l_comp, x2Lrecvsize_comp, FDataType, 
		  MPI_Arg(Nbors[2][0]), msg_xch_x2l, 
		  MPI_COMM_WORLD, &recv_x2l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2L receive error\n",myrank);
      return FAIL;
    }
  }
  // open receive buffer for x2R boundary
  if (Nbors[2][1] != MPI_PROC_NULL) {
    if (MPI_Irecv(RecvBufX2r_comp, x2Rrecvsize_comp, FDataType, 
		  MPI_Arg(Nbors[2][1]), msg_xch_x2r, 
		  MPI_COMM_WORLD, &recv_x2r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2R receive error\n",myrank);
      return FAIL;
    }
  }
  

  /////////////////////////
  // fill and send boundary data
  
  // fill and send buffer for x0R boundary
  if (Nbors[0][1] != MPI_PROC_NULL) {
    idx=0;
    for (k=0; k<x2len; k++)
      for (j=0; j<x1len; j++)
	for (i=x0len-Ng0r-NborGhosts[0][1]; i<x0len-Ng0r; i++)
	  SendBufX0l_comp[idx++] = mydata[(k*x1len + j)*x0len + i];
    if (MPI_Isend(SendBufX0l_comp, x0Rsendsize_comp, FDataType, 
		  MPI_Arg(Nbors[0][1]), msg_xch_x0l, 
		  MPI_COMM_WORLD, &send_x0l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0R send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x0L boundary
  if (Nbors[0][0] != MPI_PROC_NULL) {
    idx=0; 
    for (k=0; k<x2len; k++)
      for (j=0; j<x1len; j++)
	for (i=Ng0l; i<Ng0l+NborGhosts[0][0]; i++)
	  SendBufX0r_comp[idx++] = mydata[(k*x1len + j)*x0len + i];
    if (MPI_Isend(SendBufX0r_comp, x0Lsendsize_comp, FDataType, 
		  MPI_Arg(Nbors[0][0]), msg_xch_x0r, 
		  MPI_COMM_WORLD, &send_x0r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0L send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x1R boundary
  if (Nbors[1][1] != MPI_PROC_NULL) {
    idx=0; 
    for (k=0; k<x2len; k++)
      for (j=x1len-Ng1r-NborGhosts[1][1]; j<x1len-Ng1r; j++)
	for (i=0; i<x0len; i++)
	  SendBufX1l_comp[idx++] = mydata[(k*x1len + j)*x0len + i];
    if (MPI_Isend(SendBufX1l_comp, x1Rsendsize_comp, FDataType, 
		  MPI_Arg(Nbors[1][1]), msg_xch_x1l, 
		  MPI_COMM_WORLD, &send_x1l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1R send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x1L boundary
  if (Nbors[1][0] != MPI_PROC_NULL) {
    idx=0; 
    for (k=0; k<x2len; k++)
      for (j=Ng1l; j<Ng1l+NborGhosts[1][0]; j++)
	for (i=0; i<x0len; i++)
	  SendBufX1r_comp[idx++] = mydata[(k*x1len + j)*x0len + i];
    if (MPI_Isend(SendBufX1r_comp, x1Lsendsize_comp, FDataType, 
		  MPI_Arg(Nbors[1][0]), msg_xch_x1r, 
		  MPI_COMM_WORLD, &send_x1r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1L send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x2R boundary
  if (Nbors[2][1] != MPI_PROC_NULL) {
    idx=0; 
    for (k=x2len-Ng2r-NborGhosts[2][1]; k<x2len-Ng2r; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<x0len; i++)
	  SendBufX2l_comp[idx++] = mydata[(k*x1len + j)*x0len + i];
    if (MPI_Isend(SendBufX2l_comp, x2Rsendsize_comp, FDataType, 
		  MPI_Arg(Nbors[2][1]), msg_xch_x2l, 
		  MPI_COMM_WORLD, &send_x2l) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2R send error\n",myrank);
      return FAIL;
    }
  }
  // fill and send buffer for x2L boundary
  if (Nbors[2][0] != MPI_PROC_NULL) {
    idx=0; 
    for (k=Ng2l; k<Ng2l+NborGhosts[2][0]; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<x0len; i++)
	  SendBufX2r_comp[idx++] = mydata[(k*x1len + j)*x0len + i];
    if (MPI_Isend(SendBufX2r_comp, x2Lsendsize_comp, FDataType, 
		  MPI_Arg(Nbors[2][0]), msg_xch_x2r, 
		  MPI_COMM_WORLD, &send_x2r) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2L send error\n",myrank);
      return FAIL;
    }
  }
  


  /////////////////////////
  // wait for receives to complete

  // wait for x0L data, and update ghost cells
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Wait(&recv_x0l, &status) != 0) {
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
    if (MPI_Wait(&recv_x0r, &status) != 0) {
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
    if (MPI_Wait(&recv_x1l, &status) != 0) {
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
    if (MPI_Wait(&recv_x1r, &status) != 0) {
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
    if (MPI_Wait(&recv_x2l, &status) != 0) {
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
    if (MPI_Wait(&recv_x2r, &status) != 0) {
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
    if (MPI_Wait(&send_x0l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x0L data
  if (Nbors[0][0] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x0r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x0R wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x1R data
  if (Nbors[1][1] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x1l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x1L data
  if (Nbors[1][0] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x1r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x1R wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x2R data
  if (Nbors[2][1] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x2l, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2L wait error\n",myrank);
      return FAIL;
    }
  }
  // wait to ensure proper delivery of x2L data
  if (Nbors[2][0] != MPI_PROC_NULL) {
    if (MPI_Wait(&send_x2r, &status) != 0) {
      fprintf(stderr,"EnzoVector_Exchange_Component p%i: x2R wait error\n",myrank);
      return FAIL;
    }
  }

  //  MPI_Barrier(MPI_COMM_WORLD);
  //  printf("Exiting EnzoVector::exchange_component\n");

  return SUCCESS;


#else
  // If MPI not used, just copy one boundary to ghosts of other boundary

  // some local variables
  int i, j, k;
#ifndef MPI_INT
  int MPI_PROC_NULL = -3;
#endif
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int x2len = Nx2 + Ng2l + Ng2r;

  // first check that right and left procs are the same 
  // (should be either MPI_PROC_NULL or 0, but definitely the same)
  if (Nbors[0][0] != Nbors[0][1]) {
    fprintf(stderr,"EnzoVector::exchange_component error in x0 neighbors");
    return FAIL;
  }
  if (Nbors[1][0] != Nbors[1][1]) {
    fprintf(stderr,"EnzoVector::exchange_component error in x1 neighbors");
    return FAIL;
  }
  if (Nbors[2][0] != Nbors[2][1]) {
    fprintf(stderr,"EnzoVector::exchange_component error in x2 neighbors");
    return FAIL;
  }

  // extract this variable from the vector array; fill ghosts with opposite face data
  float *mydata = data[ivar];
    
  // Update x0L boundaries
  if (Nbors[0][0] != MPI_PROC_NULL) 
    for (k=0; k<x2len; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<Ng0l; i++)
	  mydata[(k*x1len + j)*x0len + i] = 
	    mydata[(k*x1len + j)*x0len + x0len-Ng0r-Ng0l + i];
  
  // Update x0R boundaries
  if (Nbors[0][1] != MPI_PROC_NULL) 
    for (k=0; k<x2len; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<Ng0r; i++)
	  mydata[(k*x1len + j)*x0len + x0len-Ng0r + i] = 
	    mydata[(k*x1len + j)*x0len + Ng0l + i];
  
  // Update x1L boundaries
  if (Nbors[1][0] != MPI_PROC_NULL) 
    for (k=0; k<x2len; k++)
      for (j=0; j<Ng1l; j++)
	for (i=0; i<x0len; i++)
	  mydata[(k*x1len + j)*x0len + i] = 
	    mydata[(k*x1len + x1len-Ng0r-Ng0r + j)*x0len + i];
  
  // Update x1R boundaries
  if (Nbors[1][1] != MPI_PROC_NULL) 
    for (k=0; k<x2len; k++)
      for (j=0; j<Ng1r; j++)
	for (i=0; i<x0len; i++)
	  mydata[(k*x1len + x1len-Ng1r + j)*x0len + i] = 
	    mydata[(k*x1len + Ng1l + j)*x0len + i];
  
  // Update x2L boundaries
  if (Nbors[2][0] != MPI_PROC_NULL) 
    for (k=0; k<Ng2l; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<x0len; i++)
	  mydata[(k*x1len + j)*x0len + i] = 
	    mydata[((k + x2len-Ng2r-Ng2l)*x1len + j)*x0len + i];
  
  // Update x2R boundaries
  if (Nbors[2][1] != MPI_PROC_NULL) 
    for (k=0; k<Ng2r; k++)
      for (j=0; j<x1len; j++)
	for (i=0; i<x0len; i++)
	  mydata[((k+x2len-Ng2r)*x1len + j)*x0len + i] = 
	    mydata[((k+Ng2l)*x1len + j)*x0len + i];
  
  return SUCCESS;

#endif  // end if USE_MPI
}

/******************************************************************/
