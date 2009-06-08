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
/  This routine fills in ghost cells for Enzo gravity grids based on 
/  the gravity boundary conditions (periodic or not), using neighboring 
/  grid values.  The periodic/isolating BCs must be handled through the 
/  calling routine, which should set neighbor indices to MPI_PROC_NULL 
/  for non-communicating interfaces (i.e. isolating boundaries).
/
/  Note: we perform all x0 communication first, then all x1 
/  communication, then all x2 communication.  This ensures that values 
/  at edges and corners of the box are communicated appropriately.
/
/  written by: Daniel R. Reynolds
/  date:       November, 2005
/  modified1:  Robert Harkness, August 12th 2006
/              For I/O macros.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdio.h>
/* Original includes for Enzo */
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "error.h"


int GravityBdryExchange(float *enzovec, int *layout, int *location, 
			int GravityGhosts, int edim0, int edim1, 
			int edim2, int x0l, int x0r, int x1l, int x1r, 
			int x2l, int x2r)
{

  /* If MPI not used, just return */
#ifndef USE_MPI
  return SUCCESS;
#else

  /* some local variables */
  int i, j, k, idx;
  MPI_Arg myrank;
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;

  /* Get MPI processor rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  /* set buffer sizes */
  int x0buffsize = edim1*edim2*GravityGhosts;
  int x1buffsize = edim0*edim2*GravityGhosts;
  int x2buffsize = edim0*edim1*GravityGhosts;

  /* set exchange tags */
  int msg_xch_x0l = 10000;
  int msg_xch_x0r = 10001;
  int msg_xch_x1l = 10002;
  int msg_xch_x1r = 10003;
  int msg_xch_x2l = 10004;
  int msg_xch_x2r = 10005;

  /* allocate request IDs */
  MPI_Request id_recv_x0l, id_send_x0l;
  MPI_Request id_recv_x0r, id_send_x0r;
  MPI_Request id_recv_x1l, id_send_x1l;
  MPI_Request id_recv_x1r, id_send_x1r;
  MPI_Request id_recv_x2l, id_send_x2l;
  MPI_Request id_recv_x2r, id_send_x2r;
  
  /* allocate MPI status object */
  MPI_Status status;

  /* Create file exchange buffers for x0 communication */ 
  float *x0sendbuf = new float[x0buffsize];
  float *x0recvbuf = new float[x0buffsize];


  /* Update x0L boundaries */
  {
    /* open receive buffer for x0L boundary */
    if (x0l != MPI_PROC_NULL) {
      if (MPI_Irecv(x0recvbuf, x0buffsize, DataType, x0l, 
		    msg_xch_x0l, MPI_COMM_WORLD, &id_recv_x0l) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x0L receive error\n",myrank);
	return FAIL;
      }
    }
    
    /* fill and send buffer for x0R boundary */
    if (x0r != MPI_PROC_NULL) {
      for (idx=0, k=0; k<edim2; k++)
	for (j=0; j<edim1; j++)
	  for (i=edim0-2*GravityGhosts; i<edim0-GravityGhosts; i++, idx++)
	    x0sendbuf[idx] = enzovec[(k*edim1 + j)*edim0 + i];
      
      if (MPI_Isend(x0sendbuf, x0buffsize, DataType, x0r, 
		    msg_xch_x0l, MPI_COMM_WORLD, &id_send_x0l) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x0R send error\n",myrank);
	return FAIL;
      }
    }
    
    /* wait for x0L data, and update ghost cells */
    if (x0l != MPI_PROC_NULL) {
      if (MPI_Wait(&id_recv_x0l, &status) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x0L wait error\n",myrank);
	return FAIL;
      }
      
      for (idx=0, k=0; k<edim2; k++)
	for (j=0; j<edim1; j++)
	  for (i=0; i<GravityGhosts; i++, idx++)
	    enzovec[(k*edim1 + j)*edim0 + i] = x0recvbuf[idx];      
    }
  }


  /* Update x0R boundaries */
  {
    /* open receive buffer for x0R boundary */
    if (x0r != MPI_PROC_NULL) {
      if (MPI_Irecv(x0recvbuf, x0buffsize, DataType, x0r, 
		    msg_xch_x0r, MPI_COMM_WORLD, &id_recv_x0r) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x0R receive error\n",myrank);
	return FAIL;
      }
    }
    
    /* fill and send buffer for x0L boundary */
    if (x0l != MPI_PROC_NULL) {
      for (idx=0, k=0; k<edim2; k++)
	for (j=0; j<edim1; j++)
	  for (i=GravityGhosts; i<2*GravityGhosts; i++, idx++)
	    x0sendbuf[idx] = enzovec[(k*edim1 + j)*edim0 + i];
      
      if (MPI_Isend(x0sendbuf, x0buffsize, DataType, x0l, 
		    msg_xch_x0r, MPI_COMM_WORLD, &id_send_x0r) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x0L send error\n",myrank);
	return FAIL;
      }
    }
    
    /* wait for x0R data, and update ghost cells */
    if (x0r != MPI_PROC_NULL) {
      if (MPI_Wait(&id_recv_x0r, &status) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x0R wait error\n",myrank);
	return FAIL;
      }
      
      for (idx=0, k=0; k<edim2; k++)
	for (j=0; j<edim1; j++)
	  for (i=edim0-GravityGhosts; i<edim0; i++, idx++)
	    enzovec[(k*edim1 + j)*edim0 + i] = x0recvbuf[idx];      
    }
  }


  /* Delete file exchange buffers for x0 communication */ 
  delete [] x0sendbuf;
  delete [] x0recvbuf;

  /* Create file exchange buffers for x1 communication */ 
  float *x1sendbuf = new float[x1buffsize];
  float *x1recvbuf = new float[x1buffsize];


  /* Update x1L boundaries */
  {
    /* open receive buffer for x1L boundary */
    if (x1l != MPI_PROC_NULL) {
      if (MPI_Irecv(x1recvbuf, x1buffsize, DataType, x1l, 
		    msg_xch_x1l, MPI_COMM_WORLD, &id_recv_x1l) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x1L receive error\n",myrank);
	return FAIL;
      }
    }
    
    /* fill and send buffer for x1R boundary */
    if (x1r != MPI_PROC_NULL) {
      for (idx=0, k=0; k<edim2; k++)
	for (j=edim1-2*GravityGhosts; j<edim1-GravityGhosts; j++)
	  for (i=0; i<edim0; i++, idx++)
	    x1sendbuf[idx] = enzovec[(k*edim1 + j)*edim0 + i];
      
      if (MPI_Isend(x1sendbuf, x1buffsize, DataType, x1r, 
		    msg_xch_x1l, MPI_COMM_WORLD, &id_send_x1l) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x1R send error\n",myrank);
	return FAIL;
      }
    }
    
    /* wait for x1L data, and update ghost cells */
    if (x1l != MPI_PROC_NULL) {
      if (MPI_Wait(&id_recv_x1l, &status) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x1L wait error\n",myrank);
	return FAIL;
      }
      
      for (idx=0, k=0; k<edim2; k++)
	for (j=0; j<GravityGhosts; j++)
	  for (i=0; i<edim0; i++, idx++)
	    enzovec[(k*edim1 + j)*edim0 + i] = x1recvbuf[idx];      
    }
  }


  /* Update x1R boundaries */
  {
    /* open receive buffer for x1R boundary */
    if (x1r != MPI_PROC_NULL) {
      if (MPI_Irecv(x1recvbuf, x1buffsize, DataType, x1r, 
		    msg_xch_x1r, MPI_COMM_WORLD, &id_recv_x1r) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x1R receive error\n",myrank);
	return FAIL;
      }
    }
    
    /* fill and send buffer for x1L boundary */
    if (x1l != MPI_PROC_NULL) {
      for (idx=0, k=0; k<edim2; k++)
	for (j=GravityGhosts; j<2*GravityGhosts; j++)
	  for (i=0; i<edim0; i++, idx++)
	    x1sendbuf[idx] = enzovec[(k*edim1 + j)*edim0 + i];
      
      if (MPI_Isend(x1sendbuf, x1buffsize, DataType, x1l, 
		    msg_xch_x1r, MPI_COMM_WORLD, &id_send_x1r) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x1L send error\n",myrank);
	return FAIL;
      }
    }
    
    /* wait for x1R data, and update ghost cells */
    if (x1r != MPI_PROC_NULL) {
      if (MPI_Wait(&id_recv_x1r, &status) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x1R wait error\n",myrank);
	return FAIL;
      }
      
      for (idx=0, k=0; k<edim2; k++)
	for (j=edim1-GravityGhosts; j<edim1; j++)
	  for (i=0; i<edim0; i++, idx++)
	    enzovec[(k*edim1 + j)*edim0 + i] = x1recvbuf[idx];      
    }
  }


  /* Delete file exchange buffers for x1 communication */ 
  delete [] x1sendbuf;
  delete [] x1recvbuf;

  /* Create file exchange buffers for x2 communication */ 
  float *x2sendbuf = new float[x2buffsize];
  float *x2recvbuf = new float[x2buffsize];


  /* Update x2L boundaries */
  {
    /* open receive buffer for x2L boundary */
    if (x2l != MPI_PROC_NULL) {
      if (MPI_Irecv(x2recvbuf, x2buffsize, DataType, x2l, 
		    msg_xch_x2l, MPI_COMM_WORLD, &id_recv_x2l) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x2L receive error\n",myrank);
	return FAIL;
      }
    }
    
    /* fill and send buffer for x2R boundary */
    if (x2r != MPI_PROC_NULL) {
      k=edim2-2;
      for (idx=0, k=edim2-2*GravityGhosts; k<edim2-GravityGhosts; k++)
	for (j=0; j<edim1; j++)
	  for (i=0; i<edim0; i++, idx++)
	    x2sendbuf[idx] = enzovec[(k*edim1 + j)*edim0 + i];
      
      if (MPI_Isend(x2sendbuf, x2buffsize, DataType, x2r, 
		    msg_xch_x2l, MPI_COMM_WORLD, &id_send_x2l) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x2R send error\n",myrank);
	return FAIL;
      }
    }
    
    /* wait for x2L data, and update ghost cells */
    if (x2l != MPI_PROC_NULL) {
      if (MPI_Wait(&id_recv_x2l, &status) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x2L wait error\n",myrank);
	return FAIL;
      }
      
      for (idx=0, k=0; k<GravityGhosts; k++)
	for (j=0; j<edim1; j++)
	  for (i=0; i<edim0; i++, idx++)
	    enzovec[(k*edim1 + j)*edim0 + i] = x2recvbuf[idx];      
    }
  }


  /* Update x2R boundaries */
  {
    /* open receive buffer for x2R boundary */
    if (x2r != MPI_PROC_NULL) {
      if (MPI_Irecv(x2recvbuf, x2buffsize, DataType, x2r, 
		    msg_xch_x2r, MPI_COMM_WORLD, &id_recv_x2r) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x2R receive error\n",myrank);
	return FAIL;
      }
    }
    
    /* fill and send buffer for x2L boundary */
    if (x2l != MPI_PROC_NULL) {
      for (idx=0, k=GravityGhosts; k<2*GravityGhosts; k++)
	for (j=0; j<edim1; j++)
	  for (i=0; i<edim0; i++, idx++)
	    x2sendbuf[idx] = enzovec[(k*edim1 + j)*edim0 + i];
      
      if (MPI_Isend(x2sendbuf, x2buffsize, DataType, x2l, 
		    msg_xch_x2r, MPI_COMM_WORLD, &id_send_x2r) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x2L send error\n",myrank);
	return FAIL;
      }
    }
    
    /* wait for x2R data, and update ghost cells */
    if (x2r != MPI_PROC_NULL) {
      if (MPI_Wait(&id_recv_x2r, &status) != 0) {
	fprintf(stderr,"GravityBdryExchange p%"ISYM": x2R wait error\n",myrank);
	return FAIL;
      }
      
      for (idx=0, k=edim2-GravityGhosts; k<edim2; k++)
	for (j=0; j<edim1; j++)
	  for (i=0; i<edim0; i++, idx++)
	    enzovec[(k*edim1 + j)*edim0 + i] = x2recvbuf[idx];      
    }
  }

  /* Delete file exchange buffers for x2 communication */ 
  delete [] x2sendbuf;
  delete [] x2recvbuf;
  

  return SUCCESS;

#endif
}

/******************************************************************/
