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
/  Enzo Vector operations
/
/  written by: Daniel Reynolds
/  date:       May, 2006
/  modified1:  
/
/  PURPOSE: 
/
************************************************************************/

#ifdef USE_MPI
#include <mpi.h>
#else
typedef int MPI_Request;
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "EnzoVector.h"



//  Vector Constructor
EnzoVector::EnzoVector(int N0, int N1, int N2, int G0l, int G0r, 
		       int G1l, int G1r, int G2l, int G2r, int Ns,
		       int NBx0L, int NBx0R, int NBx1L, 
		       int NBx1R, int NBx2L, int NBx2R)
{
  // set internal length specifications
  Nx0 = N0;
  Nx1 = N1;
  Nx2 = N2;
  Ng0l = G0l;
  Ng0r = G0r;
  Ng1l = G1l;
  Ng1r = G1r;
  Ng2l = G2l;
  Ng2r = G2r;
  Nspecies = Ns;
  int Nlocal = Nx0*Nx1*Nx2*Nspecies;
#ifdef USE_MPI
  if ((NBx0L == MPI_PROC_NULL) && (NBx0R == MPI_PROC_NULL) &&
      (NBx1L == MPI_PROC_NULL) && (NBx1R == MPI_PROC_NULL) && 
      (NBx2L == MPI_PROC_NULL) && (NBx2R == MPI_PROC_NULL)) {
    Nglobal = Nlocal;
  }
  else {
    MPI_Arg one = 1;
    MPI_Datatype DataType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
    MPI_Allreduce(&Nlocal, &Nglobal, one, DataType, MPI_SUM, MPI_COMM_WORLD);  
  }
#else
  Nglobal = Nlocal;
#endif

  // set parallelism information on neighbors, other info to defaults
  Nbors[0][0] = NBx0L;  Nbors[0][1] = NBx0R;
  Nbors[1][0] = NBx1L;  Nbors[1][1] = NBx1R;
  Nbors[2][0] = NBx2L;  Nbors[2][1] = NBx2R;
  NborGhosts[0][0] = 0;  NborGhosts[0][1] = 0;
  NborGhosts[1][0] = 0;  NborGhosts[1][1] = 0;
  NborGhosts[2][0] = 0;  NborGhosts[2][1] = 0;
  x0Lsendsize = x0Rsendsize = x0Lrecvsize = x0Rrecvsize = 0;
  x1Lsendsize = x1Rsendsize = x1Lrecvsize = x1Rrecvsize = 0;
  x2Lsendsize = x2Rsendsize = x2Lrecvsize = x2Rrecvsize = 0;
  SendBufX0l=NULL; RecvBufX0l=NULL; SendBufX0r=NULL; RecvBufX0r=NULL;
  SendBufX1l=NULL; RecvBufX1l=NULL; SendBufX1r=NULL; RecvBufX1r=NULL;
  SendBufX2l=NULL; RecvBufX2l=NULL; SendBufX2r=NULL; RecvBufX2r=NULL;
  SendBufX0l_comp=NULL; RecvBufX0l_comp=NULL; 
  SendBufX0r_comp=NULL; RecvBufX0r_comp=NULL;
  SendBufX1l_comp=NULL; RecvBufX1l_comp=NULL; 
  SendBufX1r_comp=NULL; RecvBufX1r_comp=NULL;
  SendBufX2l_comp=NULL; RecvBufX2l_comp=NULL; 
  SendBufX2r_comp=NULL; RecvBufX2r_comp=NULL;

  // allocate internal arrays
//  data = (float **) malloc(Nspecies * sizeof(float *));
  data = new float *[Nspecies];
  int loclen = (Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r);
  for (int idat=0; idat<Nspecies; idat++)
    data[idat] = new float[loclen];

  // set flag denoting array ownership
  owndata = true;
}


// Vector Constructor using pre-defined data 
EnzoVector::EnzoVector(int N0, int N1, int N2, int G0l, int G0r, 
		       int G1l, int G1r, int G2l, int G2r, int Ns,
		       int NBx0L, int NBx0R, int NBx1L, int NBx1R,
		       int NBx2L, int NBx2R, float **userdata)
{
  // set internal length specifications
  Nx0 = N0;
  Nx1 = N1;
  Nx2 = N2;
  Ng0l = G0l;
  Ng0r = G0r;
  Ng1l = G1l;
  Ng1r = G1r;
  Ng2l = G2l;
  Ng2r = G2r;
  Nspecies = Ns;
  int Nlocal = Nx0*Nx1*Nx2*Nspecies;
#ifdef USE_MPI
  if ((NBx0L == MPI_PROC_NULL) && (NBx0R == MPI_PROC_NULL) &&
      (NBx1L == MPI_PROC_NULL) && (NBx1R == MPI_PROC_NULL) && 
      (NBx2L == MPI_PROC_NULL) && (NBx2R == MPI_PROC_NULL)) {
    Nglobal = Nlocal;
  }
  else {
    MPI_Arg one = 1;
    MPI_Datatype DataType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
    MPI_Allreduce(&Nlocal, &Nglobal, one, DataType, MPI_SUM, MPI_COMM_WORLD);  
  }
#else
  Nglobal = Nlocal;
#endif

  // set parallelism information on neighbors, other info to defaults
  Nbors[0][0] = NBx0L;  Nbors[0][1] = NBx0R;
  Nbors[1][0] = NBx1L;  Nbors[1][1] = NBx1R;
  Nbors[2][0] = NBx2L;  Nbors[2][1] = NBx2R;
  NborGhosts[0][0] = 0;  NborGhosts[0][1] = 0;
  NborGhosts[1][0] = 0;  NborGhosts[1][1] = 0;
  NborGhosts[2][0] = 0;  NborGhosts[2][1] = 0;
  x0Lsendsize = x0Rsendsize = x0Lrecvsize = x0Rrecvsize = 0;
  x1Lsendsize = x1Rsendsize = x1Lrecvsize = x1Rrecvsize = 0;
  x2Lsendsize = x2Rsendsize = x2Lrecvsize = x2Rrecvsize = 0;
  SendBufX0l=NULL; RecvBufX0l=NULL; SendBufX0r=NULL; RecvBufX0r=NULL;
  SendBufX1l=NULL; RecvBufX1l=NULL; SendBufX1r=NULL; RecvBufX1r=NULL;
  SendBufX2l=NULL; RecvBufX2l=NULL; SendBufX2r=NULL; RecvBufX2r=NULL;
  SendBufX0l_comp=NULL; RecvBufX0l_comp=NULL; 
  SendBufX0r_comp=NULL; RecvBufX0r_comp=NULL;
  SendBufX1l_comp=NULL; RecvBufX1l_comp=NULL; 
  SendBufX1r_comp=NULL; RecvBufX1r_comp=NULL;
  SendBufX2l_comp=NULL; RecvBufX2l_comp=NULL; 
  SendBufX2r_comp=NULL; RecvBufX2r_comp=NULL;

  // set internal arrays to point at arguments
//  data = (float **) malloc(Nspecies * sizeof(float *));
  data = new float *[Nspecies];
  for (int idat=0; idat<Nspecies; idat++)
    data[idat] = userdata[idat];

  // set flag denoting array ownership
  owndata = false;
}


// Vector Constructor for EnzoVector shell only
EnzoVector::EnzoVector(int N0, int N1, int N2, int G0l, int G0r, 
		       int G1l, int G1r, int G2l, int G2r, int Ns,
		       int NBx0L, int NBx0R, int NBx1L, int NBx1R,
		       int NBx2L, int NBx2R, int empty)
{
  // set internal length specifications
  Nx0 = N0;
  Nx1 = N1;
  Nx2 = N2;
  Ng0l = G0l;
  Ng0r = G0r;
  Ng1l = G1l;
  Ng1r = G1r;
  Ng2l = G2l;
  Ng2r = G2r;
  Nspecies = Ns;
  int Nlocal = Nx0*Nx1*Nx2*Nspecies;
#ifdef USE_MPI
  if ((NBx0L == MPI_PROC_NULL) && (NBx0R == MPI_PROC_NULL) &&
      (NBx1L == MPI_PROC_NULL) && (NBx1R == MPI_PROC_NULL) && 
      (NBx2L == MPI_PROC_NULL) && (NBx2R == MPI_PROC_NULL)) {
    Nglobal = Nlocal;
  }
  else {
    MPI_Arg one = 1;
    MPI_Datatype DataType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
    MPI_Allreduce(&Nlocal, &Nglobal, one, DataType, MPI_SUM, MPI_COMM_WORLD);  
  }
#else
  Nglobal = Nlocal;
#endif

  // set parallelism information on neighbors, other info to defaults
  Nbors[0][0] = NBx0L;  Nbors[0][1] = NBx0R;
  Nbors[1][0] = NBx1L;  Nbors[1][1] = NBx1R;
  Nbors[2][0] = NBx2L;  Nbors[2][1] = NBx2R;
  NborGhosts[0][0] = 0;  NborGhosts[0][1] = 0;
  NborGhosts[1][0] = 0;  NborGhosts[1][1] = 0;
  NborGhosts[2][0] = 0;  NborGhosts[2][1] = 0;
  x0Lsendsize = x0Rsendsize = x0Lrecvsize = x0Rrecvsize = 0;
  x1Lsendsize = x1Rsendsize = x1Lrecvsize = x1Rrecvsize = 0;
  x2Lsendsize = x2Rsendsize = x2Lrecvsize = x2Rrecvsize = 0;
  SendBufX0l=NULL; RecvBufX0l=NULL; SendBufX0r=NULL; RecvBufX0r=NULL;
  SendBufX1l=NULL; RecvBufX1l=NULL; SendBufX1r=NULL; RecvBufX1r=NULL;
  SendBufX2l=NULL; RecvBufX2l=NULL; SendBufX2r=NULL; RecvBufX2r=NULL;
  SendBufX0l_comp=NULL; RecvBufX0l_comp=NULL; 
  SendBufX0r_comp=NULL; RecvBufX0r_comp=NULL;
  SendBufX1l_comp=NULL; RecvBufX1l_comp=NULL; 
  SendBufX1r_comp=NULL; RecvBufX1r_comp=NULL;
  SendBufX2l_comp=NULL; RecvBufX2l_comp=NULL; 
  SendBufX2r_comp=NULL; RecvBufX2r_comp=NULL;

  // set internal arrays to point at arguments
//  data = (float **) malloc(Nspecies * sizeof(float *));
  data = new float *[Nspecies];
  for (int idat=0; idat<Nspecies; idat++)
    data[idat] = NULL;

  // set flag denoting array ownership
  owndata = false;
}


// Vector Destructor
EnzoVector::~EnzoVector()
{
  // check data ownership, and deallocate internal arrays if necessary
  if (owndata) {
    for (int idat=0; idat<Nspecies; idat++)
      delete[] data[idat];
  }
  delete[] data;
  delete[] SendBufX0l;
  delete[] SendBufX0l_comp;
  delete[] RecvBufX0l;
  delete[] RecvBufX0l_comp;
  delete[] SendBufX0r;
  delete[] SendBufX0r_comp;
  delete[] RecvBufX0r;
  delete[] RecvBufX0r_comp;
  delete[] SendBufX1l;
  delete[] SendBufX1l_comp;
  delete[] RecvBufX1l;
  delete[] RecvBufX1l_comp;
  delete[] SendBufX1r;
  delete[] SendBufX1r_comp;
  delete[] RecvBufX1r;
  delete[] RecvBufX1r_comp;
  delete[] SendBufX2l;
  delete[] SendBufX2l_comp;
  delete[] RecvBufX2l;
  delete[] RecvBufX2l_comp;
  delete[] SendBufX2r;
  delete[] SendBufX2r_comp;
  delete[] RecvBufX2r;
  delete[] RecvBufX2r_comp;
  MPI_Request *tmp;
  tmp = (MPI_Request *) id_recv_x0l;  delete tmp;
  tmp = (MPI_Request *) id_recv_x0l_comp;  delete tmp;
  tmp = (MPI_Request *) id_send_x0l;  delete tmp;
  tmp = (MPI_Request *) id_send_x0l_comp;  delete tmp;
  tmp = (MPI_Request *) id_recv_x0r;  delete tmp;
  tmp = (MPI_Request *) id_recv_x0r_comp;  delete tmp;
  tmp = (MPI_Request *) id_send_x0r;  delete tmp;
  tmp = (MPI_Request *) id_send_x0r_comp;  delete tmp;
  tmp = (MPI_Request *) id_recv_x1l;  delete tmp;
  tmp = (MPI_Request *) id_recv_x1l_comp;  delete tmp;
  tmp = (MPI_Request *) id_send_x1l;  delete tmp;
  tmp = (MPI_Request *) id_send_x1l_comp;  delete tmp;
  tmp = (MPI_Request *) id_recv_x1r;  delete tmp;
  tmp = (MPI_Request *) id_recv_x1r_comp;  delete tmp;
  tmp = (MPI_Request *) id_send_x1r;  delete tmp;
  tmp = (MPI_Request *) id_send_x1r_comp;  delete tmp;
  tmp = (MPI_Request *) id_recv_x2l;  delete tmp;
  tmp = (MPI_Request *) id_recv_x2l_comp;  delete tmp;
  tmp = (MPI_Request *) id_send_x2l;  delete tmp;
  tmp = (MPI_Request *) id_send_x2l_comp;  delete tmp;
  tmp = (MPI_Request *) id_recv_x2r;  delete tmp;
  tmp = (MPI_Request *) id_recv_x2r_comp;  delete tmp;
  tmp = (MPI_Request *) id_send_x2r;  delete tmp;
  tmp = (MPI_Request *) id_send_x2r_comp;  delete tmp;
}


//  Vector clone
EnzoVector* EnzoVector::clone() const
{
  return new EnzoVector(Nx0, Nx1, Nx2, Ng0l, Ng0r, Ng1l, 
			Ng1r, Ng2l, Ng2r, Nspecies, 
			Nbors[0][0], Nbors[0][1], 
			Nbors[1][0], Nbors[1][1], 
			Nbors[2][0], Nbors[2][1]);
}


//  Vector clone using pre-defined data
EnzoVector* EnzoVector::clone(float **userdata) const
{
  return new EnzoVector(Nx0, Nx1, Nx2, Ng0l, Ng0r, Ng1l, 
			Ng1r, Ng2l, Ng2r, Nspecies, 
			Nbors[0][0], Nbors[0][1], 
			Nbors[1][0], Nbors[1][1], 
			Nbors[2][0], Nbors[2][1], userdata);
}


//  Write vector species to a given file (no ghosts)
int EnzoVector::write(char *outfile, int species) const
{
  if ((species<0) || (species >= Nspecies)) {
    fprintf(stderr,"EnzoVector::write ERROR, illegal species %"ISYM"\n",
	    species);
    ENZO_FAIL("EnzoVector::write error");
  }

  // append processor number to file name
  char *filename = new char[strlen(outfile)+7];
  char *tmp_str  = new char[6];
  sprintf(tmp_str,"%"ISYM,MyProcessorNumber);
  strcpy(filename, outfile);
  strcat(strcat(filename,"."),tmp_str);
  FILE *filePtr;
  filePtr = fopen(filename, "w");
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int idx;
  for (int k=0; k<Nx2; k++) 
    for (int j=0; j<Nx1; j++)
      for (int i=0; i<Nx0; i++) {
	idx = ((k+Ng2l)*x1len + j+Ng1l)*x0len + i+Ng0l;
	fprintf(filePtr,"%"ISYM" %"ISYM" %"ISYM" %24.16g\n", 
		i+1, j+1, k+1, data[species][idx]);
      }
  fclose(filePtr);
  delete[] tmp_str;
  delete[] filename;
  return SUCCESS;
}


//  Write vector species to a given file (with ghosts)
int EnzoVector::writeall(char *outfile, int species) const
{
  if ((species<0) || (species >= Nspecies)) {
    fprintf(stderr,"EnzoVector::writeall ERROR, illegal species %"ISYM"\n",
	    species);
    ENZO_FAIL("EnzoVector::writeall error");
  }

  // append processor number to file name
  char *filename = new char[strlen(outfile)+7];
  char *tmp_str  = new char[6];
  sprintf(tmp_str,"%"ISYM,MyProcessorNumber);
  strcpy(filename, outfile);
  strcat(strcat(filename,"."),tmp_str);
  FILE *filePtr;
  filePtr = fopen(filename, "w");
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int x2len = Nx2 + Ng2l + Ng2r;
  int idx;
  for (int k=0; k<x2len; k++) 
    for (int j=0; j<x1len; j++)
      for (int i=0; i<x0len; i++) {
	idx = (k*x1len + j)*x0len + i;
	fprintf(filePtr,"%"ISYM" %"ISYM" %"ISYM" %24.16g\n", 
		i+1, j+1, k+1, data[species][idx]);
      }
  fclose(filePtr);
  delete[] tmp_str;
  delete[] filename;
  return SUCCESS;
}


//  Get data array for a given species 
float* EnzoVector::GetData(int species)
{
  if ((species >= 0) && (species < Nspecies))
    return &(data[species][0]);
  else {
    return NULL;
  }
}


//  Get data array for a given species 
int EnzoVector::SetData(int species, float *NewPtr)
{
  // check that species within range 
  if ((species < 0) || (species >= Nspecies)) {
    fprintf(stderr,"EnzoVector::SetData ERROR, illegal species %"ISYM"\n",
	    species);
    ENZO_FAIL("EnzoVector::SetData error");
  }

  // check that data array not owned by vector (memory leak)
  if (owndata) 
    ENZO_FAIL("EnzoVector::SetData ERROR, vector owns data array");

  // otherwise, set the data array to point at the new memory
  data[species] = NewPtr;
  return SUCCESS;
}


//   Access dimensional data
int EnzoVector::size(int *n0, int *n1, int *n2, int *ns, int *g0l, 
		     int *g0r, int *g1l, int *g1r, int *g2l, int *g2r)
{
  *n0 = Nx0;  *n1 = Nx1;  *n2 = Nx2; *ns = Nspecies;
  *g0l = Ng0l;  *g0r = Ng0r;
  *g1l = Ng1l;  *g1r = Ng1r;
  *g2l = Ng2l;  *g2r = Ng2r;
  return SUCCESS;
}


//  Vector Copy (assumes vectors have same size)
//  [operates on ghost zones as well as active data]
int EnzoVector::copy(EnzoVector *x)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector copy ERROR: vector sizes do not match");
  for (int idat=0; idat<Nspecies; idat++)
    for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
      data[idat][i] = x->data[idat][i];
  return SUCCESS;
}


//  Single component vector Copy (assumes vectors have same size)
//  [operates on ghost zones as well as active data]
int EnzoVector::copy_component(EnzoVector *x, int c)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector copy_component ERROR: vector sizes do not match");
  if ((c<0) || (c>=Nspecies)) {
    fprintf(stderr,"copy_component ERROR: illegal component %"ISYM"\n",c);
    ENZO_FAIL("EnzoVector copy_component ERROR");
  }
  for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
    data[c][i] = x->data[c][i];
  return SUCCESS;
}


//  Vector linear sum (assumes vectors have same size)
//  [operates on ghost zones as well as active data]
int EnzoVector::linearsum(float a, EnzoVector *x, float b, EnzoVector *y)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies) || 
      (y->Nx0 != Nx0) || (y->Nx1 != Nx1) || (y->Nx2 != Nx2) ||
      (y->Ng0l != Ng0l) || (y->Ng1l != Ng1l) || (y->Ng2l != Ng2l) ||
      (y->Ng0r != Ng0r) || (y->Ng1r != Ng1r) || (y->Ng2r != Ng2r) || 
      (y->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector linearsum ERROR: vector sizes do not match");
  for (int idat=0; idat<Nspecies; idat++)
    for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
      data[idat][i] = a*x->data[idat][i] + b*y->data[idat][i];
  return SUCCESS;
}


//  Vector axpy (assumes vectors have same size)
//  [operates on ghost zones as well as active data]
int EnzoVector::axpy(float a, EnzoVector *x)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector axpy ERROR: vector sizes do not match");
  for (int idat=0; idat<Nspecies; idat++)
    for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
      data[idat][i] += a*x->data[idat][i];
  return SUCCESS;
}


//  Single component vector axpy (assumes vectors have same size)
//  [operates on ghost zones as well as active data]
int EnzoVector::axpy_component(float a, EnzoVector *x, int c)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector axpy_component ERROR: vector sizes do not match");
  if ((c<0) || (c>=Nspecies)) {
    fprintf(stderr,"axpy_component ERROR: illegal component %"ISYM"\n",c);
    ENZO_FAIL("EnzoVector axpy_component ERROR");
  }
  for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)
    data[c][i] += a*x->data[c][i];
  return SUCCESS;
}


//  Vector component scale operation, y *= a
//  [operates on ghost zones as well as active data]
int EnzoVector::scale_component(int idat, float a)
{
  if ((idat < 0) || (idat >= Nspecies)) {
    fprintf(stderr,"scale_component error: illegal var = %"ISYM"\n",idat);
    return -1.0;
  }
  for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
    data[idat][i] *= a;
  return SUCCESS;
}


//  [operates on ghost zones as well as active data]
int EnzoVector::log_component(int idat)
{
  if ((idat < 0) || (idat >= Nspecies)) {
    fprintf(stderr,"log_component error: illegal var = %"ISYM"\n",idat);
    return -1.0;
  }
  for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
    data[idat][i] = (data[idat][i] == 0.0) ? -1.2e4 : log(data[idat][i]);
  return SUCCESS;
}


//  [operates on ghost zones as well as active data]
int EnzoVector::exp_component(int idat)
{
  if ((idat < 0) || (idat >= Nspecies)) {
    fprintf(stderr,"exp_component error: illegal var = %"ISYM"\n",idat);
    return -1.0;
  }
  for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
    data[idat][i] = exp(data[idat][i]);
  return SUCCESS;
}


//  Vector scale operation, y *= a
//  [operates on ghost zones as well as active data]
int EnzoVector::scale(float a)
{
  for (int idat=0; idat<Nspecies; idat++)
    for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
      data[idat][i] *= a;
  return SUCCESS;
}


//  Vector constant operation, y(i) = a
//  [operates on ghost zones as well as active data]
int EnzoVector::constant(float a)
{
  for (int idat=0; idat<Nspecies; idat++)
    for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
      data[idat][i] = a;
  return SUCCESS;
}


//  Vector add const operation, y(i) = y(i) + a
//  [operates on ghost zones as well as active data]
int EnzoVector::addconst(float a)
{
  for (int idat=0; idat<Nspecies; idat++)
    for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
      data[idat][i] += a;
  return SUCCESS;
}


//  Vector component add const operation, y(i) = y(i) + a
//  [operates on ghost zones as well as active data]
int EnzoVector::addconst_component(int idat, float a)
{
  if ((idat < 0) || (idat >= Nspecies)) {
    fprintf(stderr,"addconst_component error: illegal var = %"ISYM"\n",idat);
    return -1.0;
  }
  for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
    data[idat][i] += a;
  return SUCCESS;
}


//  Vector absolute value operation, this(i) = |x(i)|
//  [operates on ghost zones as well as active data]
int EnzoVector::abs(EnzoVector *x)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies))
    ENZO_FAIL("EnzoVector abs ERROR: vector sizes do not match");
  for (int idat=0; idat<Nspecies; idat++)
    for (int i=0; i<((Nx0+Ng0l+Ng0r)*(Nx1+Ng1l+Ng1r)*(Nx2+Ng2l+Ng2r)); i++)  
      data[idat][i] = fabs(x->data[idat][i]);
  return SUCCESS;
}


//  Vector product operation, this(i) = x(i)*y(i)
//  [operates on active data only]
int EnzoVector::product(EnzoVector *x, EnzoVector *y)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies) || 
      (y->Nx0 != Nx0) || (y->Nx1 != Nx1) || (y->Nx2 != Nx2) ||
      (y->Ng0l != Ng0l) || (y->Ng1l != Ng1l) || (y->Ng2l != Ng2l) ||
      (y->Ng0r != Ng0r) || (y->Ng1r != Ng1r) || (y->Ng2r != Ng2r) || 
      (y->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector product ERROR: vector sizes do not match");

  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  idx = (k*x1len + j)*x0len + i;
	  data[idat][idx] = (x->data[idat][idx])*(y->data[idat][idx]);
	}
  return SUCCESS;
}


//  Vector quotient operation, this(i) = x(i)/y(i)
//  [operates on active data only; assumes y(i)!=0]
int EnzoVector::quotient(EnzoVector *x, EnzoVector *y)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies) || 
      (y->Nx0 != Nx0) || (y->Nx1 != Nx1) || (y->Nx2 != Nx2) ||
      (y->Ng0l != Ng0l) || (y->Ng1l != Ng1l) || (y->Ng2l != Ng2l) ||
      (y->Ng0r != Ng0r) || (y->Ng1r != Ng1r) || (y->Ng2r != Ng2r) || 
      (y->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector quotient ERROR: vector sizes do not match");

  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  idx = (k*x1len + j)*x0len + i;
	  data[idat][idx] = (x->data[idat][idx])/(y->data[idat][idx]);
	}
  return SUCCESS;
}



//  Vector minquotient operation, min(this(i)/y(i)), where y(i)!=0
//  [operates on active data only]
float EnzoVector::minquotient(EnzoVector *y)
{
  if ((y->Nx0 != Nx0) || (y->Nx1 != Nx1) || (y->Nx2 != Nx2) ||
      (y->Ng0l != Ng0l) || (y->Ng1l != Ng1l) || (y->Ng2l != Ng2l) ||
      (y->Ng0r != Ng0r) || (y->Ng1r != Ng1r) || (y->Ng2r != Ng2r) || 
      (y->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector minquotient ERROR: vector sizes do not match");

  float minquot = huge_number;
  float quotient, numer, denom;
  bool noneyet=true;
  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  idx = (k*x1len + j)*x0len + i;
	  numer = data[idat][idx];
	  denom = y->data[idat][idx];
	  if (denom != 0.0) {
	    quotient = numer/denom;
	    if (noneyet)  minquot = quotient;
	    else  minquot = (quotient<minquot) ? quotient : minquot;
	  }
	}
  return minquot;
}


//  Vector inverse operation, this(i) = 1.0/x(i)
//  [operates on active data only; assumes x(i)!=0]
int EnzoVector::inverse(EnzoVector *x)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector inverse ERROR: vector sizes do not match");

  float value;
  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  idx = (k*x1len + j)*x0len + i;
	  value = x->data[idat][idx];
	  data[idat][idx] = (value != 0.0) ? 1.0/value : 1.0;
	}
  return SUCCESS;
}


//  Vector constraint checking operation
//  this(i) = 0.0 where constraints true, 1.0 where false
//  Test:  if c[i] =  2.0, then x[i] must be >  0.0
//         if c[i] =  1.0, then x[i] must be >= 0.0
//         if c[i] = -1.0, then x[i] must be <= 0.0
//         if c[i] = -2.0, then x[i] must be <  0.0
//  if all constraints satisfied, returns true
//  [operates on active data only]
bool EnzoVector::constraintcheck(EnzoVector *c, EnzoVector *x)
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies) ||
      (c->Nx0 != Nx0) || (c->Nx1 != Nx1) || (c->Nx2 != Nx2) ||
      (c->Ng0l != Ng0l) || (c->Ng1l != Ng1l) || (c->Ng2l != Ng2l) ||
      (c->Ng0r != Ng0r) || (c->Ng1r != Ng1r) || (c->Ng2r != Ng2r) || 
      (c->Nspecies != Nspecies)) 
    ENZO_FAIL("EnzoVector constraintcheck ERROR: vector sizes do not match");

  bool test=TRUE;
  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  idx = (k*x1len + j)*x0len + i;
	  // check for strict positivity/negativity
	  if (fabs(c->data[idat][idx]) > 1.5) {
	    if ((c->data[idat][idx])*(x->data[idat][idx]) > 0.0) 
	      data[idat][idx] = 0.0;
	    else {
	      data[idat][idx] = 1.0;
	      test = false;
	    }
	  }
	  // check for loose positivity/negativity
	  else if (fabs(c->data[idat][idx]) > 0.5) {
	    if ((c->data[idat][idx])*(x->data[idat][idx]) >= 0.0) 
	      data[idat][idx] = 0.0;
	    else {
	      data[idat][idx] = 1.0;
	      test = false;
	    }
	  }
	}
  return test;
}


//  Vector dot product
float EnzoVector::dot(EnzoVector *x) const
{
  if ((x->Nx0 != Nx0) || (x->Nx1 != Nx1) || (x->Nx2 != Nx2) ||
      (x->Ng0l != Ng0l) || (x->Ng1l != Ng1l) || (x->Ng2l != Ng2l) ||
      (x->Ng0r != Ng0r) || (x->Ng1r != Ng1r) || (x->Ng2r != Ng2r) || 
      (x->Nspecies != Nspecies)) {
    fprintf(stderr,"EnzoVector dot ERROR: vector sizes do not match\n");
    return 0.0;
  }
  float sum=0.0, gsum;
  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  idx = (k*x1len + j)*x0len + i;
	  sum += data[idat][idx]*x->data[idat][idx];
	}

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gsum = sum;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&sum, &gsum, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  }
#else
  gsum = sum;
#endif

  return(gsum);
}


//  Vector rmsnorm
float EnzoVector::rmsnorm() const
{
  float sum=0.0, gsum;
  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  idx = (k*x1len + j)*x0len + i;
	  sum += data[idat][idx]*data[idat][idx];
	}

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gsum = sum;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&sum, &gsum, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  }
#else
  gsum = sum;
#endif

  return( sqrt(gsum/Nglobal) );
}


//   Component RMS norm,  sqrt(dot(this,this)/Nglobal)
float EnzoVector::rmsnorm_component(int idat) const
{
  if ((idat < 0) || (idat >= Nspecies)) {
    fprintf(stderr,"rmsnorm_component error: illegal var = %"ISYM"\n",idat);
    return -1.0;
  }
  float sum=0.0, gsum;
  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int k=Ng2l; k<Nx2+Ng2l; k++) 
    for (int j=Ng1l; j<Nx1+Ng1l; j++)
      for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	idx = (k*x1len + j)*x0len + i;
	sum += data[idat][idx]*data[idat][idx];
      }

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gsum = sum;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&sum, &gsum, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  }
#else
  gsum = sum;
#endif

  return( sqrt(gsum/Nglobal*Nspecies) );
}


//  Vector weighted RMS norm
float EnzoVector::wrmsnorm(EnzoVector *w) const
{
  float wl2 = this->wl2norm(w);
  return( sqrt(wl2*wl2/Nglobal) );
}


//  Vector weighted L-2 norm
float EnzoVector::wl2norm(EnzoVector *w) const
{
  if ((w->Nx0 != Nx0) || (w->Nx1 != Nx1) || (w->Nx2 != Nx2) ||
      (w->Ng0l != Ng0l) || (w->Ng1l != Ng1l) || (w->Ng2l != Ng2l) ||
      (w->Ng0r != Ng0r) || (w->Ng1r != Ng1r) || (w->Ng2r != Ng2r) || 
      (w->Nspecies != Nspecies)) {
    fprintf(stderr,"EnzoVector wl2norm ERROR: vector sizes do not match\n");
    return 0.0;
  }

  float sum=0.0, gsum, prod;
  int idx;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  idx = (k*x1len + j)*x0len + i;
	  prod = (data[idat][idx])*(w->data[idat][idx]);
	  sum += prod*prod;
	}

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gsum = sum;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&sum, &gsum, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  }
#else
  gsum = sum;
#endif

  return( sqrt(gsum) );
}


//  Vector L-1 norm
float EnzoVector::l1norm() const
{
  float sum=0.0, gsum;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++)
	  sum += fabs(data[idat][(k*x1len + j)*x0len + i]);

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gsum = sum;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&sum, &gsum, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  }
#else
  gsum = sum;
#endif

  return(gsum);
}


//  Vector infnorm
float EnzoVector::infnorm() const
{
  float max=0.0, gmax;
  float tmp;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  tmp = fabs(data[idat][(k*x1len + j)*x0len + i]);
	  max = (max > tmp) ? max : tmp;
	}

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gmax = max;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&max, &gmax, one, DataType, MPI_MAX, MPI_COMM_WORLD);
  }
#else
  gmax = max;
#endif

  return(gmax);
}


//  Component infnorm
float EnzoVector::infnorm_component(int var) const
{
  if ((var < 0) || (var >= Nspecies)) {
    fprintf(stderr,"infnorm_component error: illegal var = %"ISYM"\n",var);
    return -1.0;
  }

  float max=0.0, gmax;
  float tmp;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int k=Ng2l; k<Nx2+Ng2l; k++) 
    for (int j=Ng1l; j<Nx1+Ng1l; j++)
      for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	tmp = fabs(data[var][(k*x1len + j)*x0len + i]);
	max = (max > tmp) ? max : tmp;
      }
  
#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gmax = max;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&max, &gmax, one, DataType, MPI_MAX, MPI_COMM_WORLD);
  }
#else
  gmax = max;
#endif

  return(gmax);
}


//   Relative pointwise difference,  max(abs(this-that)/abs(this))
float EnzoVector::relative_difference(float *x, int var) const
{
  if ((var < 0) || (var >= Nspecies)) {
    fprintf(stderr,"relative_difference error: illegal var = %"ISYM"\n",var);
    return -1.0;
  }

  float max=0.0, gmax;
  float tmp, diff, val;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int k=Ng2l; k<Nx2+Ng2l; k++) 
    for (int j=Ng1l; j<Nx1+Ng1l; j++)
      for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	diff = fabs(data[var][(k*x1len + j)*x0len + i] - 
		    x[(k*x1len + j)*x0len + i]);
	val  = fabs(data[var][(k*x1len + j)*x0len + i]);
	tmp = diff/val;
	max = (max > tmp) ? max : tmp;
      }

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gmax = max;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&max, &gmax, one, DataType, MPI_MAX, MPI_COMM_WORLD);
  }
#else
  gmax = max;
#endif

  return(gmax);
}


//   Relative volumetric difference, rmsnorm((abs(this-that)/abs(this)))
float EnzoVector::relative_vol_difference(float *x, int var) const
{
  if ((var < 0) || (var >= Nspecies)) {
    fprintf(stderr,"relative_vol_difference error: illegal var = %"ISYM"\n",var);
    return -1.0;
  }

  float sum=0.0, gsum;
  float tmp, diff, val;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  for (int k=Ng2l; k<Nx2+Ng2l; k++) 
    for (int j=Ng1l; j<Nx1+Ng1l; j++)
      for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	diff = fabs(data[var][(k*x1len + j)*x0len + i] - 
		    x[(k*x1len + j)*x0len + i]);
	val  = fabs(data[var][(k*x1len + j)*x0len + i]);
	tmp = diff/val;
	sum += tmp*tmp;
      }

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gsum = sum;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&sum, &gsum, one, DataType, MPI_SUM, MPI_COMM_WORLD);
  }
  gsum = gsum;
#else
  gsum = sum;
#endif

  return(sqrt(gsum/Nglobal*Nspecies));
}


//  Vector minimum value
float EnzoVector::minval() const
{
  float min, gmin;
  float tmp;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  min = data[0][(Ng2l*x1len + Ng1l)*x0len + Ng0l];
  for (int idat=0; idat<Nspecies; idat++)
    for (int k=Ng2l; k<Nx2+Ng2l; k++) 
      for (int j=Ng1l; j<Nx1+Ng1l; j++)
	for (int i=Ng0l; i<Nx0+Ng0l; i++) {
	  tmp = data[idat][(k*x1len + j)*x0len + i];
	  min = (min < tmp) ? min : tmp;
	}

#ifdef USE_MPI
  if (Nglobal == Nx0*Nx1*Nx2*Nspecies) 
    gmin = min;
  else {
    MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
    MPI_Arg one = 1;
    MPI_Allreduce(&min, &gmin, one, DataType, MPI_MIN, MPI_COMM_WORLD);
  }
#else
  gmin = min;
#endif

  return(gmin);
}



// end of file
