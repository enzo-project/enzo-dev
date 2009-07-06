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
/  This routine runs a number of checks on the EnzoVector class.
/
/  written by: Daniel R. Reynolds
/  date:       July, 2006
/
************************************************************************/
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "EnzoVector.h"

int EnzoVector::test()
{
  
  fprintf(stdout,"p%"ISYM": Entering EnzoVector::test routine\n",MyProcessorNumber);

  // output vector size diagnostics
  fprintf(stdout,"  p%"ISYM": active dims = (%"ISYM",%"ISYM",%"ISYM")\n",
	  MyProcessorNumber, Nx0, Nx1, Nx2);
  fprintf(stdout,"  p%"ISYM": ghosts = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",
	  MyProcessorNumber, Ng0l, Ng0r, Ng1l, Ng1r, Ng2l, Ng2r);
  fprintf(stdout,"  p%"ISYM": species = %"ISYM",  nglobal = %"ISYM"\n", MyProcessorNumber, Nspecies, Nglobal);
  fprintf(stdout,"  p%"ISYM": neighbors = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n", MyProcessorNumber, Nbors[0][0], Nbors[0][1], Nbors[1][0], Nbors[1][1], Nbors[2][0], Nbors[2][1]);
  

  // make a few copies for testing
  EnzoVector *x = this->clone();
  EnzoVector *y = this->clone();
  EnzoVector *z = this->clone();


  // initialize each vector to a constant
  x->constant(3.0);
  y->constant(2.0);
  z->constant(1.0);

  // test the vector constant, write routines
  fprintf(stdout,"  p%"ISYM": writing out the constant 3.0 to test1_vec\n", MyProcessorNumber);
  x->write("test1_vec",0);

  // test the norms
  fprintf(stdout,"  p%"ISYM": x->dot(y) = %"FSYM" (%"FSYM")\n", MyProcessorNumber, x->dot(y), 6.0*Nglobal); fflush(NULL);
  fprintf(stdout,"  p%"ISYM": x->rmsnorm() = %"FSYM" (%"FSYM")\n", MyProcessorNumber, x->rmsnorm(), 3.0); fflush(NULL);
  fprintf(stdout,"  p%"ISYM": x->wrmsnorm(y) = %"FSYM" (%"FSYM")\n", MyProcessorNumber, x->wrmsnorm(y), 6.0); fflush(NULL);
  fprintf(stdout,"  p%"ISYM": x->wl2norm(y) = %"FSYM" (%"FSYM")\n", MyProcessorNumber, x->wl2norm(y), sqrt(36.0*Nglobal)); fflush(NULL);
  fprintf(stdout,"  p%"ISYM": x->l1norm() = %"FSYM" (%"FSYM")\n", MyProcessorNumber, x->l1norm(), 3.0*Nglobal); fflush(NULL);
  fprintf(stdout,"  p%"ISYM": x->infnorm() = %"FSYM" (%"FSYM")\n", MyProcessorNumber, x->infnorm(), 3.0); fflush(NULL);
  fprintf(stdout,"  p%"ISYM": x->minval() = %"FSYM" (%"FSYM")\n", MyProcessorNumber, x->minval(), 3.0); fflush(NULL);
  
  // test the vector copy
  z->copy(x);
  fprintf(stdout,"  p%"ISYM": z->copy(x), printing to test2_vec\n", MyProcessorNumber); fflush(NULL);
  z->write("test2_vec",0);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // test the axpy operation
  z->axpy(-1.0,x);
  fprintf(stdout,"  p%"ISYM": axpy test error = %"FSYM"\n", MyProcessorNumber, z->infnorm()); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // test the vector absolute value
  z->constant(-1.0);
  y->abs(z);
  x->constant(1.0);
  x->axpy(-1.0,y);
  fprintf(stdout,"  p%"ISYM": abs test error = %"FSYM"\n", MyProcessorNumber, x->infnorm()); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // test the vector scale operation
  x->constant(3.0);
  y->constant(1.0);
  y->scale(-3.0);
  y->axpy(1.0,x);
  fprintf(stdout,"  p%"ISYM": scale test error = %"FSYM"\n", MyProcessorNumber, y->infnorm()); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // test the vector linear sum operation
  x->constant(3.0);
  y->constant(2.0);
  z->linearsum(2.0,x,-3.0,y);
  fprintf(stdout,"  p%"ISYM": linearsum test error = %"FSYM"\n", MyProcessorNumber, z->infnorm()); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // test the addconst operation
  y->addconst(-2.0);
  fprintf(stdout,"  p%"ISYM": addconst test error = %"FSYM"\n", MyProcessorNumber, y->infnorm()); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  // test the product operation
  x->constant(3.0);
  y->constant(2.0);
  z->product(x,y);
  x->scale(2.0);
  x->axpy(-1.0,z);
  fprintf(stdout,"  p%"ISYM": product test error = %"FSYM"\n", MyProcessorNumber, x->infnorm()); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  // test the quotient operation
  x->quotient(z,y);
  y->constant(3.0);
  y->axpy(-1.0,x);
  fprintf(stdout,"  p%"ISYM": quotient test error = %"FSYM"\n", MyProcessorNumber, y->infnorm()); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  // test the minquotient operation
  x->constant(3.0);
  y->constant(2.0);
  fprintf(stdout,"  p%"ISYM": minquotient test error = %"FSYM"\n", MyProcessorNumber, x->minquotient(y)-1.5); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // test the inverse operation
  x->constant(4.0);
  y->inverse(x);
  z->constant(0.25);
  z->axpy(-1.0,y);
  fprintf(stdout,"  p%"ISYM": inverse test error = %"FSYM"\n", MyProcessorNumber, z->infnorm()); fflush(NULL);
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  // test the communication routine
  float *xdata = x->GetData(0);
  int i, j, k, idx;
  for (k=0; k<Nx2+Ng2l+Ng2r; k++) 
    for (j=0; j<Nx1+Ng1l+Ng1r; j++) 
      for (i=0; i<Nx0+Ng0l+Ng0r; i++) {
	idx = (k*(Nx1+Ng1l+Ng1r) + j)*(Nx0+Ng0l+Ng0r) + i;
	xdata[idx] = 10000*k + 100*j + i;
      }
  fprintf(stdout,"  p%"ISYM": writing initial x to x_init\n", MyProcessorNumber);
  x->writeall("x_init",0);
  x->exchange();
  fprintf(stdout,"  p%"ISYM": writing communicated x to x_comm\n", MyProcessorNumber);
  x->writeall("x_comm",0);

  fprintf(stdout,"\n"); fflush(NULL);

  return SUCCESS;
}

/******************************************************************/
