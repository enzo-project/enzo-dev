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
/  Gray Flux-Limited Diffusion Implicit Problem linear Newton system 
/  solution function
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  August 13, 2007, by John Hayes; implemented function calls
/              for 2D and 1D versions of MatrixEntries and SetNewtonBCs.
/
/  PURPOSE: Solves the linear Newton system J(u)*s = b.  For the Gray 
/           FLD Problem (one radiation group) without advection, the 
/           problem may be solved via the following Schur-complement 
/           formulation.
/
/           The Jacobian may be written in the operator form: 
/                   [ L_ee  L_en    L_eE        ]   
/               J = [ L_ne  L_nn    L_nE        ] = [ M U ]
/                   [ L_Ee  L_En  (L_EE + D_EE) ]   [ L D ]
/           where for the fluid energy e, chemical species n, and 
/           radiation group r,
/               L_ee = local Jacobian of e wrt e
/               L_en = local Jacobian of e wrt n_j, j=1:Nchem
/               L_eE = local Jacobian of e wrt E
/               L_ne = local Jacobian of n_i wrt e, i=1:Nchem
/               L_nn = local Jacobian of n_i wrt n_j, i,j=1:Nchem
/               L_nE = local Jacobian of n_i wrt E, i=1:Nchem
/               L_Ee = local Jacobian of E wrt e
/               L_En = local Jacobian of E wrt n_j, j=1:Nchem
/               L_EE = local Jacobian of E wrt E
/               D_EE = diffusive Jacobian of E wrt E
/                 M = [ L_ee  L_en ]
/                     [ L_ne  L_nn ]
/                 U = [ L_eE L_nE ]^T
/                 L = [ L_Ee L_En ]
/                 D = (L_EE + D_EE).
/            The Schur complement formulation provides that
/                 Ji = [ I -Mi*U ][ Mi 0  ][   I   0 ]
/                      [ 0    I  ][ 0  Pi ][ -L*Mi I ]
/            where we use 'i' to denote the inverse, e.g. Ai = A^{-1}, 
/            and where the Schur complement is formed as P = D-L*Mi*U.  
/            Therefore, the solve J*s = b, where 
/            s = (s_e s_n s_E)^T = (s_m s_E)^T  (other vectors similarly 
/            condense e and n into m) may be broken down into the stages:
/                 (1) Solve for c_m:  M*c_m = b_m             (local)
/                 (2) Solve for y_m:  M*y_m = U               (local)
/                 (3) Construct:  y_E = L_EE - L*y_m          (local)
/                     note: as this matrix is diagonal, store in y_E 
/                 (4) Update: b_E = b_E - L*c_m               (local)
/                 (5) Construct: P = D_EE + I*y_E             (local)
/                 (6) Solve for s_E:  P*s_E = b_E             (nonlocal)
/                 (7) Set s_m = c_m - y_m*s_E.                (local)
/            We note that all of the steps are completely local except 
/            for the step (6), which requires one scalar-valued diffusive 
/            solve, for which we use the HYPRE library's PFMG solver 
/            (extensible to vector-valued diffusive systems).
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"


int gFLDProblem::lsolve(EnzoVector *s, EnzoVector *b, 
			EnzoVector *u, float delta)
{
//   if (debug)  printf("Entering gFLDProblem::lsolve routine\n");

  // check that the gFLDProblem has been set up
  if (!prepared) 
    ENZO_FAIL("lsolve error: gFLDProblem not yet prepared");

  // in case MPI is not included
#ifndef MPI_INT
  int MPI_COMM_WORLD = 0;
#endif
  
  // have b communicate neighbor information and enforce BCs
//   if (b->exchange() == FAIL) 
//      ENZO_FAIL("lsolve error: vector exchange failure");
  if (this->EnforceBoundary(b,1) == FAIL) 
    ENZO_FAIL("lsolve error: EnforceBoundary failure");

  // check that b matches local vector size (including ghosts, etc)
  int ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  int vsz[4];
  b->size(&vsz[0], &vsz[1], &vsz[2], &vsz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (vsz[0] != LocDims[0]) 
    ENZO_FAIL("lsolve error: x0 vector dims do not match");
  if (vsz[1] != LocDims[1]) 
    ENZO_FAIL("lsolve error: x1 vector dims do not match");
  if (vsz[2] != LocDims[2]) 
    ENZO_FAIL("lsolve error: x2 vector dims do not match");
  if (vsz[3] != (2+Nchem)) 
    ENZO_FAIL("lsolve error: nspecies do not match");
  if ((vsz[0]+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("lsolve error: x0 vector sizes do not match");
  if ((vsz[1]+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("lsolve error: x1 vector sizes do not match");
  if ((vsz[2]+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("lsolve error: x2 vector sizes do not match");

//   if (debug)
//     printf("gFLDProblem::lsolve -- creating Schur comp. temp vector\n");

  // create temporary vector yvec for Schur complement correction
  EnzoVector *yvec = tmp1;
  yvec->constant(0.0);

//   if (debug)  printf("gFLDProblem::lsolve -- performing steps 1-2\n");

  //////////////////////////////////////////////////////////////
  // steps (1) and (2): local solves 
  //          c_m = Mi*b_m   (store c_m in b_m) 
  //          y_m = Mi*U
  int ix, iy, iz, idx, irow, icol, size=Nchem+1, two=2;
  float *Lblock, *tmp;
  float *M = new float[size*size];
  float *bvec = new float[2*size];
  float *xvec = new float[2*size];
  float D;
  for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
    for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
      for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;

	// set up local multiple RHS matrix system 
	//     M*[c_m, y_m] = [b_m, U]
	for (irow=0; irow<size; irow++) {
	  for (icol=0; icol<size; icol++) {

	    // set matrix element with correct Jacobian component
	    Lblock = L[irow+1]->GetData(icol+1);
	    M[icol*size+irow] = Lblock[idx];
	  }

	  // set up rhs vector bvec (1st column is b_m, 2nd is U)
	  //    1st column contains b_m
	  tmp = b->GetData(irow+1);
	  bvec[irow] = tmp[idx];
	  //    2nd column contains U
	  tmp = (L[irow+1])->GetData(0);
	  bvec[size+irow] = tmp[idx];
	}


	// solve dense local systems
	//   without chemistry, this is a pair of scalar problems
	if (Nchem == 0) {
	  xvec[0] = bvec[0]/M[0];
	  xvec[1] = bvec[1]/M[0];
	}
	else if (Nchem == 1) {
	  D = M[0]*M[3] - M[1]*M[2];
	  xvec[0] = (M[3]*bvec[0] - M[2]*bvec[1])/D;
	  xvec[1] = (M[0]*bvec[1] - M[1]*bvec[0])/D;
	  xvec[2] = (M[3]*bvec[2] - M[2]*bvec[3])/D;
	  xvec[3] = (M[0]*bvec[3] - M[1]*bvec[2])/D;
	}
	else {
	  if (this->BlockSolve(M, xvec, bvec, &size, &two) != SUCCESS) 
	    ENZO_FAIL("lsolve: Error in BlockSolve routine");
	}

	// extract solution components to appropriate locations
	for (irow=0; irow<size; irow++) {
	  //    put c_m back in b_m
	  tmp = b->GetData(irow+1);
	  tmp[idx]= xvec[irow];

	  //    put Mi*U into y_m vector (all of yvec except species 0)
	  tmp = yvec->GetData(irow+1);
	  tmp[idx] = xvec[size+irow];
	}
      }
    }
  }
  delete[] M;
  delete[] xvec;
  delete[] bvec;
//   if (debug)  printf("gFLDProblem::lsolve -- performing steps 3-4\n");

  //////////////////////////////////////////////////////////////
  // steps (3) and (4): Construct local update for P and adjust rhs b_E
  //     (3) construct:  y_E = L_EE - L*y_m
  //     (4) update:     b_E = b_E - L*c_m  (note: c_m stored in b_m)
  float *y_E = yvec->GetData(0);
  float *b_E = b->GetData(0);
  float *Ldiag = (L[0])->GetData(0);
  float *b_m, *y_m;
  for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
    for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
      for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	
 	// set y_E with L_EE at first
	y_E[idx] = Ldiag[idx];
	
	// iterate over other fluid energy, chemistry 
	// (radiation in irow 0)
	for (irow=1; irow<(2+Nchem); irow++) {

	  // L_E* is contained in L[0]
	  Lblock = (L[0])->GetData(irow);

	  // update y_E component
	  //   Mi*U is contained in y_m (the rest of yvec)
	  y_m = yvec->GetData(irow);
	  y_E[idx] -= Lblock[idx]*y_m[idx];

	  // update rhs vector b_E
	  //   store c_m in b_m
	  b_m = b->GetData(irow);
	  b_E[idx] -= Lblock[idx]*b_m[idx];
	}
      }
    }
  }


//   if (debug)  printf("gFLDProblem::lsolve -- performing step 5\n");

  //////////////////////////////////////////////////////////////
  // step (5): Construct (locally) the Schur complement matrix 
  //    P = D-L*Mi*U.  In the code, this corresponds to the 
  //    matrix  P = D_EE + I*y_E;  we then scale the system 
  //    via  P = diag(P)^{-1}*P,  b_E = diag(P)^{-1}*b_E
  float *s_E = s->GetData(0);

  //       communicate yvec to spread local corrections to neighbors
  if (yvec->exchange_component(0) == FAIL) 
    ENZO_FAIL("lsolve error: vector exchange_component error");
    
#ifdef USE_HYPRE

  //       set matrix values over grid
  float *u_E = u->GetData(0);
  float *u0_E = U0->GetData(0);
  int Nx = (SolvIndices[0][1]-SolvIndices[0][0]+1);
  int Ny = (SolvIndices[1][1]-SolvIndices[1][0]+1);
  int Nz = (SolvIndices[2][1]-SolvIndices[2][0]+1);
  for (ix=0; ix<stSize*Nx*Ny*Nz; ix++)  Ptmpvec[ix]=0.0;
  //       NOTE: Temp and OpacityE are still valid from nlresid
  //
  if (this->MatrixEntries(Ptmpvec, u_E, u0_E, Temp, OpacityE, y_E) != SUCCESS) 
    ENZO_FAIL("lsolve: Error in MatrixEntries routine");
  if (this->SetNewtonBCs(Ptmpvec, b_E) != SUCCESS) 
    ENZO_FAIL("lsolve: Error in SetNewtonBCs routine");
//   if (debug)  printf("lsolve: calling HYPRE_StructMatrixSetBoxValues\n");
  Eint32 zed=0;
  Eint32 one=1;
  Eint32 entries[7] = {0, 1, 2, 3, 4, 5, 6};
  Eint32 ilower[3] = {SolvIndices[0][0],SolvIndices[1][0],SolvIndices[2][0]};
  Eint32 iupper[3] = {SolvIndices[0][1],SolvIndices[1][1],SolvIndices[2][1]};
  HYPRE_StructMatrixSetBoxValues(P, ilower, iupper, stSize, entries, Ptmpvec); 

  //       assemble matrix
//   if (debug)  printf("lsolve: calling HYPRE_StructMatrixAssemble\n");
  HYPRE_StructMatrixAssemble(P);

  // set symmetry of matrix
  //  HYPRE_StructMatrixSetSymmetric(P, one);

//   if (debug)  printf("gFLDProblem::lsolve -- performing step 6\n");

  //////////////////////////////////////////////////////////////
  // step (6): Solve the (nonlocal) Schur complement system,
  // i.e. solve for s_E:  P*s_E = b_E

  //       re-scale delta to relative residual and not actual
  delta /= b->rmsnorm();
  delta = min(delta, 1.0e-6);

  //       insert rhs, sol vectors into HYPRE vectors x and b
  ilower[0] = SolvIndices[0][0];
  iupper[0] = SolvIndices[0][1];
  int xBuff, yBuff, zBuff;
  xBuff = ghXl-SolvOff[0];
  yBuff = (ghYl-SolvOff[1])-SolvIndices[1][0];
  zBuff = (ghZl-SolvOff[2])-SolvIndices[2][0];
  int Zbl, Ybl;
//   if (debug)  printf("lsolve: calling HYPRE_StructVectorSetBoxValues\n");
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper[1] = iy;
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	HYPREbuff[ix] = b_E[Zbl+Ybl+xBuff+ix];
      HYPRE_StructVectorSetBoxValues(rhsvec, ilower, iupper, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	HYPREbuff[ix] = 0.0;
      HYPRE_StructVectorSetBoxValues(solvec, ilower, iupper, HYPREbuff);
    }
  }

  //       assemble vectors
//   if (debug)  printf("lsolve: calling HYPRE_StructVectorAssemble\n");
  HYPRE_StructVectorAssemble(solvec);
  HYPRE_StructVectorAssemble(rhsvec);

  //       set up the solver [PCG] and preconditioner [PFMG]
  //          create the solver & preconditioner
  HYPRE_StructSolver solver;
  HYPRE_StructSolver preconditioner;
//   if (debug)  printf("lsolve: calling HYPRE_StructPCGCreate\n");
  HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
//   if (debug)  printf("lsolve: calling HYPRE_StructPFMGCreate\n");
  HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &preconditioner);

  //          set preconditioner options
//   if (debug)  printf("lsolve: calling HYPRE_StructPFMGSet*\n");
  HYPRE_StructPFMGSetMaxIter(preconditioner, sol_maxit/5);
//   HYPRE_StructPFMGSetMaxIter(preconditioner, 10);
//   HYPRE_StructPFMGSetRelChange(preconditioner, sol_relch);
  HYPRE_StructPFMGSetRelaxType(preconditioner, sol_rlxtype);
  HYPRE_StructPFMGSetNumPreRelax(preconditioner, sol_npre);
  HYPRE_StructPFMGSetNumPostRelax(preconditioner, sol_npost);
  HYPRE_StructPFMGSetPrintLevel(preconditioner, sol_printl);
  HYPRE_StructPFMGSetLogging(preconditioner, sol_log);
//    if (delta != 0.0)   
//     HYPRE_StructPFMGSetTol(preconditioner, Eflt64(delta*10));
//   if (sol_zeroguess)  HYPRE_StructPFMGSetZeroGuess(preconditioner);

  //          set solver options
//   if (debug)  printf("lsolve: calling HYPRE_StructPCGSet*\n");
  if (rank > 1) {
    HYPRE_StructPCGSetMaxIter(solver, sol_maxit);
    HYPRE_StructPCGSetPrecond(solver, 
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSolve,  
		     (HYPRE_PtrToStructSolverFcn) HYPRE_StructPFMGSetup, 
		      preconditioner);
  }
  else {    // ignore pfmg preconditioner for 1D tests (bug); increase CG its
    HYPRE_StructPCGSetMaxIter(solver, sol_maxit*500);
  }
  if (delta != 0.0)   HYPRE_StructPCGSetTol(solver, Eflt64(delta));
  HYPRE_StructPCGSetup(solver, P, rhsvec, solvec);

  //       solve the linear system
//   if (debug)  printf("lsolve: calling HYPRE_StructPCGSolve\n");
  HYPRE_StructPCGSolve(solver, P, rhsvec, solvec);

  //       extract solver & preconditioner statistics
  Eflt64 finalresid=1.0;
  Eint32 Sits=0;
  Eint32 Pits=0;
//   if (debug)  printf("lsolve: calling HYPRE_StructPCGGet*\n");
  HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &finalresid);
  HYPRE_StructPCGGetNumIterations(solver, &Sits);
//   if (debug)  printf("lsolve: calling HYPRE_StructPFMGGet*\n");
  HYPRE_StructPFMGGetNumIterations(preconditioner, &Pits);
  totIters += Sits;
  if (debug)
    printf("   HYPRE resid = %g (tol = %g), PCG = %i, PFMG = %i\n",
	   finalresid,delta,Sits,Pits);

  //       extract values from solution vector 
//   if (debug)  printf("lsolve: calling HYPRE_StructVectorGetBoxValues\n");
  for (iz=SolvIndices[2][0]; iz<=SolvIndices[2][1]; iz++) {
    Zbl = (iz+zBuff)*ArrDims[0]*ArrDims[1];
    ilower[2] = iz;  iupper[2] = iz;
    for (iy=SolvIndices[1][0]; iy<=SolvIndices[1][1]; iy++) {
      Ybl = (iy+yBuff)*ArrDims[0];
      ilower[1] = iy;  iupper [1] = iy;
      HYPRE_StructVectorGetBoxValues(solvec, ilower, iupper, HYPREbuff);
      for (ix=0; ix<=SolvIndices[0][1]-SolvIndices[0][0]; ix++) 
	s_E[Zbl+Ybl+xBuff+ix] = HYPREbuff[ix];
    }
  }

  //       destroy HYPRE matrix, vector and solver structures
  HYPRE_StructPCGDestroy(solver);
  HYPRE_StructPFMGDestroy(preconditioner);

#else  // ifdef USE_HYPRE

  ENZO_FAIL("gFLDProblem_lsolve ERROR: module requires USE_HYPRE to be set!");
  
#endif


//   if (debug)  printf("gFLDProblem::lsolve -- performing step 7\n");

  //////////////////////////////////////////////////////////////
  // step (7) Set s_m = c_m - y_m*s_E  (Note: c_m stored in b_m)
  float *s_m;
  for (iz=ghZl; iz<ghZl+vsz[2]; iz++) {
    for (iy=ghYl; iy<ghYl+vsz[1]; iy++) {
      for (ix=ghXl; ix<ghXl+vsz[0]; ix++) {
	idx = (iz*ArrDims[1] + iy)*ArrDims[0] + ix;
	for (irow=1; irow<(2+Nchem); irow++) {
	  b_m = b->GetData(irow);
	  y_m = yvec->GetData(irow);
	  s_m = s->GetData(irow);
	  s_m[idx] = b_m[idx] - y_m[idx]*s_E[idx];
	}
      }
    }
  }
//   if (debug)  printf("gFLDProblem::lsolve -- finished!\n");

  // return success
  return SUCCESS;
}
#endif
