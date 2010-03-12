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
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified:   
/
/  PURPOSE: This class defines problem-specific functions for an
/           implicit free-streaming radiation solve.
/
************************************************************************/

#ifdef TRANSFER
#ifndef FS_IMPLICIT_PROBLEM_DEFINED__
#define FS_IMPLICIT_PROBLEM_DEFINED__

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EnzoVector.h"
#include "ImplicitProblemABC.h"


class FSProb : public virtual ImplicitProblemABC {

 private:

  // overall time spent in solver
  float FStime;

  // HYPRE Struct-specific data
  Eint32 mattype;                // HYPRE matrix type for solve
  Eint32 stSize;                 // stencil size
#ifdef USE_HYPRE
  HYPRE_StructGrid grid;         // HYPRE grid object for setup
  HYPRE_StructStencil stencil;   // stencil object
#endif

  // HYPRE Solver-specific data
  Eflt32 sol_tolerance;          // desired solver tolerance
  Eint32 sol_maxit;              // maximum number of iterations
  Eint32 sol_rlxtype;            // relaxation type:
                                 //    0,1 -> weighted Jacobi
                                 //    2,3 -> red-black Gauss-Seidel
  Eint32 sol_npre;               // num. pre-relaxation sweeps
  Eint32 sol_npost;              // num. post-relaxation sweeps
  Eint32 sol_printl;             // print output level
  Eint32 sol_log;                // amount of logging
  Eint32 SolvIndices[3][2];      // L/R edge indices of subdomain in global mesh
                                 // Note: these INCLUDE Dirichlet zones, even 
                                 //   though those are not included as active 
                                 //   data in the vectors or physics routines.
  int SolvOff[3];                // offset between HYPRE mesh and active mesh; 
                                 //   typically 0, but will be 1 for inclusion 
                                 //   of Dirichlet zones in HYPRE grid.

  // HYPRE interface temporary data
#ifdef USE_HYPRE
  HYPRE_StructMatrix J;          // holds free-streaming matrix
  HYPRE_StructVector rhsvec;     // holds free-streaming rhs vector
  HYPRE_StructVector solvec;     // holds free-streaming solution vector
#endif
  Eflt64 *matentries;            // holds matrix entries
  Eflt64 *rhsentries;            // linear system right-hand side
  Eflt64 *HYPREbuff;             // holds contiguous sections of solution

  // HYPRE solver diagnostics
  int totIters;                  // total MG iterations for solves

  // General problem grid information
  bool OnBdry[3][2]; // denotes if proc owns piece of boundary
  int rank;          // Rank of problem
  int layout[3];     // number of procs in each dim (1-based)
  int location[3];   // location of this proc in each dim (0-based)
  int NBors[3][2];   // process IDs of L/R neighbors in each dim
  int LocDims[3];    // implicit problem local dims (no ghost or bdry cells)
  int ArrDims[3];    // local array sizes (includes ghost and bdry cells)
  int GhDims[3][2];  // ghost cells at each face (includes Dirichlet bdry zones)
  int GlobDims[3];   // implicit problem global dimensions (active cells only)
  float dx[3];             // mesh size in each dimension
  float EdgeVals[3][2];    // L/R edges of this proc's subdomain
  float *BdryVals[3][2];   // boundary condition values 

  // time-stepping related data
  int initial_guess;   // method for initial guess calculation
  float maxdt;         // maximum FS radiation time step size
  float tnew;          // new time
  float told;          // old time
  float dt;            // time step size
  float dt_suggest;    // suggested time step size for next iteration
  float theta;         // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  float kappa0;        // background opacity
  int kappa_h2on;      // spatially dependent opacity (1=on, 0=off)
  int LimType;         // flux limiter formulation:
                       //    0 -> standard Levermore-Pomraning limiter (LP, 1981)
                       //    1 -> rational approx. to LP limiter (LP, 1981)
                       //    2 -> Reynolds approx. to LP limiter
                       //    3 -> no limiter
                       //    4 -> ZEUS limiter (like 1, but no 'albedo')
  EnzoVector *U0;      // old time-level state
  EnzoVector *extsrc;  // temporary vector holding external forcing sources
  EnzoVector *sol;     // linear system solution
  EnzoVector *kappa;   // Spatially dependent opacity

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient (new time)
  FLOAT a0;            // cosmology expansion coefficient (old time)
  FLOAT adot;          // time-derivative of a (new time)
  FLOAT adot0;         // time-derivative of a (old time)
  float aUnits;        // expansion parameter scaling
  float EScale;        // radiation energy density scaling factor
  float EUnits;        // radiation energy density unit conversion factor (new)
  float EUnits0;       // radiation energy density unit conversion factor (old)
  float LenUnits;      // length scaling factor (new)
  float LenUnits0;     // length scaling factor (old)
  float TimeUnits;     // time scaling factor (new)
  float TimeUnits0;    // time scaling factor (old)
  float DenUnits;      // density scaling factor (new)
  float DenUnits0;     // density scaling factor (old)

  // ionization parameters
  Eflt64 NGammaDot;     // ionization strength (photons/sec)
  float EtaRadius;     // ionization source radius
  float EtaCenter[3];  // ionization source location

  // private computation routines
  int EnforceBoundary(EnzoVector *vec);
  int SetupSystem(Eflt64 *mat, Eflt64 *rhs, float *rhsnorm, 
		  float *E, float *E0, float *eta, float *opacity);
  int RadiationSource(float *Efsrc);
  int InitialGuess(EnzoVector *Ef, EnzoVector *Ef0, EnzoVector *Efsrc);


 public:

  // boundary type in each dimension, face:
  //    0->periodic
  //    1->dirichlet
  //    2->neumann
  int BdryType[3][2];

  ///////////////////////////////////////
  // Module-Specific Routines

  // Constructor
  FSProb();
  
  // Destructor
  ~FSProb();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Re-sets emissivity source magnitude (look in FSProb_Initialize.C)
  int SetEmissivity(float NGammaDot);
  
  // Problem setup/solver
  int Evolve(HierarchyEntry *ThisGrid, float deltat);

  // Write module parameters to file
  int WriteParameters(FILE *fptr);

  // Problem debug (for output upon failure)
  int Dump(EnzoVector *ucur);
  
  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int BdryConst, float *BdryData);

  int ComputeOpacityLW(float *H2Density);

};


#endif
#endif
