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
/  Single-Group, Multi-species, Gray Flux-Limited Diffusion Implicit 
/  Problem Class
/
/  written by: Daniel Reynolds
/  date:       August, 2006
/  modified1:  June 12, 2007 by John Hayes; added MarshakParms pointer
/
/  PURPOSE: This class defines problem-specific functions for an 
/           implicit gray flux-limited diffusion solve.
/
/           The variables are stored in the following order: 
/              0 -> radiation energy density
/              1 -> fluid energy correction
/              2:Nspecies+1 -> chemical species (Nspecies may be 0)
/
************************************************************************/

#ifdef TRANSFER
#ifndef FLD_NONLINEAR_PROBLEM_DEFINED__
#define FLD_NONLINEAR_PROBLEM_DEFINED__

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
#include "InexactNewton.h"
#include "NonlinearProblemABC.h"


class gFLDProblem : public virtual NonlinearProblemABC {

 private:
  
  // flag denoted problem preparedness
  bool prepared;

  // overall time spent in solver
  float RTtime;
  
  // HYPRE Struct-specific data
  Eint32 mattype;                // HYPRE matrix type for solve
  Eint32 stSize;                 // stencil size
#ifdef USE_HYPRE
  HYPRE_StructGrid grid;         // HYPRE grid object for setup
  HYPRE_StructStencil stencil;   // stencil object
#endif

  // HYPRE Solver-specific data
  Eint32 sol_zeroguess;          // use a zero initial guess
  Eint32 sol_maxit;              // maximum number of iterations
  Eint32 sol_relch;              // relative change stopping criteria
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
  HYPRE_StructMatrix P;          // holds Schur complement matrix
  HYPRE_StructVector rhsvec;     // holds Schur complement rhs vector
  HYPRE_StructVector solvec;     // holds Schur complement solution vector
#endif
  Eflt64 *Ptmpvec;               // holds Schur complement matrix entries
  Eflt64 *HYPREbuff;             // holds contiguous sections of rhs/sol

  // HYPRE solver diagnostics
  int totIters;                  // total MG iterations for solves

  // Inexact Newton solver-specific data
  InexactNewtonSolver *INSolve;  // Inexact Newton solver
  int approx_jac;                // param to approximate the local jacobian: 
                                 //    0 -> use the analytical jac (default)
                                 //    1 -> approximate the jac
  int initial_guess;             // parameter for setting the initial guess:
                                 //    0 -> use previous time step
                                 //    1 -> full fwd Euler (local sources)
                                 //    2 -> full fwd Euler (all rhs)
                                 //    3 -> partial fwd Euler (local sources)
                                 //    4 -> partial fwd Euler (all rhs)
                                 //    5 -> local analytic guess (Rx only)
  int newt_linesearch;           // use a linesearch in the Newton method
                                 //    0 -> no linesearch
                                 //  !=0 -> linesearch (i.e. damped Newton)
  int newt_maxit;                // maximum number of iterations
  int newt_norm;                 // norm for convergence measurement
  float newt_INconst;            // Inexact-Newton constant
  float newt_tol;                // Newton tolerance
  float newt_MinLinesearch;      // minimum allowed line-search length

  // General problem grid information
  bool OnBdry[3][2]; // denotes if proc owns piece of boundary
  int rank;          // Rank of self-gravity problem
  int layout[3];     // number of procs in each dim (1-based)
  int location[3];   // location of this proc in each dim (0-based)
  int NBors[3][2];   // process IDs of L/R neighbors in each dim
  int LocDims[3];    // implicit problem local dims (no ghost or bdry cells)
  int ArrDims[3];    // local array sizes (includes ghost and bdry cells)
  int GhDims[3][2];  // ghost cells at each face (includes Dirichlet bdry zones)
  int GlobDims[3];   // implicit problem global dimensions (active cells only)
  float dx[3];             // mesh size in each dimension
  float EdgeVals[3][2];    // L/R edges of this proc's subdomain
  float *EBdryVals[3][2];  // boundary values for radiation BCs
  float *FBdryVals[3][2];  // additional boundary values for mixed BCs

  // time-stepping related data
  float initdt;        // initial radiation time step size
  float maxdt;         // maximum radiation time step size
  float mindt;         // minimum radiation time step size
  float dtfac[3];      // desired relative change in fields per step
  float dtnorm;        // norm choice for computing relative change:
                       //    0 -> max pointwise norm (default)
                       //   >0 -> rms p-norm over entire domain
  int dtsolver;        // enable solver heuristic-based time step selection
  float tnew;          // new time
  float told;          // old time
  float dt;            // time step size
  float theta;         // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  int LimType;         // flux limiter formulation:
                       //    0 -> standard Levermore-Pomraning limiter (LP, 1981)
                       //    1 -> rational approx. to LP limiter (LP, 1981)
                       //    2 -> Larsen n=2 limiter
                       //    3 -> no limiter
                       //    4 -> ZEUS limiter (like 1, but no 'albedo')
  EnzoVector *sol;     // solution vector
  EnzoVector *U0;      // old time-level state
  EnzoVector *rhs;     // current time-level rhs
  EnzoVector *rhs0;    // old time-level rhs
  EnzoVector *extsrc;  // temporary vector holding external forcing sources
  EnzoVector *tmp1;    // temporary 
  EnzoVector *tmp2;    // temporary 
  EnzoVector *tmp3;    // temporary 

  int AnalyticChem;    // use analytical reaction solver instead of theta-method
  
  // problem defining data
  float *UTypVals;     // constants giving typical values for each species
  int Nchem;           // number of chemical species (non-negative integer)
  int Model;           // model choice, 0=>decoupled ODE test case
                       //               1=>case B HII recomb, no emissivity
                       //               2=>case A HII recomb, with emissivity
                       //               3=>coupled ODE test case
                       //               4=>isothermal ionization, pt. src. rad.
                       //               5=>isothermal model, pt. source radiation
                       //              10=>ZEUS couplings, constant opacity
                       //              20=>Marshak-style (Su & Olson) coupling
                       //                  specific heat ~ T^3; 
                       //                  constant opacity
  float *MarshakParms; // parameters used to configure Marshak-type problems
  float IonizationParms[5];  // parameters for configuring ionization problems

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient
  FLOAT adot;          // time-derivative of a
  float aUnits;        // expansion parameter scaling
  float ErScale;       // radiation energy density scaling factor
  float ecScale;       // specific energy correction scaling factor
  float NiScale;       // species density scaling factor
  float ErUnits;       // radiation energy density unit conversion factor
  float ecUnits;       // specific energy correction unit conversion factor
  float NiUnits;       // species density unit conversion factor

  float DenUnits;      // density scaling factor
  float LenUnits;      // length scaling factor
  float TimeUnits;     // time scaling factor
  float TempUnits;     // temperature scaling factor
  float VelUnits;      // velocity scaling factor
  double MassUnits;    // mass scaling factor

  // chemistry constants
  float HFrac;         // Fraction of matter composed of Hydrogen

  // opacity constants
  //   opacities computed as C0 * (rho/C1)^C2 & (T/C3)^C4
  float PlanckOpacityC0;
  float PlanckOpacityC1;
  float PlanckOpacityC2;
  float PlanckOpacityC3;
  float PlanckOpacityC4;
  float EnergyOpacityC0;
  float EnergyOpacityC1;
  float EnergyOpacityC2;
  float EnergyOpacityC3;
  float EnergyOpacityC4;

  // storage for integrals over radiation spectrum (set during initialization)
  float hnu0_HI;            // HI ionization threshold (eV)
  float hnu0_HeI;           // HeI ionization threshold (eV)
  float hnu0_HeII;          // HeII ionization threshold (eV)
  int ESpectrum;            // integer flag determining spectrum choice
                            //   1 -> 1e5 black body spectrum
                            //   0 -> simple power law spectrum
                            //  -1 -> monochromatic spectrum
  float intSigE;            // int_{nu0}^{inf} sigma_E(nu) d nu
  float intSigESigHI;       // int_{nu0}^{inf} sigma_E(nu)*sigma_HI(nu) d nu
  float intSigESigHeI;      // int_{nu0}^{inf} sigma_E(nu)*sigma_HeI(nu) d nu
  float intSigESigHeII;     // int_{nu0}^{inf} sigma_E(nu)*sigma_HeII(nu) d nu
  float intSigESigHInu;     // int_{nu0}^{inf} sigma_E(nu)*sigma_HI(nu)/nu d nu
  float intSigESigHeInu;    // int_{nu0}^{inf} sigma_E(nu)*sigma_HeI(nu)/nu d nu
  float intSigESigHeIInu;   // int_{nu0}^{inf} sigma_E(nu)*sigma_HeII(nu)/nu d nu

  // linear solver/Jacobian arrays
  EnzoVector **L;    // local Jacobian components 

  // access to Enzo data
  float *vx;         // x0-directional velocity
  float *vy;         // x1-directional velocity
  float *vz;         // x2-directional velocity
  float *rho;        // density
  float *eh;         // fluid energy (total or internal, depending on DualEnergyFormalism)

  // stored arrays for increased efficiency
  float *FluidEnergyCorrection;    // gas energy correction
  float *Temp;      // gas temperature
  float *OpacityP;  // Planck mean opacity
  float *OpacityE;  // Energy mean opacity

  // private computation routines
  int LocRHS(EnzoVector *locrhs, float time, EnzoVector *u);
  int ComputeRHS(EnzoVector *rhsval, float time, EnzoVector *u);
  int ComputeTemperature(float *Temperature, float time, 
			 FLOAT a, EnzoVector *u);
  int MatrixEntries(Eflt64 *matentries, float *Eg, float *Eg0, 
		    float *Temperature, float *OpacityE, float *adjvec);
  int SetNewtonBCs(Eflt64 *matentries, float *rhsentries);
  int DiffRHS(float *drhs, float *Eg, float *Eg0, float *Temperature, 
	      float *OpacityE);
  int LocalRHS(float *Egrhs, float *ecrhs, float *HIrhs, float *HeIrhs, 
	       float *HeIIrhs, float *Egsrc, float *ecsrc, float *HIsrc, 
	       float *HeIsrc, float *HeIIsrc, float *time, float *ec, 
	       float *Eg, float *Temperature, float *kappaP, 
	       float *kappaE, float *nHI, float *nHeI, float *nHeII);
  int LocalJac(float *Egjac_Eg, float *Egjac_ec, float *Egjac_HI, 
	       float *Egjac_HeI, float *Egjac_HeII, float *ecjac_Eg, 
	       float *ecjac_ec, float *ecjac_HI, float *ecjac_HeI, 
	       float *ecjac_HeII, float *HIjac_Eg, float *HIjac_ec, 
	       float *HIjac_HI, float *HIjac_HeI, float *HIjac_HeII, 
	       float *HeIjac_Eg, float *HeIjac_ec, float *HeIjac_HI, 
	       float *HeIjac_HeI, float *HeIjac_HeII, float *HeIIjac_Eg, 
	       float *HeIIjac_ec, float *HeIIjac_HI, float *HeIIjac_HeI, 
	       float *HeIIjac_HeII, float *time, float *ec, float *Eg, 
	       float *nHI, float *nHeI, float *nHeII);
  int BlockSolve(float *Amat, float *xvec, float *bvec, int *N, int *M);
  int Opacity(float *OpacityP, float *OpacityE, float *time, float *n_HI, 
	      float *n_HeI, float *n_HeII, float *Temperature);
  int RadiationSource(float *Egsrc, float *time, float *Eg, float *ec, 
		      float *n_HI, float *n_HeI, float *n_HeII, 
		      float *Temperature);
  int GasEnergySource(float *ecsrc, float *time, float *Eg, float *ec, 
		      float *n_HI, float *n_HeI, float *n_HeII, 
		      float *Temperature);
  int ChemistrySource(float *HIsrc, float *HeIsrc, float *HeIIsrc, 
		      float *time, float *Eg, float *ec, float *n_HI, 
		      float *n_HeI, float *n_HeII, float *Temperature);
  float RadiationSpectrum(float nu);
  float CrossSections(float nu, int species);
  int ComputeRadiationIntegrals();
  int AnalyticInitGuess(EnzoVector *u, float dt);
  int AnalyticInitGuess2(EnzoVector *u, float dt);
  int AnalyticResid(EnzoVector *u0, EnzoVector *u, EnzoVector *fu, float dt);


 public:

  // boundary type in each dimension, face:
  //    0->periodic
  //    1->dirichlet
  //    2->neumann
  int BdryType[3][2];

  ///////////////////////////////////////
  // FLD-Specific Routines

  // Constructor
  gFLDProblem();
  
  // Destructor
  ~gFLDProblem();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem Evolver
  int Evolve(HierarchyEntry *ThisGrid, float deltat);
  
  // Write module parameters to file
  int WriteParameters(FILE *fptr);

  // Problem debug (for output upon failure)
  int Dump(EnzoVector *ucur);
  
  // Problem Boundary Condition setup (called once or at each time step, 
  //    must be called for each locally-owned external face separately)
  int SetupBoundary(int Dimension, int Face, int BdryConst, float *BdryData);

  // Enforce boundary conditions onto a vector
  int EnforceBoundary(EnzoVector *vec, int flag);

  // Update boundary conditions (for time-dependent BCs)
  int UpdateBoundary(EnzoVector *vec, float time, int flag);

  // Fill in initial guess for time-evolved solution
  int InitialGuess(EnzoVector *uvec);
  
  // Return the maximum rad-hydro time step size
  float ComputeTimeStep(EnzoVector *uold, EnzoVector *unew, 
			int NewtIts, float FStep, float FResid);

  ///////////////////////////////////////
  // Nonlinear Solver Interface Routines
  
  // Problem-defining nonlinear residual operations (called repeatedly)
  int nlresid(EnzoVector *fu, EnzoVector *u);
  
  // Problem-specific Linear system setup function, sets up the 
  //   linear Newton system matrix J(u) given an updated state u
  //   (called once per Newton iteration)
  int lsetup(EnzoVector *u);
  
  // Problem-specific Linear solver function 
  //   solves J(u)*s = b to tolerance delta
  //   (called once per Newton iteration)
  int lsolve(EnzoVector *s, EnzoVector *b, EnzoVector *u, float delta);
  
};


#endif
#endif
