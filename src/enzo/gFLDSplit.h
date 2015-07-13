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
/  Single-Group, Multi-species, Gray Flux-Limited Diffusion Implicit 
/  Problem Class
/
/  written by: Daniel Reynolds
/  date:       July 2009
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
#ifndef FLD_SPLIT_PROBLEM_DEFINED__
#define FLD_SPLIT_PROBLEM_DEFINED__

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
#include "ImplicitProblemABC.h"


class gFLDSplit : public virtual ImplicitProblemABC {

 private:
  
  // overall time spent in solver and components
  float RTtime;
  float HYPREtime;
  float ChemTime;
  
  // HYPRE Struct-specific data
  Eint32 mattype;                // HYPRE matrix type for solve
  Eint32 stSize;                 // stencil size
#ifdef USE_HYPRE
  HYPRE_StructGrid grid;         // HYPRE grid object for setup
  HYPRE_StructStencil stencil;   // stencil object
#endif

  // HYPRE Solver-specific data
  float  sol_tolerance;          // desired solver tolerance
  Eint32 sol_maxit;              // maximum number of iterations
  Eint32 sol_rlxtype;            // relaxation type:
                                 //    0,1 -> weighted Jacobi
                                 //    2,3 -> red-black Gauss-Seidel
  Eint32 sol_npre;               // num. pre-relaxation sweeps
  Eint32 sol_npost;              // num. post-relaxation sweeps
  Eint32 sol_printl;             // print output level
  Eint32 sol_log;                // amount of logging
  int    Krylov_method;          // flag denoting which outer solver to use:
                                 //    0 => PCG
                                 //    1 => BiCGStab (default)
                                 //    2 => GMRES
  Eint32 SolvIndices[3][2];      // L/R edge indices of subdomain in global mesh
                                 // Note: these INCLUDE Dirichlet zones, even 
                                 //   though those are not included as active 
                                 //   data in the vectors or physics routines.
  int SolvOff[3];                // offset between HYPRE mesh and active mesh; 
                                 //   typically 0, but will be 1 for inclusion 
                                 //   of Dirichlet zones in HYPRE grid.

  // HYPRE interface temporary data
#ifdef USE_HYPRE
  HYPRE_StructMatrix P;          // holds radiation matrix
  HYPRE_StructVector rhsvec;     // holds radiation rhs vector
  HYPRE_StructVector solvec;     // holds radiation solution vector
#endif
  Eflt64 *matentries;            // holds radiation matrix entries
  Eflt64 *rhsentries;            // linear system rhs entries
  Eflt64 *HYPREbuff;             // holds contiguous sections of rhs/sol

  // HYPRE solver diagnostics
  int totIters;                  // total MG iterations for solves

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
  float *BdryVals[3][2];   // boundary values for radiation BCs

  // time-stepping related data
  int initial_guess;   // parameter for setting the initial guess:
  float initdt;        // initial radiation time step size
  float maxdt;         // maximum radiation/chemistry/heating time step size
  float mindt;         // minimum radiation/chemistry/heating time step size
  float maxsubcycles;  // max subcycle factor for rad time step within hydro step
  float maxchemsub;    // max subcycle factor for chem time step within rad step
  float dtfac[3];      // desired relative change in fields per step
  float dtnorm;        // norm choice for computing relative change:
                       //    0 -> max pointwise norm (default)
                       //   >0 -> rms p-norm over entire domain
  float dtgrowth;      // time step growth factor (1 < dtgrowth < 10)
  float tnew;          // new time
  float told;          // old time
  float dt;            // time step size
  float dtrad;         // radiation time step size (subcycled)
  float dtchem;        // chemistry/gas time step size (subcycled)
  float theta;         // implicitness parameter (1->BE, 0.5->CN, 0->FE)
  EnzoVector *sol;     // solution vector
  EnzoVector *U0;      // old time-level state
  EnzoVector *extsrc;  // temporary vector holding external forcing sources
  
  // problem defining data
  int Nchem;           // number of chemical species (non-negative integer)
  int Model;           // model choice:
                       //    1 => case B HII recombination rates
                       //    4 => isothermal, case B HII recombination
                       //   10 => LTE couplings (ZEUS), constant opacity
  float NGammaDot;     // ionization strength (photons/sec)
  float EtaRadius;     // ionization source radius
  float EtaCenter[3];  // ionization source location

  // cosmology and scaling constants
  FLOAT a;             // cosmology expansion coefficient
  FLOAT a0;            // cosmology expansion coefficient (old time)
  FLOAT adot;          // time-derivative of a
  FLOAT adot0;         // time-derivative of a (old time)
  float aUnits;        // expansion parameter scaling
  float ErScale;       // radiation energy density scaling factor
  float ecScale;       // specific energy correction scaling factor
  float NiScale;       // species density scaling factor
  bool  autoScale;     // flag to enable/disable automatic scaling factors
  bool  StartAutoScale;  // flag to turn begin automatic scaling in a run
  float ErUnits;       // radiation energy density unit conversion factor
  float ecUnits;       // specific energy correction unit conversion factor
  float NiUnits;       // species density unit conversion factor
  float ErUnits0;      // radiation energy density unit conversion factor
  float NiUnits0;      // species density unit conversion factor

  float DenUnits;      // density scaling factor
  float LenUnits;      // length scaling factor
  float TimeUnits;     // time scaling factor
  float VelUnits;      // velocity scaling factor
  float DenUnits0;     // density scaling factor
  float LenUnits0;     // length scaling factor

  // chemistry constants
  float HFrac;         // Fraction of matter composed of Hydrogen

  // opacity constants
  //   opacity computed as C0 * (rho/C1)^C2
  float EnergyOpacityC0;
  float EnergyOpacityC1;
  float EnergyOpacityC2;

  // storage for integrals over radiation spectrum (set during initialization)
  float hnu0_HI;            // HI ionization threshold (eV)
  float hnu0_HeI;           // HeI ionization threshold (eV)
  float hnu0_HeII;          // HeII ionization threshold (eV)
  int ESpectrum;            // integer flag determining spectrum choice, 
                            // negative values imply monochromatic SED
                            //   1 -> 1e5 black body spectrum
                            //   0 -> simple power law spectrum
                            //  -1 -> monochromatic spectrum @ hnu0_HI
                            //  -2 -> monochromatic spectrum @ hnu0_HeI
                            //  -3 -> monochromatic spectrum @ hnu0_HeII
  float intSigE;            // int_{nu0}^{inf} sigma_E(nu) d nu
  float intSigESigHI;       // int_{nu0}^{inf} sigma_E(nu)*sigma_HI(nu) d nu
  float intSigESigHeI;      // int_{nu0}^{inf} sigma_E(nu)*sigma_HeI(nu) d nu
  float intSigESigHeII;     // int_{nu0}^{inf} sigma_E(nu)*sigma_HeII(nu) d nu
  float intSigESigHInu;     // int_{nu0}^{inf} sigma_E(nu)*sigma_HI(nu)/nu d nu
  float intSigESigHeInu;    // int_{nu0}^{inf} sigma_E(nu)*sigma_HeI(nu)/nu d nu
  float intSigESigHeIInu;   // int_{nu0}^{inf} sigma_E(nu)*sigma_HeII(nu)/nu d nu

  // access to Enzo data
  float *vx;         // x0-directional velocity
  float *vy;         // x1-directional velocity
  float *vz;         // x2-directional velocity
  float *rho;        // density
  float *eh;         // fluid energy (total or internal, see DualEnergyFormalism)

  // stored arrays for increased efficiency
  float *Temperature;              // gas 'Temperature' (for black-body radiation)
  float *Temperature0;             // gas 'Temperature' (old time)
  float *FluidEnergyCorrection;    // gas energy correction
  float *OpacityE;                 // Energy mean opacity

  // private computation routines
  int EnforceBoundary(EnzoVector *u);
  int ComputeTemperature(float *Temperature, EnzoVector *u);
  int SetupSystem(Eflt64 *matentries, Eflt64 *rhsentries, float *rhsnorm, 
		  float *Eg0, float *Eg, float *OpacityE, float *Temp, 
		  float *Temp0, float *eta);
  int Opacity(float *OpacityE, float *time, EnzoVector *u);
  int RadiationSource(float *Egsrc, float *time);
  int GasEnergySource(float *ecsrc, float *time);
  int ChemistrySource(float *HIsrc, float *HeIsrc, float *HeIIsrc, float *time);
  float RadiationSpectrum(float nu);
  float CrossSections(float nu, int species);
  int ComputeRadiationIntegrals();
  int AnalyticInitGuess(EnzoVector *u, float dt);
  int AnalyticChemistry(EnzoVector *u0, EnzoVector *u, EnzoVector *src, float dt);
  int FillRates(EnzoVector *u, EnzoVector *u0, float *phHI, float *phHeI, 
		float *phHeII, float *PhotoGamma, float *dissH2I);
  int RadStep(HierarchyEntry *ThisGrid, int eta_set);
  int ChemStep(HierarchyEntry *ThisGrid, float thisdt, float tcur);
  int ChemBounds(HierarchyEntry *ThisGrid);
  

 public:

  // boundary type in each dimension, face:
  //    0->periodic
  //    1->dirichlet
  //    2->neumann
  int BdryType[3][2];

  ///////////////////////////////////////
  // FLD-Specific Routines

  // Constructor
  gFLDSplit();
  
  // Destructor
  ~gFLDSplit();

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

  // Fill in initial guess for time-evolved solution
  int InitialGuess(EnzoVector *uvec);
  
  // Return the maximum rad-hydro time step size
  float ComputeTimeStep(EnzoVector *uold, EnzoVector *unew, int flag);

};


#endif
#endif
