#ifndef __typedefs_h_
#define __typedefs_h_
/***********************************************************************
/
/  MISCELANEOUS TYPEDEFS AND ENUMERATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include "CloudyCoolingData.h"
#include "CoolData.h"
#include "RateData.h"
#include "RadiationFieldData.h"
#include "TestProblemData.h"
#include "IndividualStarData.h"
#include "StellarYieldsData.h"

/* These are the different types of baryon fields. */

#ifdef SMALL_INTS
typedef int field_type;
typedef int boundary_type;
typedef int gravity_boundary_type;
typedef int interpolation_type;
typedef int hydro_method;
typedef int star_type;
typedef int enum_type;
typedef int staggering;
typedef int fieldtype;
typedef int mhd_ct_method;
typedef int forcing_type;
#endif

#ifdef LARGE_INTS
typedef long_int field_type;
typedef long_int boundary_type;
typedef long_int gravity_boundary_type;
typedef long_int interpolation_type;
typedef long_int hydro_method;
typedef long_int star_type;
typedef long_int enum_type;
typedef long_int staggering;
typedef long_int fieldtype;
typedef int mhd_ct_method;
typedef long_int forcing_type;
#endif

const field_type 
  Density         = 0,
  TotalEnergy     = 1,
  InternalEnergy  = 2,
  Pressure        = 3,
  Velocity1       = 4,
  Velocity2       = 5,
  Velocity3       = 6,
  ElectronDensity = 7,
  HIDensity       = 8,
  HIIDensity      = 9,
  HeIDensity      = 10,
  HeIIDensity     = 11,
  HeIIIDensity    = 12,
  HMDensity       = 13,
  H2IDensity      = 14,
  H2IIDensity     = 15,
  DIDensity       = 16,
  DIIDensity      = 17,
  HDIDensity      = 18,
  SNColour        = 19,
  Metallicity     = 20,
  ExtraType0      = 21,
  ExtraType1      = 22,
  kphHI           = 23,
  PhotoGamma      = 24,
  kphHeI          = 25,
  gammaHeI        = 26,
  kphHeII         = 27,
  gammaHeII       = 28,
  kdissH2I        = 29,
  GravPotential   = 30,
  Acceleration0   = 31,
  Acceleration1   = 32,
  Acceleration2   = 33,
  RadPressure0    = 34,
  RadPressure1    = 35,
  RadPressure2    = 36,
  Emissivity0     = 37,

/* these pseudo-fields are used to access grid data 
   the "g" prefix is to avoid namespace conflict */

  gParticlePosition     = 37,
  gParticleVelocity     = 38,
  gParticleMass         = 39,
  gParticleAcceleration = 40,
  gParticleNumber       = 41,
  gParticleType         = 42,
  gParticleAttribute    = 43,
  gPotentialField       = 44,
  gAccelerationField    = 45,
  gGravitatingMassField = 46,
  gFlaggingField        = 47,
  gVelocity             = 48,

  Bfield1               = 49,
  Bfield2               = 50,
  Bfield3               = 51,
  PhiField              = 52,
  Phi_pField            = 53,
  DebugField            = 54, 

  DrivingField1         = 55, 
  DrivingField2         = 56, 
  DrivingField3         = 57,

  AccelerationField1    = 58, 
  AccelerationField2    = 59, 
  AccelerationField3    = 60,

  Galaxy1Colour          = 61,
  Galaxy2Colour          = 62,
/* these are required for Sam Skillman's Shock/Cosmic ray models. */
  Mach            = 63,
  PreShockTemperature = 64,
  PreShockDensity = 65,  

/* these are required for Simon Glover's chemistry (which also needs some of the
   other fields, which are used for MultiSpecies) */
  CIDensity       = 66,
  CIIDensity      = 67, 
  OIDensity       = 68, 
  OIIDensity      = 69,
  SiIDensity      = 70,
  SiIIDensity     = 71,
  SiIIIDensity    = 72,
  CHIDensity      = 73,
  CH2IDensity     = 74,
  CH3IIDensity    = 75,
  C2IDensity      = 76,
  COIDensity      = 77,
  HCOIIDensity    = 78,
  OHIDensity      = 79,
  H2OIDensity     = 80,
  O2IDensity      = 81,

  MBHColour       = 82,
  ForbiddenRefinement = 83,

/* FLD radiation module stuff (D. Reynolds) */ 
  RadiationFreq0  = 84,
  RadiationFreq1  = 85,
  RadiationFreq2  = 86,
  RadiationFreq3  = 87,
  RadiationFreq4  = 88,
  RadiationFreq5  = 89,
  RadiationFreq6  = 90,
  RadiationFreq7  = 91,
  RadiationFreq8  = 92,
  RadiationFreq9  = 93,

  /* Number of ray segments for ray tracing load balancing */
  RaySegments     = 94,

  /* Metals from Type Ia SNe */
  MetalSNIaDensity = 95,
  MetalSNIIDensity = 96,

  /* Cosmic Ray Energy Density */
  CRDensity = 97,

  /* Chemical Evolution Element Densities */
  LiDensity = 98,
  BeDensity = 99,
  BDensity = 100,
  CDensity = 101,
  NDensity = 102,
  ODensity = 103,
  FDensity = 104,
  NeDensity = 105,
  NaDensity = 106,
  MgDensity = 107,
  AlDensity = 108,
  SiDensity = 109,
  PDensity = 110,
  SDensity = 111,
  ClDensity = 112,
  ArDensity = 113,
  KDensity = 114,
  CaDensity = 115,
  ScDensity = 116,
  TiDensity = 117,
  VDensity = 118,
  CrDensity = 119,
  MnDensity = 120,
  FeDensity = 121,
  CoDensity = 122,
  NiDensity = 123,
  CuDensity = 124,
  ZnDensity = 125,
  GaDensity = 126,
  GeDensity = 127,
  AsDensity = 128,
  SeDensity = 129,
  BrDensity = 130,
  KrDensity = 131,
  RbDensity = 132,
  SrDensity = 133,
  YDensity = 134,
  ZrDensity = 135,
  NbDensity = 136,
  MoDensity = 137,
  TcDensity = 138,
  RuDensity = 139,
  RhDensity = 140,
  PdDensity = 141,
  AgDensity = 142,
  CdDensity = 143,
  InDensity = 144,
  SnDensity = 145,
  SbDensity = 146,
  TeDensity = 147,
  IDensity = 148,
  XeDensity = 149,
  CsDensity = 150,
  BaDensity = 151,
  LaDensity = 152,
  CeDensity = 153,
  PrDensity = 154,
  NdDensity = 155,
  PmDensity = 156,
  SmDensity = 157,
  EuDensity = 158,
  GdDensity = 159,
  TbDensity = 160,
  DyDensity = 161,
  HoDensity = 162,
  ErDensity = 163,
  TmDensity = 164,
  YbDensity = 165,
  LuDensity = 166,
  HfDensity = 167,
  TaDensity = 168,
  WDensity = 169,
  ReDensity = 170,
  OsDensity = 171,
  IrDensity = 172,
  PtDensity = 173,
  AuDensity = 174,
  HgDensity = 175,
  TlDensity = 176,
  PbDensity = 177,
  BiDensity = 178,

  PeHeatingRate = 179,
  OTLWkdissH2I  = 180,

  FieldUndefined  = 181;

/*
enum field_type {Density, TotalEnergy, InternalEnergy, Pressure,
		 Velocity1, Velocity2, Velocity3, 
		 ElectronDensity, HIDensity, HIIDensity,  HeIDensity, 
		 HeIIDensity, HeIIIDensity, HMDensity, H2IDensity, 
		 H2IIDensity, DIDensity, DIIDensity, HDIDensity,
                 Metallicity, ExtraType0, ExtraType1, GravPotential,
		 Acceleration0, Acceleration1,Acceleration2,
                 FieldUndefined};
*/

#define FieldTypeIsDensity(A) ((((A) >= TotalEnergy && (A) <= Velocity3) || ((A) >= kphHI && (A) <= kdissH2I) || ((A) >= RadiationFreq0 && (A) <= RaySegments) || ((A) >= Bfield1 && (A) <= AccelerationField3) || ((A) == PeHeatingRate) || ((A) == OTLWkdissH2I)) ? FALSE : TRUE)
#define FieldTypeIsRadiation(A) ((((A) >= kphHI && (A) <= kdissH2I) || ((A) >= RadiationFreq0 && (A) <= RadiationFreq9) || ((A) == PeHeatingRate) || ((A) == OTLWkdissH2I)) ? TRUE : FALSE)
#define FieldTypeNoInterpolate(A) (((((A) >= Mach) && ((A) <= PreShockDensity)) || ((A) == GravPotential) || ((A) == RaySegments)) || ((A) == OTLWkdissH2I) || ((A) == PeHeatingRate) ? TRUE : FALSE)

/* Different stochastic forcing types */
const forcing_type
  None       = 0,
  Peak       = 1,
  Parabolic  = 2,
  Band       = 3;

/* These are the different types of fluid boundary conditions. */

const boundary_type
  reflecting        = 0,
  outflow           = 1,
  inflow            = 2,
  periodic          = 3,
  shearing          = 4,
  BoundaryUndefined = 5;

// enum boundary_type {reflecting, outflow, inflow, periodic, BoundaryUndefined};

/* These are the different types of gravity boundary conditions. */

const gravity_boundary_type
  TopGridPeriodic  = 0,
  TopGridIsolated  = 1,
  SubGridIsolated  = 2,
  GravityUndefined = 3;

// enum gravity_boundary_type {TopGridPeriodic, TopGridIsolated, 
// 				    SubGridIsolated, GravityUndefined};

/* Interpolation types. */

const interpolation_type
  ThirdOrderA            = 0,
  SecondOrderA           = 1,
  SecondOrderB           = 2,
  SecondOrderC           = 3,
  FirstOrderA            = 4,
  InterpolationUndefined = 5;


// enum interpolation_type {ThirdOrderA, SecondOrderA, SecondOrderB, SecondOrderC,
// 			 FirstOrderA, InterpolationUndefined};

/* Hydrodynamics methods. */

const hydro_method
  PPM_DirectEuler      = 0,
  PPM_LagrangeRemap    = 1,
  Zeus_Hydro           = 2,
  HD_RK                = 3,
  MHD_RK               = 4,
  NoHydro              = 5, 
  MHD_Li             = 6,
  HydroMethodUndefined = 7;

// enum hydro_method {PPM_DirectEuler, PPM_LagrangeRemap, Zeus_Hydro};

const enum_type iHI = 0, iHeI = 1, iHeII = 2, iH2I = 3, iHII = 4;
const enum_type Cartesian = 0, Spherical = 1, Cylindrical = 2;
const enum_type PLM = 0, PPM = 1, CENO = 2, WENO3 = 3, WENO5 = 4, ZERO = 5;
const enum_type FluxReconstruction = 0, HLL = 1, Marquina = 2,
  LLF = 3, HLLC = 4, TwoShock = 5, HLLD = 6;
const enum_type Neumann = 0, Dirichlet = 1;
const enum_type Isotropic = 1, Beamed = -2, Episodic = -3;

/* Stanford RK MUSCL solvers support */ 
//enum {Cartesian, Spherical, Cylindrical};
//enum {PLM, PPM, CENO, WENO3, WENO5};
//enum {FluxReconstruction, HLL, Marquina, LLF, HLLC};

/* These are the different types of poisson cleaining boundary conditions. */
//enum{Neumann, Dirichlet};

const mhd_ct_method 
  CT_None = 0,
  CT_BalsaraSpicer = 1,
  CT_Athena_LF = 2,
  CT_Athena_Switch = 3,
  CT_Biermann = 4;

/* Definitions for streaming format */

const staggering VERTEX_CENTERED = 0, CELL_CENTERED = 1;
const fieldtype SCALAR = 1, VECTOR = 3;

/* Star particle types */

const star_type
  PopIII = PARTICLE_TYPE_SINGLE_STAR,
  PopII = PARTICLE_TYPE_CLUSTER,
  NormalStar = PARTICLE_TYPE_STAR,
  SimpleSource = PARTICLE_TYPE_SIMPLE_SOURCE,
  BlackHole = PARTICLE_TYPE_BLACK_HOLE,
  IndividualStar = PARTICLE_TYPE_INDIVIDUAL_STAR,
  IndividualStarWD = PARTICLE_TYPE_INDIVIDUAL_STAR_WD,
  IndividualStarRemnant = PARTICLE_TYPE_INDIVIDUAL_STAR_REMNANT,
  IndividualStarUnresolved = PARTICLE_TYPE_INDIVIDUAL_STAR_UNRESOLVED,
  PopIII_CF = PARTICLE_TYPE_COLOR_STAR, // Non-radiating PopIII
  MBH = PARTICLE_TYPE_MBH;

/* Define a float/int union. */

union float_int {
  long_int ival;
  PINT IVAL;
  float fval;
  FLOAT FVAL;
};

//struct ParticleMoveList {
//  int FromGrid;
//  int ToGrid[6];
//  int NumberToMove[6];
//  float_int *Pointer[6];
//};

struct particle_data {
  FLOAT pos[MAX_DIMENSION];
  float vel[MAX_DIMENSION];
  float mass;
  float attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  PINT  id;
  int   type;
  int   grid;
  int   proc;
};

#include "StarBuffer.h"
struct star_data {
  StarBuffer data;
  int grid;
  int proc;
};

struct hilbert_data {
  double hkey;
  int grid_num;
};

// Used in DepositParticleMassFlaggingField.C
struct two_int {
  int grid;
  int proc;
};

#endif
