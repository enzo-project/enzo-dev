#ifndef __macros_and_parameters_h_
#define __macros_and_parameters_h_
/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/
#ifdef USE_PYTHON
#ifndef ENZO_PYTHON_IMPORTED
#define PY_ARRAY_UNIQUE_SYMBOL enzo_ARRAY_API
#define NO_IMPORT_ARRAY 
#include <Python.h>
#include "numpy/arrayobject.h"
#endif
#endif

#ifdef ECUDA
#ifdef LARGE_INTS
// CUDA hates LARGE_INTS, and who can blame it?
#error "Sorry, you need to be using 32 bit integers with CUDA because of #define int!"
#endif // LARGE_INTS
#ifdef CONFIG_BFLOAT_8
#error "Sorry, you need to be using 32 bit precision with CUDA because of #define float!"
#endif
#endif

/* Modifiable Parameters */

#define MAX_NUMBER_OF_TASKS             16384

#define MAX_NUMBER_OF_NODES              1024

#define MAX_TASKS_PER_NODE                   __max_cpu_per_node

#define MAX_NUMBER_OF_BARYON_FIELDS          __max_baryons  /* must be at least 6 */

#define MAX_NUMBER_OF_SUBGRIDS               __max_subgrids

#define MAX_DEPTH_OF_HIERARCHY             50

#define MAX_LINE_LENGTH                   512

#define MAX_NAME_LENGTH                   512

#define MAX_GRID_TAG_SIZE         16
#define MAX_TASK_TAG_SIZE         16
#define MAX_GROUP_TAG_SIZE        16
#define MAX_CYCLE_TAG_SIZE        16
#define GRID_TAG_FORMAT        "4.4"
#define TASK_TAG_FORMAT        "4.4"
#define GROUP_TAG_FORMAT       "8.8"
#define CYCLE_TAG_FORMAT       "4.4"
#define MAX_COUNTERS              40

#define MEMORY_POOL_SIZE  __memory_pool_size

#define DEFAULT_GHOST_ZONES                 3  /* at least 3 */

#define MAX_NUMBER_OF_OUTPUT_REDSHIFTS    500

#define GRAVITY_BUFFER_SIZE                 3

#define MAX_FLAGGING_METHODS                9

#define MAX_STATIC_REGIONS               1000

#define MAX_REFINE_REGIONS               150

#ifdef WINDS 
#define MAX_NUMBER_OF_PARTICLE_ATTRIBUTES  7
#else
#define MAX_NUMBER_OF_PARTICLE_ATTRIBUTES  4
#endif

#define MAX_TIME_ACTIONS                   10

#define MAX_CUBE_DUMPS                     50

#define MAX_CUBE_DUMPS                     50

#define MAX_MOVIE_FIELDS                    6

#define MAX_POTENTIAL_ITERATIONS            80

#define MAX_ENERGY_BINS                    10

#define ROOT_PROCESSOR                      0

#define VERSION                             2.0  /* current version number */

/* Unmodifiable Parameters */

#define MAX_DIMENSION                       3  /* must be 3! */

#include "message.h"

/* We have two possibilities for throwing an exception.
   You can throw ENZO_FAIL, which takes no arguments and is safe to use with a
   semicolon, or you can throw ENZO_VFAIL which must NOT have a semicolon and
   comes enclosed in brackets.  You can supply format strings to ENZO_VFAIL.  */
#ifdef CONFIG_THROW_ABORT
#define ENZO_FAIL(A) raise(SIGABRT);
#define ENZO_VFAIL(A, ...) raise(SIGABRT);
#else
#define ENZO_FAIL(A) throw(EnzoFatalException(A, __FILE__, __LINE__));
#define ENZO_VFAIL(format, ...) {snprintf(current_error, 254, format, ##__VA_ARGS__); throw(EnzoFatalException(current_error, __FILE__, __LINE__));}
#endif

/* Fortran name generator (cpp blues) */

#if defined(SUN_OLD)
#define FORTRAN_NAME(NAME) NAME/**/_
#endif

#if defined(IRIS4) || defined(CONVEX) || defined(COMPAQ) || defined(SUN) || defined(LINUX) || defined(IA64) || defined(CRAYX1) || defined(XT3)
#define FORTRAN_NAME(NAME) NAME##_
#endif

#if defined(SPP) || defined(SP2) || defined(BGL)
#define FORTRAN_NAME(NAME) NAME
#endif

#ifdef CONFIG_PFLOAT_16
#define PFORTRAN_NAME(NAME) NAME##_c
#else
#define PFORTRAN_NAME(NAME) FORTRAN_NAME(NAME)
#endif

/* Precision-related definitions. */

typedef long long long_int;
typedef long double long_double;
typedef unsigned int unsigned_int;
typedef unsigned long long int unsigned_long_int;

/* Previously in hdf4.h */

typedef float        float32;
typedef double       float64;
typedef long double  float128;

/* Macro definitions for portability */

typedef void           *VOIDP;
typedef int            Eint32;
typedef long long int  Eint64;
typedef float          Eflt32;
typedef double         Eflt64;
typedef long double    Eflt128;
typedef long long int  Elong_int;

typedef int            MPI_Arg;

typedef int            HDF5_hid_t;

/* HDF5 definitions */

#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I8 H5T_STD_I64BE
#define HDF5_FILE_R4 H5T_IEEE_F32BE
#define HDF5_FILE_R8 H5T_IEEE_F64BE
#define HDF5_FILE_B8 H5T_STD_B8BE

#define HDF5_I4  H5T_NATIVE_INT
#define HDF5_I8  H5T_NATIVE_LLONG
#define HDF5_R4  H5T_NATIVE_FLOAT
#define HDF5_R8  H5T_NATIVE_DOUBLE
#define HDF5_R16 H5T_NATIVE_LDOUBLE

/* Precision-dependent definitions */

#if defined(INITS32)
#define inits_type float32
#endif

#if defined(INITS64)
#define inits_type float64
#endif

#define ByteDataType MPI_BYTE

#ifdef SMALL_INTS
#define Eint int
#define Eunsigned_int unsigned_int
#define ISYM "d"
#define IntDataType MPI_INT
#define HDF5_INT HDF5_I4
#define HDF5_FILE_INT HDF5_FILE_I4
#define nint(A) ( (int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) abs((int) (A))
#define ENPY_INT NPY_INT
#define enpy_int npy_int
#endif

#ifdef LARGE_INTS
#define int long_int // CUDA doesn't like this, and who can blame it?
#define Eint long_int
#define Eunsigned_int unsigned_long_int
#define ISYM "lld"
#define IntDataType MPI_LONG_LONG_INT
#define HDF5_INT HDF5_I8
#define HDF5_FILE_INT HDF5_FILE_I8
#define nint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) labs((long_int) (A))
#define ENPY_INT NPY_LONG
#define enpy_int npy_long
#endif

#ifdef CONFIG_BFLOAT_4
#define BFLOAT_EPSILON 1e-6f
#define Eflt float
#define FSYM "f"
#define ESYM "e"
#define FloatDataType MPI_FLOAT
#ifdef COMPACT_IO
#define HDF5_REAL HDF5_R4
#define HDF5_FILE_REAL HDF5_FILE_R4
#else
#define HDF5_REAL HDF5_R4
#define HDF5_FILE_REAL HDF5_FILE_R8
#endif
#ifdef USE_PYTHON
#define ENPY_BFLOAT NPY_FLOAT
#define enpy_bfloat npy_float
#endif
#endif

#ifdef CONFIG_BFLOAT_8
#define BFLOAT_EPSILON 1e-12f
#define Eflt double
#define FSYM "lf"
#define ESYM "le"
#define FloatDataType MPI_DOUBLE
#define float32 TEMP_HOLD_NAME
#define float double
#define TEMP_HOLD_NAME float32
#define HDF5_REAL HDF5_R8
#define HDF5_FILE_REAL HDF5_FILE_R8
#ifdef USE_PYTHON
#define ENPY_BFLOAT NPY_DOUBLE
#define enpy_bfloat npy_double
#endif
#endif

#ifdef CONFIG_PFLOAT_4
#define PFLOAT_EPSILON 1e-6f
#define FLOAT Eflt32
#define PEXP expf
#define PSYM "f"
#define GSYM "g"
#define GOUTSYM ".8g"
#define MY_MPIFLOAT MPI_FLOAT
#define FLOATDataType MPI_FLOAT
#define HDF5_PREC HDF5_R4
#define HDF5_FILE_PREC HDF5_R4
#ifdef USE_PYTHON
#define ENPY_PFLOAT NPY_FLOAT
#define enpy_pfloat npy_float
#endif
#endif

#ifdef CONFIG_PFLOAT_8
#define PFLOAT_EPSILON 1e-12f
#define FLOAT double
#define PEXP exp
#define PSYM "lf"
#define GSYM "g"
#define GOUTSYM ".14g"
#define MY_MPIFLOAT MPI_DOUBLE
#define FLOATDataType MPI_DOUBLE
#define HDF5_PREC HDF5_R8
#define HDF5_FILE_PREC HDF5_R8
#ifdef USE_PYTHON
#define ENPY_PFLOAT NPY_DOUBLE
#define enpy_pfloat npy_double
#endif
#endif

#ifdef CONFIG_PFLOAT_16
#define PFLOAT_EPSILON 1e-16f
#define FLOAT long_double
#define PEXP expl
#define PSYM "Lf"
#define GSYM "g"
#define GOUTSYM ".21Lg"
#define MY_MPIFLOAT MPI_LONG_DOUBLE
#define FLOATDataType MPI_LONG_DOUBLE
#define HDF5_PREC HDF5_R16
#define HDF5_FILE_PREC HDF5_R16
#ifdef USE_PYTHON
#define ENPY_PFLOAT NPY_LONGDOUBLE
#define enpy_pfloat npy_longdouble
#endif
#endif

/* Definitions for controlling the integer type for particle IDs
   (8-byte needed for >2 billion particle simulations) */

#ifdef CONFIG_PINT_4
#define PINT Eint32
#define PINTDataType MPI_INT
#define HDF5_PINT HDF5_I4
#define HDF5_FILE_PINT HDF5_FILE_I4
#define PISYM "d"
#define ENPY_PINT NPY_INT
#endif

#ifdef CONFIG_PINT_8
#define PINT Eint64
#define PINTDataType MPI_LONG_LONG_INT
#define HDF5_PINT HDF5_I8
#define HDF5_FILE_PINT HDF5_FILE_I8
#define PISYM "lld"
#define ENPY_PINT NPY_LONG
#endif

/* Standard definitions (well, fairly standard) */

#ifndef NULL
#define NULL      0
#endif

#ifdef FAIL
#undef FAIL
#endif
#define FAIL      0
#define SUCCESS   1

#ifndef FALSE
#define FALSE     0
#define TRUE      1
#endif

/* Not-so standard definitions */
#ifndef HDF_FAIL
#define HDF_FAIL -1
#endif

#define FLOAT_UNDEFINED  -99999.0
#define INT_UNDEFINED    -99999

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define POW(X,Y) pow((double) (X), (double) (Y))
#define COS(X) cos((double) (X))
#define SIN(X) sin((double) (X))

/* Macros for grid indices (with and without ghost zones, and
   vertex-centered data) */

#define GRIDINDEX(A,B,C) ((((C)+GridStartIndex[2])*GridDimension[1]+((B)+GridStartIndex[1]))*GridDimension[0]+GridStartIndex[0]+(A))
#define GRIDINDEX_NOGHOST(A,B,C) (((C)*GridDimension[1]+(B))*GridDimension[0]+A)
#define VCGRIDINDEX(A,B,C) (((C)*(GridEndIndex[1]-GridStartIndex[1]+2) + (B)) * (GridEndIndex[0]-GridStartIndex[0]+2) + A)

/* Definitions for FastFourierTransform and related routines */

#define FFT_FORWARD  +1
#define FFT_INVERSE  -1
#define REAL_TO_COMPLEX    0
#define COMPLEX_TO_COMPLEX 1

/* Definitions for grid::RestoreEnergyConsistency */

#define ENTIRE_REGION  0
#define ONLY_BOUNDARY  1

/* Definitions for grid::ZeroSolutionUnderSubgrid */

#define ZERO_ALL_FIELDS          0
#define ZERO_UNDER_SUBGRID_FIELD 1

/* Definitions for grid::CommunicationSend/ReceiveRegion and 
   grid::DepositPositions */

#define INTERPOLATED_FIELDS              -8
#define PARTICLE_MASS_FLAGGING_FIELD     -7
#define MASS_FLAGGING_FIELD              -6
#define ACCELERATION_FIELDS              -5
#define POTENTIAL_FIELD                  -4
#define GRAVITATING_MASS_FIELD           -3
#define GRAVITATING_MASS_FIELD_PARTICLES -2
#define ALL_FIELDS   -1

#define NEW_AND_OLD   0
#define NEW_ONLY      1
#define OLD_ONLY      2

/* Definitions for grid::ComputeAccelerationField */

#define PARTICLES  0
#define GRIDS      1
#define ZEUS_GRIDS 2

#define DIFFERENCE_TYPE_STAGGERED 0
#define DIFFERENCE_TYPE_NORMAL    1

/* Definitions for CommunicationTranspose */

#define NORMAL_ORDER      0
#define TRANSPOSE_FORWARD 1
#define TRANSPOSE_REVERSE 2

/* Definitions for CommunicationTransferParticles */

#define COPY_IN   0
#define COPY_OUT  1

/* Definitions for CommunicationCollectParticles */

#define SUBGRIDS_LOCAL 0
#define SUBGRIDS_GLOBAL 1
#define SIBLINGS_ONLY 2
#define ALLGRIDS 3

/* Definitions for CommunicationDirection */

#define COMMUNICATION_SEND_RECEIVE 0
#define COMMUNICATION_POST_RECEIVE 1
#define COMMUNICATION_SEND         2
#define COMMUNICATION_RECEIVE      3

/* MPI Tags */

#define MPI_TRANSPOSE_TAG 10
#define MPI_SENDREGION_TAG 11
#define MPI_FLUX_TAG 12
#define MPI_TRANSFERPARTICLE_TAG 13
#define MPI_SENDPARTFLOAT_TAG 14
#define MPI_SENDPARTINT_TAG 15
#define MPI_PHOTON_TAG 16
#define MPI_PHOTONGROUP_TAG 17
#define MPI_DENSEST_TAG 18
#define MPI_KEEPTRANSPORTING_TAG 19
#define MPI_NPHOTON_TAG 20
#define MPI_SENDSTAR_TAG 21
#define MPI_SENDPMFLAG_TAG 22
#define MPI_SENDPART_TAG 23
#define MPI_SENDMARKER_TAG 24
#define MPI_SGMARKER_TAG 25

// There are 5 tags related to this (1000-1004)
#define MPI_SENDPARTFIELD_TAG 1000

/* Definitions for CommunicationBufferedSend. */

#define BUFFER_IN_PLACE -1

/* Particle types (note: gas is a conceptual type) */

#define NUM_PARTICLE_TYPES 11

#define PARTICLE_TYPE_RESET         -1
#define PARTICLE_TYPE_GAS            0
#define PARTICLE_TYPE_DARK_MATTER    1
#define PARTICLE_TYPE_STAR           2
#define PARTICLE_TYPE_TRACER         3
#define PARTICLE_TYPE_MUST_REFINE    4
#define PARTICLE_TYPE_SINGLE_STAR    5
#define PARTICLE_TYPE_BLACK_HOLE     6
#define PARTICLE_TYPE_CLUSTER        7
#define PARTICLE_TYPE_MBH            8
#define PARTICLE_TYPE_COLOR_STAR     9
#define PARTICLE_TYPE_SIMPLE_SOURCE 10

/* Star particle handling */

#define NORMAL_STAR	0
#define UNIGRID_STAR    1
#define KRAVTSOV_STAR   2
#define POP3_STAR	3
#define SINK_PARTICLE	4
#define STAR_CLUSTER    5
#define INSTANT_STAR    7
#define SPRINGEL_HERNQUIST_STAR 8
#define MBH_PARTICLE    9
#define COLORED_POP3_STAR  10
#define H2REG_STAR     11
#define UNIGRID_STAR_MOM 12


#define STARMAKE_METHOD(A) (StarParticleCreation >> (A) & 1)
#define STARFEED_METHOD(A) (StarParticleFeedback >> (A) & 1)

/* Feedback modes */

#define TO_DELETE -99999  // TO_DELETE is for "type", not for "FeedbackFlag"
#define NO_FEEDBACK 0
#define ACCRETION 1
#define SUPERNOVA 2
#define CONT_SUPERNOVA 3
#define FORMATION 4
#define STROEMGREN 5
#define DEATH 6
#define MBH_THERMAL 7
#define MBH_JETS 8
#define COLOR_FIELD 9

/* Sink particle accretion modes */

#define LOCAL_ACCRETION 0  // accrete a fraction of the mass in the cell
#define BONDI_ACCRETION 1
#define RADIAL_ACCRETION 2 // accrete based on a radial profile of dM/dt

/* Correcting Bondi accretion rate */

#define BONDI_ACCRETION_CORRECT_ANALYTIC -1   // by using R^-1.5 profile
#define BONDI_ACCRETION_CORRECT_NUMERICAL -2  // by stepping outwards

/* Type of metal cooling */

#define JHW_METAL_COOLING 1
#define CEN_METAL_COOLING 2
#define CLOUDY_METAL_COOLING 3

/* Definitions for grid::ComputeHeat, grid::ComputeConductionTimeStep, and grid::ConductHeat */

#define ELT(a,b,c) ( (a) + GridDimension[0]*((b) + GridDimension[1]*(c)) )

/* Streaming format parameters */

#define ALL_PARTICLES 1
#define NON_DM_PARTICLES 2
#define NON_DM_PARTICLES_MERGED_LEVEL 3
#define NON_DM_PARTICLES_MERGED_ALL 4
#define TEMPERATURE_FIELD 1000

/* Maximum number of leafs per parent in radiation source tree. */

#define MAX_LEAF 2

/* Number of entries in the Pop III IMF lookup table */

#define IMF_TABLE_ENTRIES 1000

#ifdef USE_MPI
#define MPI_INSTRUMENTATION
#else /* USE_MPI */
//#undef MEM_TRACE
#undef MPI_TRACE
#endif

//#ifndef OLD_HDF5 /* prior to HDF5-1.6.5 */
//#define hssize_t hsize_t
//#endif

#ifdef TIME_MESSAGING
#define PROCS_PER_NODE 8
#define TIME_MSG(A) \
  time (&rawtime); \
  timeinfo = localtime(&rawtime); \
  if (MyProcessorNumber % (NumberOfProcessors/PROCS_PER_NODE) == 0) \
    printf("%s :: P%d :: %s", A, MyProcessorNumber, asctime(timeinfo));
#else
#define TIME_MSG(A) ;
#endif

#endif
