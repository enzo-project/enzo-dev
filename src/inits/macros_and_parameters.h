/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

/* Modifiable Parameters */

#define MAX_LINE_LENGTH          512

#define VERSION                  2.0  /* current version number */

#define MAX_SPECIES               10

#define MAX_PS_NAMES               2

#define MAX_LEVELS                10 /* Used by AutomaticSubgridGeneration */

#define NO_SHIFT_FOR_LARS

#define MAX_DIMENSION              3  /* don't change this */

#define PARTICLE_TYPE_GAS          0
#define PARTICLE_TYPE_DARK_MATTER  1
#define PARTICLE_TYPE_STAR         2
#define PARTICLE_TYPE_TRACER       3
#define PARTICLE_TYPE_MUST_REFINE  4

/* Fortran name generator */

#if defined(SUN_OLD)
#define FORTRAN_NAME(NAME) NAME/**/_
#endif

#if defined(IRIS4) || defined(CONVEX) || defined(COMPAQ) || defined(SUN) || defined(LINUX) || defined(IA64) || defined(CRAYX1) || defined(XT3)
#define FORTRAN_NAME(NAME) NAME##_
#endif

#if defined(SPP) || defined(SP2)
#define FORTRAN_NAME(NAME) NAME
#endif

/* Precision-related definitions. */

/* Define the long int type. */

typedef long long long_int;

/* Define the long double type. */

typedef long double long_double;

/* Macro definitions for portability */

typedef unsigned int   Eunsigned_int;
typedef long long int  Elong_int;
typedef void           *VOIDP;
typedef int            Eint32;
typedef long long int  Eint64;
typedef float          Eflt32;
typedef double         Eflt64;
typedef long double    Eflt128;

typedef int            MPI_Arg;

/* HDF5 definitions */


#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I8 H5T_STD_I64BE
#define HDF5_FILE_R4 H5T_IEEE_F32BE
#define HDF5_FILE_R8 H5T_IEEE_F64BE
#define HDF5_FILE_B8 H5T_STD_B8BE


/* #ifdef HDF5_LE */
/* #define HDF5_FILE_I4 H5T_STD_I32LE */
/* #define HDF5_FILE_I8 H5T_STD_I64LE */
/* #define HDF5_FILE_R4 H5T_IEEE_F32LE */
/* #define HDF5_FILE_R8 H5T_IEEE_F64LE */
/* #define HDF5_FILE_B8 H5T_STD_B8LE */
/* #endif */

#define HDF5_I4 H5T_NATIVE_INT
#define HDF5_I8 H5T_NATIVE_LLONG
#define HDF5_R4 H5T_NATIVE_FLOAT
#define HDF5_R8 H5T_NATIVE_DOUBLE
#define HDF5_R16 H5T_NATIVE_LDOUBLE

/* Precision-dependent definitions */

#ifdef SMALL_INTS
#define Eint int
#define ISYM "d"
#define HDF5_INT HDF5_I4
#define HDF5_FILE_INT HDF5_FILE_I4
#define nint(A) ( (int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) abs((int) (A))
#endif

#ifdef LARGE_INTS
#define Eint long_int
#define int long_int
#define ISYM "lld"
#define HDF5_INT HDF5_I8
#define HDF5_FILE_INT HDF5_FILE_I8
#define nint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) labs((long_int) (A))
#endif

#ifdef CONFIG_BFLOAT_4
#define Eflt float
#define FLOAT float
#define FSYM "f"
#define GSYM "g"
#define GOUTSYM ".7g"
#define HDF5_REAL HDF5_R4
#define HDF5_FILE_REAL HDF5_FILE_R4
#endif /* r4 */

#ifdef CONFIG_BFLOAT_8
#define Eflt double
#define FLOAT double
#define FSYM "lf"
#define GSYM "g"
#define GOUTSYM ".14g"
#define HDF5_REAL HDF5_R8
#define HDF5_FILE_REAL HDF5_FILE_R8
#endif /* r8 */


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

#define HDF_FAIL -1

#define FLOAT_UNDEFINED  -99999.0
#define INT_UNDEFINED    -99999

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define POW(X,Y) pow((double) (X), (double) (Y))

/* Definitions for FastFourierTransform and related routines */

#define FFT_FORWARD  +1
#define FFT_INVERSE  -1
#define REAL_TO_COMPLEX    0
#define COMPLEX_TO_COMPLEX 1

/* RH debug */

#define RH_D TRUE

/* RH memory trace */

#define RH_M TRUE

/* RH flow trace */

#define RH_F TRUE

#ifndef OLD_HDF5
#  define hssize_t hsize_t
#endif

