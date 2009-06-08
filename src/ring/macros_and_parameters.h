/***********************************************************************
/
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

/* Precision-related definitions. */

typedef long long long_int;
typedef long double long_double;
typedef unsigned int unsigned_int;
typedef unsigned long long int unsigned_long_int;

/* Macro definitions for portability */

typedef void           *VOIDP;
typedef int            Eint32;
typedef long long int  Eint64;
typedef float          Eflt32;
typedef double         Eflt64;
typedef long double    Eflt128;
typedef long long int  Elong_int;

typedef int            MPI_Arg;


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
/* #endif */ /* HDF5_LE */

#define HDF5_I4 H5T_NATIVE_INT
#define HDF5_I8 H5T_NATIVE_LLONG
#define HDF5_R4 H5T_NATIVE_FLOAT
#define HDF5_R8 H5T_NATIVE_DOUBLE
#define HDF5_R16 H5T_NATIVE_LDOUBLE

#ifdef SMALL_INTS
#define Eint int
#define Eunsigned_int unsigned_int
#define ISYM "d"
#define IntDataType MPI_INT
#define HDF5_INT HDF5_I4
#define HDF5_FILE_INT HDF5_FILE_I4
#endif

#ifdef LARGE_INTS
#define int long_int
#define Eint long_int
#define Eunsigned_int unsigned_long_int
#define ISYM "lld"
#define IntDataType MPI_LONG_LONG_INT
#define HDF5_INT HDF5_I8
#define HDF5_FILE_INT HDF5_FILE_I8
#endif

#define Eflt double
#define FSYM "lf"
#define PSYM "lf"
#define GSYM "g"

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define nint(A) int((A) + 0.5*sign(A))
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define POW(X,Y) pow((double) (X), (double) (Y))

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

#ifndef OLD_HDF5 /* prior to HDF5-1.6.5 */
#define hssize_t hsize_t
#endif
