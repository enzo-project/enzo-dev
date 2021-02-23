/****************************************************************
/
/ STELLAR YIELDS DATA STRUCTURE AND FUNCTIONS
/
/ written by: Andrew Emerick
/ date:       March, 2016
/ modified1:
/ date     :
/
/ PURPOSE:
/
****************************************************************/

// may not need any includes


// i_m = mass, i_z = metallicity, i_y = yield

#define YIELD_INDEX(i_m, i_z, i_y, Nm, Nz) (i_m + (i_z + i_y*Nz)*Nm)

#define YIELD_INDEX_DM(Nm,Nz) (1)
#define YIELD_INDEX_DZ(Nm,Nz) (Nm)
#define YIELD_INDEX_DY(Nm,Nz) (Nm*Nz)

struct StellarYieldsDataType
{
  int Nm; // Number of mass bins
  int Nz; // Number of metallicity bins
  int Ny; // Number of Yields

  int dm; // Next mass   (fixed Z, yield)
  int dz; // Next metallicity (fixed M, yield)
  int dy; // Next yield (field M, Z)

  int size;

  float *M;
  float *Z;

  // arrays are 1D representations of multi-D arrays
  float  *Mtot;       // 2D  (NmxNz)
  float  *Metal_Mtot; // 2D  (NmxNz)
  float  *Yields;     // 3D  (NmxNzxNy)

  // indexes are:
  //
  //    2D (NmxNz)    :   i + j*Nm
  //    3D (NmxNzxNy) :   i + (j + k*Nz)*Nm
  //
  //    where i iterates over mass, j over Z, k over yield
};

struct MetalMixingExperimentDataType
{
  int NumberOfEvents;

  float *time;

  float *xpos;
  float *ypos;
  float *zpos;

  float *M_ej;
  float *E_ej;

  int *anums;
  float **yield;     // MAX_STELLAR_ABUND length for each action
};
