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

struct StellarYieldsDataType
{
  int NumberOfMassBins;
  int NumberOfMetallicityBins;
  int NumberOfYields;

  float *M;
  float *Z;

  int *atomic_number;

  float  **Mtot;
  float  **Metal_Mtot;
  float ***Yields;
};


