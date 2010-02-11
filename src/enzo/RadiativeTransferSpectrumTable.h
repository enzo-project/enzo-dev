/***********************************************************************
/
/  SPECTRUM TABLE FOR RADIATIVE TRANSFER
/
***********************************************************************/

struct RadiativeTransferSpectrumTableType
{

  int NumberOfColumnDensityBins;

  float *columndensity_table;
  float *fractionphotons_table[3];
  float *meanenergy_table;

};
