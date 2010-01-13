/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS FOR ANALYSIS ROUTINES
/
************************************************************************/

#define MAX_PROFILES 200   //############# I found it finally!!!

struct AnalyzeClusterParameters {

  /* Radial profile parameters (in Mpc). */

  FLOAT rinner; 
  FLOAT router;
  int npoints;

  /* Virial overdensity parameter (times critical density). */

  float virial_dens; 

  /* Fraction of rvir within which the mean velocity should be computed. */

  float MeanVelocityVirialFraction;

  /* The maximum temperature of "cold gas" in K.  If the second parameter
     is > 0, then the first is set to Tvir * fraction. */

  float ColdTemperatureCutoff; 
  float ColdTemperatureCutoffVirialFraction;

  /* Normalization of the M-T-z relation (f_T in Bryan & Norman 1998). */

  float VirialTemperatureNormalization;

  /* Lower cut-off used in computing dense gas and cold gas fractions,
     as well as L_dense (in solar masses/Mpc^3). */

  float LowerDensityCutoff;

  /* Upper density cutoff used in computing L_dens (Msolar/Mpc^3). */

  float UpperDensityCutoff;

  /* Boolean: if set to TRUE (1), then computes disk information. */

  int   ComputeDiskInformation; 

  /* Disk image size in pixels, and size of disk image in terms of the virial
     radius. */

  int   DiskImageSize; 
  float DiskRadius;

  /* Lower and Upper energy cutoff in luminosity look-up table computation,
     in keV.  The last variable is the file name of the look-up table.
     (if NULL, then rho^2 * T is used instead). */

  float XrayLowerCutoffkeV;
  float XrayUpperCutoffkeV;
  char  *XrayTableFileName;

  /* Flag indicating if the special clumping factor computation should be
     called. */

  int ComputeClumpingFactor;

  /* Any extra information to be placed in the AnalyzeCluster file */

  char *MetaData;

  /* For Vertical profile */

  float DiskRadiusCutoff;
  int LinearProfileRadiusForVertical;

  /* For printing out global profile values (averaged gas surface densities 
     on the disk, averaged SFR surface densities on the disk, etc.) */ 
  
  int PrintGlobalProfileValues;

};
