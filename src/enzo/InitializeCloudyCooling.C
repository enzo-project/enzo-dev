/***********************************************************************
/
/  INITIALIZE CLOUDY COOLING
/
/  written by: Britton Smith
/  date:       November, 2005
/  modified1:  May, 2009
/              Converted Cloudy table format from ascii to hdf5.
/
/  PURPOSE:  Read in heating, cooling, and mean molecular weight values 
/            from file.
/  Rank = 1: interpolate over temperature.
/  Rank = 2: interpolate over density and temperature.
/  Rank = 3: interpolate over density, metallicity, and temperature.
/  Rank = 4: interpolate over density, metallicity, electron fraction, 
/            and temperature.
/  Rank = 5: interpolate over density, metallicity, electron fraction, 
/            redshift (for Haardt Madau background), and temperature.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include "hdf5.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

#define SMALL_LOG_VALUE -99.0

/**************************** Functions Prototypes ******************************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

// Initialize Cloudy Cooling
int InitializeCloudyCooling(FLOAT Time)
{

  FLOAT a = 1, dadt;
  int q, w;
  float64 *temp_data;
  long_int temp_int;
  long_int *temp_int_arr;
  char parameter_name[MAX_LINE_LENGTH];

  // Initialize things needed even if cloudy cooling is not used.

    CloudyCoolingData.CloudyCoolingGridParameters = new float*[CLOUDY_COOLING_MAX_DIMENSION];
    CloudyCoolingData.CloudyCoolingGridDimension = new int[CLOUDY_COOLING_MAX_DIMENSION];
    for (q = 0;q < CLOUDY_COOLING_MAX_DIMENSION;q++) {
      CloudyCoolingData.CloudyCoolingGridDimension[q] = 0; 
    }

  // Zero arrays if cloudy cooling not used.

  if (MetalCooling != CLOUDY_METAL_COOLING) {
    CloudyCoolingData.CloudyCoolingGridRank = 0;
    return SUCCESS;
  }

  if (debug) {
    fprintf(stderr,"Initializing Cloudy cooling.\n");
    fprintf(stderr,"CloudyCoolingGridFile: %s.\n",CloudyCoolingData.CloudyCoolingGridFile);
    fprintf(stderr,"IncludingCloudyHeating: %"ISYM".\n",CloudyCoolingData.IncludeCloudyHeating);
    fprintf(stderr,"CMBTemperatureFloor: %"ISYM".\n",CloudyCoolingData.CMBTemperatureFloor);
  }

  /* If using cosmology, compute the expansion factor and get units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  if (ComovingCoordinates) {

    if (CosmologyComputeExpansionFactor(Time, &a, &dadt) 
	== FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      return FAIL;
    }

    aUnits = 1.0/(1.0 + InitialRedshift);

  }

  /* Get conversion units. */

  double tbase1 = TimeUnits;
  double xbase1 = LengthUnits/(a*aUnits);
  double dbase1 = DensityUnits * POW(a*aUnits, 3);
  double mh = 1.67e-24;
  double CoolUnit = (POW(aUnits,5) * POW(xbase1,2) * POW(mh,2)) /
                    (POW(tbase1,3) * dbase1);

  // Read cooling data in from hdf5 file.

  hid_t       file_id, dset_id, attr_id; 
  herr_t      status;
  herr_t      h5_error = -1;

  if (debug) fprintf(stderr,"Reading Cloudy data from %s.\n", CloudyCoolingData.CloudyCoolingGridFile);
  file_id = H5Fopen(CloudyCoolingData.CloudyCoolingGridFile, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open cooling dataset and get grid dimensions.

  dset_id =  H5Dopen(file_id, "/Cooling");
  if (dset_id == h5_error) {
    fprintf(stderr,"Can't open Cooling in %s.\n",CloudyCoolingData.CloudyCoolingGridFile);
    return FAIL;
  }

  // Grid rank.
  attr_id = H5Aopen_name(dset_id, "Rank");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open Rank attribute in Cooling dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8, &temp_int);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to read Rank attribute in Cooling dataset.\n");
    return FAIL;
  }
  CloudyCoolingData.CloudyCoolingGridRank = (int) temp_int;
  if (debug) fprintf(stderr,"Cloudy cooling grid rank: %"ISYM".\n",CloudyCoolingData.CloudyCoolingGridRank);
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Rank attribute in Cooling dataset.\n");
    return FAIL;
  }

  // Grid dimension.
  temp_int_arr = new long_int[CloudyCoolingData.CloudyCoolingGridRank];
  attr_id = H5Aopen_name(dset_id, "Dimension");
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to open Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  status = H5Aread(attr_id, HDF5_I8,temp_int_arr);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to read Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  if (debug) fprintf(stderr,"Cloudy cooling grid dimensions:");
  for (q = 0;q < CloudyCoolingData.CloudyCoolingGridRank;q++) {
    CloudyCoolingData.CloudyCoolingGridDimension[q] = (int) temp_int_arr[q];
    if (debug) fprintf(stderr," %"ISYM,CloudyCoolingData.CloudyCoolingGridDimension[q]);
  }
  if (debug) fprintf(stderr,".\n");
  status = H5Aclose(attr_id);
  if (attr_id == h5_error) {
    fprintf(stderr,"Failed to close Dimension attribute in Cooling dataset.\n");
    return FAIL;
  }
  delete [] temp_int_arr;

  // Read Cooling data.
  CloudyCoolingData.CloudyDataSize = 1;
  for (q = 0;q < CloudyCoolingData.CloudyCoolingGridRank;q++) {
    CloudyCoolingData.CloudyDataSize *= CloudyCoolingData.CloudyCoolingGridDimension[q];
  }
  temp_data = new float64[CloudyCoolingData.CloudyDataSize];

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
  if (debug) fprintf(stderr,"Reading Cloudy Cooling dataset.\n");
  if (status == h5_error) {
    fprintf(stderr,"Failed to read Cooling dataset.\n");
    return FAIL;
  }

  CloudyCoolingData.CloudyCooling = new float[CloudyCoolingData.CloudyDataSize];
  for (q = 0;q < CloudyCoolingData.CloudyDataSize;q++) {
    CloudyCoolingData.CloudyCooling[q] = temp_data[q] > 0 ? (float) log10(temp_data[q]) : (float) SMALL_LOG_VALUE;

    // Convert to code units.
    CloudyCoolingData.CloudyCooling[q] -= log10(CoolUnit);
  }
  delete [] temp_data;

  status = H5Dclose(dset_id);
  if (status == h5_error) {
    fprintf(stderr,"Failed to close Cooling dataset.\n");
    return FAIL;
  }

  // Read Heating data.
  if (CloudyCoolingData.IncludeCloudyHeating > 0) {

    temp_data = new float64[CloudyCoolingData.CloudyDataSize];

    dset_id =  H5Dopen(file_id, "/Heating");
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open Heating in %s.\n",CloudyCoolingData.CloudyCoolingGridFile);
      return FAIL;
    }

    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
    if (debug) fprintf(stderr,"Reading Cloudy Heating dataset.\n");
    if (status == h5_error) {
      fprintf(stderr,"Failed to read Heating dataset.\n");
      return FAIL;
    }

    CloudyCoolingData.CloudyHeating = new float[CloudyCoolingData.CloudyDataSize];
    for (q = 0;q < CloudyCoolingData.CloudyDataSize;q++) {
      CloudyCoolingData.CloudyHeating[q] = temp_data[q] > 0 ? (float) log10(temp_data[q]) : (float) SMALL_LOG_VALUE;

      // Convert to code units.
      CloudyCoolingData.CloudyHeating[q] -= log10(CoolUnit);
    }
    delete [] temp_data;

    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close Heating dataset.\n");
      return FAIL;
    }
  }

  // Read in grid parameters.
  for (q = 0;q < CloudyCoolingData.CloudyCoolingGridRank;q++) {

    if (q < CloudyCoolingData.CloudyCoolingGridRank - 1) {
      sprintf(parameter_name,"/Parameter%"ISYM,(q+1));
    }
    else {
      sprintf(parameter_name,"/Temperature");
    }

    temp_data = new float64[CloudyCoolingData.CloudyCoolingGridDimension[q]];

    dset_id =  H5Dopen(file_id, parameter_name);
    if (dset_id == h5_error) {
      fprintf(stderr,"Can't open %s in %s.\n",parameter_name,CloudyCoolingData.CloudyCoolingGridFile);
      return FAIL;
    }

    status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
    if (debug) fprintf(stderr,"Reading Cloudy %s dataset.\n",parameter_name);
    if (status == h5_error) {
      fprintf(stderr,"Failed to read %s dataset.\n",parameter_name);
      return FAIL;
    }

    CloudyCoolingData.CloudyCoolingGridParameters[q] = new float[CloudyCoolingData.CloudyCoolingGridDimension[q]];
    for (w = 0;w < CloudyCoolingData.CloudyCoolingGridDimension[q];w++) {
      if (q < CloudyCoolingData.CloudyCoolingGridRank - 1) {
	CloudyCoolingData.CloudyCoolingGridParameters[q][w] = (float) temp_data[w];
      }
      else {
	// convert temeperature to log
	CloudyCoolingData.CloudyCoolingGridParameters[q][w] = (float) log10(temp_data[w]);
      }
    }
    delete [] temp_data;

    status = H5Dclose(dset_id);
    if (status == h5_error) {
      fprintf(stderr,"Failed to close %s dataset.\n",parameter_name);
      return FAIL;
    }

    if (debug) fprintf(stderr,"%s: %"GSYM" to %"GSYM" (%"ISYM" steps).\n",parameter_name,
		       CloudyCoolingData.CloudyCoolingGridParameters[q][0],
		       CloudyCoolingData.CloudyCoolingGridParameters[q][CloudyCoolingData.CloudyCoolingGridDimension[q]-1],
		       CloudyCoolingData.CloudyCoolingGridDimension[q]);

  } // for (q = 0;q < CloudyCoolingData.CloudyCoolingGridRank;q++)

  status = H5Fclose (file_id);

  if (CloudyCoolingData.CloudyCoolingGridRank > CLOUDY_COOLING_MAX_DIMENSION) {
    fprintf(stderr,"Error: rank of Cloudy cooling data must be less than or equal to %"ISYM".\n",
	    CLOUDY_COOLING_MAX_DIMENSION);
    return FAIL;
  }

  return SUCCESS;
}
