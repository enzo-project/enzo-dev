/***********************************************************************
/
/  WRITES THE NEW OUTPUT INTO THE DATABASE FILE
/
/  written by: Matthew Turk
/  date:       May, 2011
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine writes the parameter file in the argument and sets parameters
//   based on it.
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifdef USE_SQLITE
#include "sqlite3.h"
#endif
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"

const char creation_query[] = 
"CREATE TABLE IF NOT EXISTS simulation_outputs ("\
"dset_uuid TEXT PRIMARY KEY, "\
"output_type TEXT NOT NULL, "\
"pf_path TEXT NOT NULL, "\
"creation_time INTEGER NOT NULL, "\
"last_seen_time INTEGER NOT NULL, "\
"simulation_uuid TEXT NOT NULL, "\
"redshift REAL, time REAL, "\
"topgrid0 INTEGER, topgrid1 INTEGER, topgrid2 INTEGER)";
 
/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int  GetUnits(float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MAssUnits, FLOAT Time);


int UpdateLocalDatabase(TopGridData &MetaData, int CurrentTimeID,
                        char *dset_uuid, char *Filename)
{
 
#ifdef USE_SQLITE
  if((DatabaseLocation == NULL) ||
     (MyProcessorNumber != ROOT_PROCESSOR)) return SUCCESS;

  int retval = 0;

  static sqlite3 *handle = NULL;
  static int DatabaseInitialized = 0;
  if (DatabaseInitialized == 0) {
    retval = sqlite3_open(DatabaseLocation, &handle);
    if(retval){
      fprintf(stderr, "Database connection failed to: %s\n", DatabaseLocation);
      /* Seems kind of silly to kill the simulation over such a trivial thing ...
       */
      DatabaseInitialized = -1;
      return SUCCESS;
    }

    retval = sqlite3_exec(handle, creation_query,0,0,0);
    if(retval){
      fprintf(stderr, "Failed executing table query.");
      DatabaseInitialized = -1;
      return SUCCESS;
    }
  } else if (DatabaseInitialized == -1) {
    return SUCCESS;
  }
  /* Compute some quantities. */
  
  /* Compute the current redshift if we're using cosmology. */
  FLOAT a, dadt, CurrentRedshift;
  if (ComovingCoordinates) {
      CosmologyComputeExpansionFactor(MetaData.Time, &a, &dadt);
      CurrentRedshift = (1 + InitialRedshift)/a - 1;
  } else {
      CurrentRedshift = -1.;
  }
  
  /* Get the current simulation time in seconds. */
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits,  MetaData.Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
  float this_time = MetaData.Time * TimeUnits;
  
  /* Top Grid Dimensions. */
  int topgrid[3], dim;
  for (dim = 0; dim < 3; dim++) topgrid[dim] = 0;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      topgrid[dim] = MetaData.TopGridDims[dim];
  }
  
  char *errmsg;
  char *Fullpath = realpath(Filename, NULL);
  char insertion_query[1024];
  snprintf(insertion_query, 1023,
           "INSERT INTO simulation_outputs VALUES ('%s',"
           "'EnzoStaticOutput',"
           "'%s', %"ISYM", %"ISYM", '%s', %"GOUTSYM", %"ESYM","
           "%"ISYM", %"ISYM", %"ISYM")",
           dset_uuid, Fullpath, CurrentTimeID, CurrentTimeID,
           MetaData.SimulationUUID,
           CurrentRedshift, this_time,
           topgrid[0], topgrid[1], topgrid[2]);
  free(Fullpath);
  retval = sqlite3_exec(handle, insertion_query, 0,0,&errmsg);
  if(retval){
    fprintf(stderr, "Failed to insert the current output.\n");
    fprintf(stderr, "ERR: %s\n", errmsg);
    return SUCCESS;
  }
 
#endif
  return SUCCESS;
}
