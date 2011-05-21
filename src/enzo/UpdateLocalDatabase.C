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

const char creation_query[] = 
"CREATE TABLE IF NOT EXISTS enzo_outputs ("\
"dset_uuid TEXT PRIMARY KEY, "\
"pf_path TEXT NOT NULL, "\
"creation_time INTEGER NOT NULL, "\
"simulation_uuid TEXT)";
 
/* function prototypes */
 
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
  char *errmsg;
  char *Fullpath = realpath(Filename, NULL);
  char insertion_query[1024];
  snprintf(insertion_query, 1023,
           "INSERT INTO enzo_outputs VALUES ('%s',"
           "'%s', %"ISYM", '%s')",
           dset_uuid, Fullpath, CurrentTimeID, MetaData.SimulationUUID);
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
