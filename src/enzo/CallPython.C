/***********************************************************************
/
/  INITIALIZE PYTHON INTERFACE AND START INTERPRETER
/
/  written by: Matthew Turk
/  date:       September, 2008
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "TopGridData.h"

#ifdef USE_PYTHON
#ifndef ENZO_PYTHON_IMPORTED
#define PY_ARRAY_UNIQUE_SYMBOL enzo_ARRAY_API
#include <Python.h>
#include "numpy/arrayobject.h"
#define ENZO_PYTHON_IMPORTED
#endif
#endif

int ExposeDataHierarchy(TopGridData *MetaData, HierarchyEntry *Grid, 
		       int &GridID, FLOAT WriteTime, int reset, int ParentID, int level);
void ExposeGridHierarchy(int NumberOfGrids);

void ExportParameterFile(TopGridData *MetaData, FLOAT CurrentTime);

int CallPython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level)
{
#ifndef USE_PYTHON
    return SUCCESS;
#else
    if (LevelArray[level+1] != NULL) return SUCCESS;
    FLOAT CurrentTime;
    int num_grids, start_index;
    num_grids = start_index = 0;

    PyDict_Clear(grid_dictionary);
    PyDict_Clear(old_grid_dictionary);

    LevelHierarchyEntry *Temp2 = LevelArray[0];
    /* Count the grids */
    /* I think there is a better idiom for this somewhere
       but I couldn't find it, and I think this works     */
    for (int lc = 0; LevelArray[lc] != NULL; lc++){
        Temp2 = LevelArray[lc];
        while (Temp2 != NULL) {
            num_grids++; Temp2 = Temp2->NextGridThisLevel;
        }
    }
    ExposeGridHierarchy(num_grids);
    Temp2 = LevelArray[0];
    while (Temp2->NextGridThisLevel != NULL)
        Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
    CurrentTime = LevelArray[level]->GridData->ReturnTime();
    if (ExposeDataHierarchy(MetaData, Temp2->GridHierarchyEntry, start_index,
                CurrentTime, 1, 0, 0) == FAIL) {
        fprintf(stderr, "Error in ExposeDataHierarchy\n");
        return FAIL;
    }

  ExportParameterFile(MetaData, CurrentTime);

  PyRun_SimpleString("import user_script\nuser_script.main()\n");

  NumberOfPythonCalls++;

  PyDict_Clear(grid_dictionary);
  PyDict_Clear(hierarchy_information);
  PyDict_Clear(yt_parameter_file);
  PyDict_Clear(conversion_factors);
  return SUCCESS;
#endif
}

