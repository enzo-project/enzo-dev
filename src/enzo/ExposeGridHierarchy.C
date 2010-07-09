/***********************************************************************
/
/  GENERATE GRID HIERARCHY HANDLERS
/
/  written by: Matthew Turk
/  date:       September, 2008
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

// This function writes out the data hierarchy (TopGrid)

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

void ExposeGridHierarchy(int NumberOfGrids)
{

#ifdef USE_PYTHON
#undef int

  /* This function just fills the dictionaries */

  // We'll set up some common variables
  npy_intp flat_dimensions[2];
  PyArrayObject *temp_array;

  PyDict_Clear(hierarchy_information);
  flat_dimensions[0] = (npy_intp) NumberOfGrids;

  /* Here are the arrays yt expects from the data file:
    GridDimensions[:] = harray[:,0:3]
    GridStartIndices[:] = harray[:,3:6]
    GridEndIndices[:] = harray[:,6:9]
    GridLeftEdge[:] = harray[:,9:12]
    GridRightEdge[:] = harray[:,12:15]
    GridLevels[:] = harray[:,15:16]
    GridTimes[:] = harray[:,16:17]
    GridNumberOfParticles[:] = harray[:,17:18] */

  int counter=0;

  //fprintf(stderr, "counter: %d\n", counter++); // 0
  flat_dimensions[1] = 3;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_INT);
  PyDict_SetItemString(hierarchy_information, "GridDimensions", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 1
  flat_dimensions[1] = 3;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_INT);
  PyDict_SetItemString(hierarchy_information, "GridStartIndices", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 2
  flat_dimensions[1] = 3;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_INT);
  PyDict_SetItemString(hierarchy_information, "GridEndIndices", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 3
  flat_dimensions[1] = 3;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_PFLOAT);
  PyDict_SetItemString(hierarchy_information, "GridLeftEdge", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 4
  flat_dimensions[1] = 3;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_PFLOAT);
  PyDict_SetItemString(hierarchy_information, "GridRightEdge", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 5
  flat_dimensions[1] = 1; /* a bit iffy */
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_INT);
  PyDict_SetItemString(hierarchy_information, "GridLevels", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 6
  flat_dimensions[1] = 1;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_PFLOAT);
  PyDict_SetItemString(hierarchy_information, "GridTimes", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 7
  flat_dimensions[1] = 1;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_PFLOAT);
  PyDict_SetItemString(hierarchy_information, "GridOldTimes", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 8
  flat_dimensions[1] = 1;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_INT);
  PyDict_SetItemString(hierarchy_information, "GridProcs", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 9
  flat_dimensions[1] = 1;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_INT);
  PyDict_SetItemString(hierarchy_information, "GridNumberOfParticles", (PyObject *) temp_array);
  Py_DECREF(temp_array);

  //fprintf(stderr, "counter: %d\n", counter++); // 10
  flat_dimensions[1] = 1;
  temp_array = (PyArrayObject *) PyArray_SimpleNew(2, flat_dimensions, ENPY_INT);
  PyDict_SetItemString(hierarchy_information, "GridParentIDs", (PyObject *) temp_array);
  Py_DECREF(temp_array);

#endif
  //return SUCCESS;
}
