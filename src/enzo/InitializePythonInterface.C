#ifdef CONFIG_PYTHON_ENABLED
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

#ifdef CONFIG_PYTHON_ENABLED
#define PY_ARRAY_UNIQUE_SYMBOL enzo_ARRAY_API
#include <Python.h>
#include "numpy/arrayobject.h"
#define ENZO_PYTHON_IMPORTED
#endif

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

int  CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
		       float *TemperatureUnits, float *TimeUnits,
		       float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


static PyMethodDef _EnzoModuleMethods[] = {
  {NULL, NULL, 0, NULL}
};

int InitializePythonInterface(int argc, char *argv[])
{
#undef int

  //Py_SetProgramName(argv[0]);
  Py_SetProgramName("embed_enzo");

  Py_Initialize();

  PySys_SetArgv(argc, argv);
  PyRun_SimpleString("import sys\nsys.path.insert(0,'.')\nsys._parallel = True\n");
  PyObject *enzo_module, *enzo_module_dict; 
  enzo_module = Py_InitModule("enzo", _EnzoModuleMethods);
  enzo_module_dict = PyModule_GetDict(enzo_module);
  if(enzo_module == NULL){fprintf(stderr, "Failed on Enzo Module!\n");return FAIL;}
  if(enzo_module_dict == NULL){fprintf(stderr, "Failed on Dict!\n");return FAIL;}
  PyDict_SetItemString(enzo_module_dict, "grid_data", grid_dictionary);
  PyDict_SetItemString(enzo_module_dict, "old_grid_data", old_grid_dictionary);
  PyDict_SetItemString(enzo_module_dict, "hierarchy_information", hierarchy_information);
  PyDict_SetItemString(enzo_module_dict, "yt_parameter_file", yt_parameter_file);
  PyDict_SetItemString(enzo_module_dict, "conversion_factors", conversion_factors);
  import_array1(FAIL);
  npy_intp flat_dimensions[1];
  flat_dimensions[0] = (npy_intp) 5;
  PyObject *temp = PyArray_ZEROS(1, flat_dimensions, NPY_INT, 0);
  fprintf(stderr, "Completed initialization\n");
  return SUCCESS;
}


#define TEMP_PYINT(A) Py_XDECREF(temp_int); temp_int = PyLong_FromLong((long) A);
#define TEMP_PYFLOAT(A) Py_XDECREF(temp_float); temp_float = PyFloat_FromDouble((double) A);
#define TEMP_PYSTRING(A) Py_XDECREF(temp_string); temp_string = PyString_FromString(A);

int ExportParameterFile(TopGridData &MetaData, FLOAT CurrentTime)
{
  /* We need: */

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
        VelocityUnits = 1;

  if (ComovingCoordinates)
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, CurrentTime) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }

  FLOAT a, dadt, FinalRedshift, CurrentRedshift;
  if (CosmologyComputeExpansionFactor(MetaData.StopTime, &a, &dadt) == FAIL) {
    fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
    return FAIL;
  }
  FinalRedshift = (1 + InitialRedshift)/a - 1;

  /* Compute the current redshift (for information only). */

  CosmologyComputeExpansionFactor(CurrentTime, &a, &dadt);
  CurrentRedshift = (1 + InitialRedshift)/a - 1;


  PyObject *temp_int = NULL;
  PyObject *temp_float = NULL;
  PyObject *temp_string = NULL;

  if (ComovingCoordinates) {
    TEMP_PYFLOAT(CurrentRedshift);
    PyDict_SetItemString(yt_parameter_file, "CosmologyCurrentRedshift", temp_float);

    TEMP_PYFLOAT(ComovingBoxSize);
    PyDict_SetItemString(yt_parameter_file, "CosmologyComovingBoxSize", temp_float);

    TEMP_PYFLOAT(OmegaMatterNow);
    PyDict_SetItemString(yt_parameter_file, "CosmologyOmegaMatterNow",temp_float);

    TEMP_PYFLOAT(HubbleConstantNow);
    PyDict_SetItemString(yt_parameter_file, "CosmologyHubbleConstantNow",temp_float);

    TEMP_PYFLOAT(InitialRedshift);
    PyDict_SetItemString(yt_parameter_file, "CosmologyInitialRedshift",temp_float);
  }

  TEMP_PYFLOAT(DensityUnits);
  PyDict_SetItemString(yt_parameter_file, "DensityUnits", temp_float);

  TEMP_PYFLOAT(LengthUnits);
  PyDict_SetItemString(yt_parameter_file, "LengthUnits", temp_float);

  TEMP_PYFLOAT(TemperatureUnits);
  PyDict_SetItemString(yt_parameter_file, "TemperatureUnits", temp_float);

  TEMP_PYFLOAT(TimeUnits);
  PyDict_SetItemString(yt_parameter_file, "TimeUnits", temp_float);

  /*TEMP_PYSTRING(MetaData.MetaDataString);
  PyDict_SetItemString(yt_parameter_file, "MetaDataString", temp_string);*/

  TEMP_PYINT(HydroMethod);
  PyDict_SetItemString(yt_parameter_file, "HydroMethod", temp_int);

  TEMP_PYINT(DualEnergyFormalism);
  PyDict_SetItemString(yt_parameter_file, "DualEnergyFormalism", temp_int);

  TEMP_PYFLOAT(CurrentTime);
  PyDict_SetItemString(yt_parameter_file, "InitialTime", temp_float);

  TEMP_PYINT(ComovingCoordinates);
  PyDict_SetItemString(yt_parameter_file, "ComovingCoordinates", temp_int);

  TEMP_PYINT(MultiSpecies);
  PyDict_SetItemString(yt_parameter_file, "MultiSpecies", temp_int);

  TEMP_PYINT(MetaData.TopGridRank);
  PyDict_SetItemString(yt_parameter_file, "TopGridRank", temp_int);

  TEMP_PYINT(RefineBy);
  PyDict_SetItemString(yt_parameter_file, "RefineBy", temp_int);

  TEMP_PYINT(NumberOfPythonCalls);
  PyDict_SetItemString(yt_parameter_file, "NumberOfPythonCalls", temp_int);

  PyObject *tgd_tuple, *tgd0, *tgd1, *tgd2;
  
  /* Construct a tuple */
  tgd0 = PyLong_FromLong((long) DomainLeftEdge[0]);
  tgd1 = PyLong_FromLong((long) DomainLeftEdge[1]);
  tgd2 = PyLong_FromLong((long) DomainLeftEdge[2]);
  tgd_tuple = PyTuple_Pack(3, tgd0, tgd1, tgd2);
  PyDict_SetItemString(yt_parameter_file, "DomainLeftEdge", tgd_tuple);
  Py_XDECREF(tgd_tuple); Py_XDECREF(tgd0); Py_XDECREF(tgd1); Py_XDECREF(tgd2);

  tgd0 = PyLong_FromLong((long) DomainRightEdge[0]);
  tgd1 = PyLong_FromLong((long) DomainRightEdge[1]);
  tgd2 = PyLong_FromLong((long) DomainRightEdge[2]);
  tgd_tuple = PyTuple_Pack(3, tgd0, tgd1, tgd2);
  PyDict_SetItemString(yt_parameter_file, "DomainRightEdge", tgd_tuple);
  Py_XDECREF(tgd_tuple); Py_XDECREF(tgd0); Py_XDECREF(tgd1); Py_XDECREF(tgd2);

  tgd0 = PyLong_FromLong((long) MetaData.TopGridDims[0]);
  tgd1 = PyLong_FromLong((long) MetaData.TopGridDims[1]);
  tgd2 = PyLong_FromLong((long) MetaData.TopGridDims[2]);
  tgd_tuple = PyTuple_Pack(3, tgd0, tgd1, tgd2);
  PyDict_SetItemString(yt_parameter_file, "TopGridDimensions", tgd_tuple);
  Py_XDECREF(tgd_tuple); Py_XDECREF(tgd0); Py_XDECREF(tgd1); Py_XDECREF(tgd2);

  /* Do the conversion factors */
  for (int dim = 0; dim < MAX_NUMBER_OF_BARYON_FIELDS; dim++) {
    if (DataLabel[dim]) {
      if (strstr(DataLabel[dim], "Density") != NULL) {
        TEMP_PYFLOAT(DensityUnits);
	    PyDict_SetItemString(conversion_factors, DataLabel[dim], temp_float);
        }
      if (strstr(DataLabel[dim], "velocity") != NULL) {
        TEMP_PYFLOAT(VelocityUnits);
	    PyDict_SetItemString(conversion_factors, DataLabel[dim], temp_float);
        }
      if (strstr(DataLabel[dim], "_kph") != NULL)  {
        TEMP_PYFLOAT(1./TimeUnits);
	    PyDict_SetItemString(conversion_factors, DataLabel[dim], temp_float);
        }
    }
  }


  Py_XDECREF(temp_int); Py_XDECREF(temp_float); Py_XDECREF(temp_string);
  return SUCCESS;
}

int CallPython(TopGridData &MetaData, FLOAT CurrentTime)
{
  ExportParameterFile(MetaData, CurrentTime);

  PyRun_SimpleString("import user_script\nuser_script.main()\n");

  NumberOfPythonCalls++;

  PyDict_Clear(grid_dictionary);
  PyDict_Clear(hierarchy_information);
  PyDict_Clear(yt_parameter_file);
  PyDict_Clear(conversion_factors);
  return SUCCESS;
}

#endif
