#ifdef USE_PYTHON
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

#ifndef ENZO_PYTHON_IMPORTED
#define PY_ARRAY_UNIQUE_SYMBOL enzo_ARRAY_API
#define ENZO_PYTHON_IMPORTED
#endif

#include <Python.h>
#include "numpy/arrayobject.h"

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "TopGridData.h"

int  GetUnits(float *DensityUnits, float *LengthUnits,
		       float *TemperatureUnits, float *TimeUnits,
		       float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int ExposeDataHierarchy(TopGridData *MetaData, HierarchyEntry *Grid, 
		       int &GridID, FLOAT WriteTime, int reset, int ParentID, int level);
void ExposeGridHierarchy(int NumberOfGrids);

static PyObject *_parameterFindingError;

static PyObject *PyGetEnzoParameter(PyObject *obj, PyObject *args)
{

    char *parameter_name;

    if (!PyArg_ParseTuple(args, "s", &parameter_name))
        return PyErr_Format(_parameterFindingError,
                    "FindBindingEnergy: Invalid parameters.");

    fprintf(stderr, "Looking for %s\n", parameter_name);
#include "InitializePythonInterface_finderfunctions.inc"

    return Py_None;
    
}

static PyMethodDef _EnzoModuleMethods[] = {
  {"get_parameter", PyGetEnzoParameter, METH_VARARGS},
  {NULL, NULL, 0, NULL}
};

int enzo_import_numpy()
{
    import_array1(FAIL);
    return SUCCESS;
}

int InitializePythonInterface(int argc, char *argv[])
{
#undef int
  static int PythonInterpreterInitialized = 0;
  PyObject *modstring, *mod;
  if(PythonInterpreterInitialized == 0){
    Py_SetProgramName("embed_enzo");
    Py_InitializeEx(0);
    if (!Py_IsInitialized()) ENZO_FAIL("Couldn't initialize Python!");
    PySys_SetArgv(argc, argv);
    int rv = enzo_import_numpy();
    if (rv == FAIL) {
        _import_array();
        PyErr_PrintEx(0);
        ENZO_FAIL("Couldn't import numpy.  Dying!");
    }
    PyRun_SimpleString("import sys\nsys.path.insert(0,'.')\nsys._parallel = True\n");
    PyRun_SimpleString("import gc\n");
    PythonInterpreterInitialized = 1;
  }
  static int PythonEnzoModuleInitialized = 0;
  if (PythonEnzoModuleInitialized == 0){
    PyObject *enzo_module, *enzo_module_dict; 
    enzo_module = Py_InitModule("enzo", _EnzoModuleMethods);
    enzo_module_dict = PyModule_GetDict(enzo_module);
    if(enzo_module == NULL) ENZO_FAIL("Failed on Enzo Module!");
    if(enzo_module_dict == NULL) ENZO_FAIL("Failed on Dict!");
    PyDict_SetItemString(enzo_module_dict, "grid_data", grid_dictionary);
    PyDict_SetItemString(enzo_module_dict, "old_grid_data", old_grid_dictionary);
    PyDict_SetItemString(enzo_module_dict, "hierarchy_information", hierarchy_information);
    PyDict_SetItemString(enzo_module_dict, "yt_parameter_file", yt_parameter_file);
    PyDict_SetItemString(enzo_module_dict, "conversion_factors", conversion_factors);
    PyDict_SetItemString(enzo_module_dict, "my_processor", my_processor);
    if (PyRun_SimpleString("import user_script\n")) ENZO_FAIL("Importing user_script failed!");
    if(debug)fprintf(stdout, "Completed Python interpreter initialization\n");
    PythonEnzoModuleInitialized = 1;
  }
  return SUCCESS;
}

int FinalizePythonInterface()
{
  Py_Finalize();
  if(debug)fprintf(stdout, "Completed Python interpreter finalization.\n");
  return SUCCESS;
}

#define TEMP_PYINT(A) Py_XDECREF(temp_int); temp_int = PyLong_FromLong((long) A);
#define TEMP_PYFLOAT(A) Py_XDECREF(temp_float); temp_float = PyFloat_FromDouble((double) A);
#define TEMP_PYSTRING(A) Py_XDECREF(temp_string); temp_string = PyString_FromString(A);

void ExportParameterFile(TopGridData *MetaData, FLOAT CurrentTime, FLOAT OldTime, 
			 float dtFixed)
{
  /* We need: */

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
    VelocityUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, CurrentTime);

  PyObject *temp_int = NULL;
  PyObject *temp_float = NULL;
  PyObject *temp_string = NULL;

  if (ComovingCoordinates) {
    FLOAT a, dadt, FinalRedshift, CurrentRedshift;
    CosmologyComputeExpansionFactor(MetaData->StopTime, &a, &dadt);

    FinalRedshift = (1 + InitialRedshift)/a - 1;

    /* Compute the current redshift (for information only). */

    CosmologyComputeExpansionFactor(CurrentTime, &a, &dadt);
    CurrentRedshift = (1 + InitialRedshift)/a - 1;

    TEMP_PYFLOAT(CurrentRedshift);
    PyDict_SetItemString(yt_parameter_file, "CosmologyCurrentRedshift", temp_float);

    TEMP_PYFLOAT(ComovingBoxSize);
    PyDict_SetItemString(yt_parameter_file, "CosmologyComovingBoxSize", temp_float);

    TEMP_PYFLOAT(OmegaMatterNow);
    PyDict_SetItemString(yt_parameter_file, "CosmologyOmegaMatterNow",temp_float);

    TEMP_PYFLOAT(OmegaLambdaNow);
    PyDict_SetItemString(yt_parameter_file, "CosmologyOmegaLambdaNow",temp_float);

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

  /*TEMP_PYSTRING(MetaData->MetaDataString);
  PyDict_SetItemString(yt_parameter_file, "MetaDataString", temp_string);*/

  TEMP_PYINT(HydroMethod);
  PyDict_SetItemString(yt_parameter_file, "HydroMethod", temp_int);

  TEMP_PYINT(DualEnergyFormalism);
  PyDict_SetItemString(yt_parameter_file, "DualEnergyFormalism", temp_int);

  TEMP_PYFLOAT(CurrentTime);
  PyDict_SetItemString(yt_parameter_file, "InitialTime", temp_float);

  TEMP_PYFLOAT(MetaData->StopTime);
  PyDict_SetItemString(yt_parameter_file, "StopTime", temp_float);

  TEMP_PYFLOAT(OldTime);
  PyDict_SetItemString(yt_parameter_file, "OldTime", temp_float);

  TEMP_PYFLOAT(dtFixed);
  PyDict_SetItemString(yt_parameter_file, "dtFixed", temp_float);

  TEMP_PYINT(ComovingCoordinates);
  PyDict_SetItemString(yt_parameter_file, "ComovingCoordinates", temp_int);

  TEMP_PYINT(MultiSpecies);
  PyDict_SetItemString(yt_parameter_file, "MultiSpecies", temp_int);

  TEMP_PYINT(MetaData->TopGridRank);
  PyDict_SetItemString(yt_parameter_file, "TopGridRank", temp_int);

  TEMP_PYINT(RefineBy);
  PyDict_SetItemString(yt_parameter_file, "RefineBy", temp_int);

  TEMP_PYINT(NumberOfPythonCalls);
  PyDict_SetItemString(yt_parameter_file, "NumberOfPythonCalls", temp_int);

  TEMP_PYINT(NumberOfPythonSubcycleCalls);
  PyDict_SetItemString(yt_parameter_file, "NumberOfPythonSubcycleCalls", temp_int);

  TEMP_PYINT(NumberOfPythonTopGridCalls);
  PyDict_SetItemString(yt_parameter_file, "NumberOfPythonTopGridCalls", temp_int);

  TEMP_PYFLOAT(CurrentMaximumDensity);
  PyDict_SetItemString(yt_parameter_file, "CurrentMaximumDensity", temp_float);

  TEMP_PYFLOAT(AngularVelocity);
  PyDict_SetItemString(yt_parameter_file, "AngularVelocity", temp_float);

  TEMP_PYFLOAT(VelocityGradient);
  PyDict_SetItemString(yt_parameter_file, "VelocityGradient", temp_float);

  TEMP_PYFLOAT(MetaData->dtDataDump);
  PyDict_SetItemString(yt_parameter_file, "dtDataDump", temp_float);

  TEMP_PYFLOAT(MetaData->TimeLastDataDump);
  PyDict_SetItemString(yt_parameter_file, "TimeLastDataDump", temp_float);

  TEMP_PYINT(MetaData->WroteData);
  PyDict_SetItemString(yt_parameter_file, "WroteData", temp_int);

  PyObject *tgd_tuple, *tgd0, *tgd1, *tgd2;
  
  /* Construct a tuple */
  tgd0 = PyFloat_FromDouble((double) DomainLeftEdge[0]);
  tgd1 = PyFloat_FromDouble((double) DomainLeftEdge[1]);
  tgd2 = PyFloat_FromDouble((double) DomainLeftEdge[2]);
  tgd_tuple = PyTuple_Pack(3, tgd0, tgd1, tgd2);
  PyDict_SetItemString(yt_parameter_file, "DomainLeftEdge", tgd_tuple);
  Py_XDECREF(tgd_tuple); Py_XDECREF(tgd0); Py_XDECREF(tgd1); Py_XDECREF(tgd2);

  tgd0 = PyFloat_FromDouble((double) DomainRightEdge[0]);
  tgd1 = PyFloat_FromDouble((double) DomainRightEdge[1]);
  tgd2 = PyFloat_FromDouble((double) DomainRightEdge[2]);
  tgd_tuple = PyTuple_Pack(3, tgd0, tgd1, tgd2);
  PyDict_SetItemString(yt_parameter_file, "DomainRightEdge", tgd_tuple);
  Py_XDECREF(tgd_tuple); Py_XDECREF(tgd0); Py_XDECREF(tgd1); Py_XDECREF(tgd2);

  tgd0 = PyLong_FromLong((long) MetaData->TopGridDims[0]);
  tgd1 = PyLong_FromLong((long) MetaData->TopGridDims[1]);
  tgd2 = PyLong_FromLong((long) MetaData->TopGridDims[2]);
  tgd_tuple = PyTuple_Pack(3, tgd0, tgd1, tgd2);
  PyDict_SetItemString(yt_parameter_file, "TopGridDimensions", tgd_tuple);
  Py_XDECREF(tgd_tuple); Py_XDECREF(tgd0); Py_XDECREF(tgd1); Py_XDECREF(tgd2);

  tgd0 = PyLong_FromLong((long) MetaData->LeftFaceBoundaryCondition[0]);
  tgd1 = PyLong_FromLong((long) MetaData->LeftFaceBoundaryCondition[1]);
  tgd2 = PyLong_FromLong((long) MetaData->LeftFaceBoundaryCondition[2]);
  tgd_tuple = PyTuple_Pack(3, tgd0, tgd1, tgd2);
  PyDict_SetItemString(yt_parameter_file, "LeftFaceBoundaryCondition", tgd_tuple);
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
  return;
}


#endif
