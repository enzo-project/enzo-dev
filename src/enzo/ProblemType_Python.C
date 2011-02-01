/***********************************************************************
/
/  PYTHON PROBLEM TYPE INITIALIZER
/
/  written by: Matthew Turk
/  date:       January, 2011
/
/  PURPOSE:
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#ifdef USE_PYTHON

#ifndef ENZO_PYTHON_IMPORTED
#define PY_ARRAY_UNIQUE_SYMBOL enzo_ARRAY_API
#define ENZO_PYTHON_IMPORTED
#endif

#include <Python.h>
#include "numpy/arrayobject.h"
#include "ProblemType_Python.h"

namespace{
    EnzoProblemType_creator_concrete<ProblemType_Python>
        python_initializer("PythonInitializer");
}

char *argv[] = {"enzo.exe", "-d", "Something"};

ProblemType_Python::ProblemType_Python() : EnzoProblemType() {
}

ProblemType_Python::~ProblemType_Python() { }

int ProblemType_Python::InitializeSimulation(FILE *pftr, FILE *Outfptr,
    HierarchyEntry &TopGrid, TopGridData &MetaData)
{
    int rv;
#undef int
    if (PythonInterpreterInitialized == 0){
      Py_SetProgramName("embed_enzo");
      Py_Initialize();
      PySys_SetArgv(3, argv);
      import_array1(FAIL);
      PyRun_SimpleString("import sys\nsys.path.insert(0,'.')\nsys._parallel = True\n");
      PythonInterpreterInitialized = 1;
    }

    initproblemtype_handler();
    print_hello();
    rv = create_problem_instance(this); 
}

#endif
#endif
