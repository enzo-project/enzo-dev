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
#define PY_ARRAY_UNIQUE_SYMBOL enzo_ARRAY_API
#include <Python.h>
#include "ProblemType_Python.h"

namespace{
    EnzoProblemType_creator_concrete<ProblemType_Python>
        python_initializer("PythonInitializer");
}

ProblemType_Python::ProblemType_Python() : EnzoProblemType() {
}

ProblemType_Python::~ProblemType_Python() { }

int ProblemType_Python::InitializeSimulation(FILE *pftr, FILE *Outfptr,
    HierarchyEntry &TopGrid, TopGridData &MetaData)
{
    int rv;
    //InitializePythonInterface(0, NULL);
    //import_numpy();
    //import_array1(FAIL);

#undef int
    Py_SetProgramName("embed_enzo");
    Py_Initialize();
    import_enzolib__problemtype_handler();
    print_hello();
    rv = create_problem_instance(this); 
}

#endif
#endif
