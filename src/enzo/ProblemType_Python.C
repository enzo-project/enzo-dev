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

void ProblemType_Python::SetField(PythonGrid *grid,
        int FieldIndex, float *data, int FieldType) {

    if (grid->BaryonField[FieldIndex] != NULL) {
        delete grid->BaryonField[FieldIndex];
    } else {
        /* We may not want to do this once we move to more types
           of field generation */
        grid->NumberOfBaryonFields++;
    }
    grid->BaryonField[FieldIndex] = data;
    grid->FieldType[FieldIndex] = FieldType;
    fprintf(stderr, "Seting %"ISYM" to %"ISYM"\n",
                FieldIndex, FieldType);
}

float* ProblemType_Python::GetField(PythonGrid *grid, int FieldIndex) {
    if (FieldIndex > MAX_NUMBER_OF_BARYON_FIELDS) return NULL;
    return grid->BaryonField[FieldIndex]; /* Should be NULL if not allocated */
}

void ProblemType_Python::GetGridInformation(
        PythonGrid *grid, int *ActiveDimensions,
        FLOAT *GridLeftEdge, FLOAT *GridRightEdge) {
        /*fprintf(stderr, "SI %"ISYM" %"ISYM" %"ISYM"\n",
                grid->GridStartIndex[0],
                grid->GridStartIndex[1],
                grid->GridStartIndex[2]);
        fprintf(stderr, "EI %"ISYM" %"ISYM" %"ISYM"\n",
                grid->GridEndIndex[0],
                grid->GridEndIndex[1],
                grid->GridEndIndex[2]);
        fprintf(stderr, "LE %"GSYM" %"GSYM" %"GSYM"\n",
                grid->GridLeftEdge[0],
                grid->GridLeftEdge[1],
                grid->GridLeftEdge[2]);
        fprintf(stderr, "RE %"GSYM" %"GSYM" %"GSYM"\n",
                grid->GridRightEdge[0],
                grid->GridRightEdge[1],
                grid->GridRightEdge[2]);
                */
    for (int i = 0; i < MAX_DIMENSION; i++) {
        ActiveDimensions[i] = grid->GridEndIndex[i] - grid->GridStartIndex[i] + 1;
        GridLeftEdge[i] = grid->GridLeftEdge[i];
        GridRightEdge[i] = grid->GridRightEdge[i];
    }
    return;
}

// All methods must go above this method

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
    fprintf(stderr, "ABOUT TO CREATE\n");
    PythonGrid *pgrid = static_cast<PythonGrid*> (TopGrid.GridData);
    fprintf(stderr, "Dimensions: %"ISYM" %"ISYM" %"ISYM"\n",
        pgrid->GridDimension[0],
        pgrid->GridDimension[1],
        pgrid->GridDimension[2]);
    rv = create_problem_instance(this, pgrid);

    return SUCCESS;
}

#endif
#endif
