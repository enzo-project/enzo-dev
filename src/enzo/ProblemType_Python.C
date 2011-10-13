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

#include <Python.h>
#include "numpy/arrayobject.h"
#include "ProblemType_Python.h"

void ExportParameterFile(TopGridData *MetaData, FLOAT CurrentTime);
int InitializePythonInterface(int argc, char **argv);

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
    for (int i = 0; i < MAX_DIMENSION; i++) {
        ActiveDimensions[i] = grid->GridEndIndex[i] - grid->GridStartIndex[i] + 1;
        GridLeftEdge[i] = grid->GridLeftEdge[i];
        GridRightEdge[i] = grid->GridRightEdge[i];
    }
    return;
}

int ProblemType_Python::InitializeSimulation(
    HierarchyEntry &TopGrid, TopGridData &MetaData)
{
    return SUCCESS;
}

// All methods must go above this method

int ProblemType_Python::InitializeSimulation(FILE *pftr, FILE *Outfptr,
    HierarchyEntry &TopGrid, TopGridData &MetaData)
{
    int rv;

    char *argv[] = {"enzo"};
    InitializePythonInterface(1, argv);

    initproblemtype_handler();
    PythonGrid *pgrid = static_cast<PythonGrid*> (TopGrid.GridData);
    pgrid->Level = 0;

    ExportParameterFile(&MetaData, MetaData.Time);

    rv = create_problem_instance(this, pgrid);

    return SUCCESS;
}

#endif
#endif
