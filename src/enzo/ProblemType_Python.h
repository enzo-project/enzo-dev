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
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

class ProblemType_Python;

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "TopGridData.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ProblemType.h"
#include "EventHooks.h"

class PythonGrid : private grid {
    friend class ProblemType_Python;
    public:
        int Level;
};

/* 
This next bit of sanitization is required because of Cython's difficulty
with C++ code mixed with public APIs.  It may some day become obsolete.
*/
#ifdef __cplusplus
#undef __cplusplus
#define __reset_cplusplus
#endif
#include "python_bridge/problemtype_handler.h"
#ifdef __reset_cplusplus
#define __cplusplus
#endif

#define PyArray_NOOWNDATA(obj) (((PyArrayObject *)(obj))->flags &= ~NPY_OWNDATA)

class ProblemType_Python : public EnzoProblemType
{
    public:
    ProblemType_Python();
    ~ProblemType_Python();
    /*PyObject *problem_creator_wrapper;*/

    virtual int InitializeSimulation(FILE *fptr, FILE *Outfptr,
            HierarchyEntry &TopGrid, TopGridData &MetaData);

    void SetField(PythonGrid *grid, int FieldIndex, float *data, int FieldType);
    float* GetField(PythonGrid *grid, int FieldIndex);

    void GetGridInformation(PythonGrid *grid,
        int ActiveDimensions[MAX_DIMENSION],
        FLOAT GridLeftEdge[MAX_DIMENSION],
        FLOAT GridRightEdge[MAX_DIMENSION]);

};

#endif /* USE_PYTHON */
#endif
