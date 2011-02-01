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
#include "enzolib/enzolib/problemtype_handler_api.h"

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
};

class ProblemType_Python : public EnzoProblemType
{
    public:
    ProblemType_Python();
    ~ProblemType_Python();
    /*PyObject *problem_creator_wrapper;*/

    virtual int InitializeSimulation(FILE *fptr, FILE *Outfptr,
            HierarchyEntry &TopGrid, TopGridData &MetaData);

};

#endif /* USE_PYTHON */
#endif
