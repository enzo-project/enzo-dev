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
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
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

class ProblemType_Python;

class PythonGrid : private grid {
    friend class ProblemType_Python;
};

class ProblemType_Python : public EnzoProblemType
{
    public:
    ProblemType_Python();
    ~ProblemType_Python();
    int AddDataLabel(const char *FieldName);

    //private:
    int DataLabelCount;

};

#endif
