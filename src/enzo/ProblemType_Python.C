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
#include "ProblemType_Python.h"

ProblemType_Python::ProblemType_Python() : EnzoProblemType() { 
    this->DataLabelCount = 0;
}

ProblemType_Python::~ProblemType_Python() { }

int ProblemType_Python::AddDataLabel(const char *FieldName) {
    /* We allocate a new copy of FieldName */
    /* Include NUL-terminator */
    int slen = strlen(FieldName) + 1;
    char *fcopy = new char[slen];
    strncpy(fcopy, FieldName, slen);
    DataLabel[this->DataLabelCount] = fcopy;
    return DataLabelCount++;
}

#endif
