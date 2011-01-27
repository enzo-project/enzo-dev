/***********************************************************************
/
/  GENERAL PROBLEM TYPE FOR MANIPULATION
/
/  written by: Matthew Turk
/  date:       January, 2011
/
/  PURPOSE:
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include "ProblemType_General.h"

ProblemType_General::ProblemType_General() : EnzoProblemType() { 
    this->DataLabelCount = 0;
}

ProblemType_General::~ProblemType_General() { }

int ProblemType_General::AddDataLabel(const char *FieldName) {
    /* We allocate a new copy of FieldName */
    /* Include NUL-terminator */
    int slen = strlen(FieldName) + 1;
    char *fcopy = new char[slen];
    strncpy(fcopy, FieldName, slen);
    DataLabel[this->DataLabelCount] = fcopy;
    return DataLabelCount++;
}

#endif
