/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Null Implicit Problem Class (shell only)
/
/  written by: Daniel Reynolds
/  date:       June, 2009
/  modified:   
/
/  PURPOSE: This class defines the required ImplicitProblemABC functions
/           in case a user does not wish to solve any implicit problem.
/
************************************************************************/
#include "NullProblem.h"


// Constructor
NullProblem::NullProblem() { }

// Destructor
NullProblem::~NullProblem() { }

// Problem Initializer
int NullProblem::Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData) { 
  return SUCCESS; }
  
// Problem setup/solver
int NullProblem::Evolve(HierarchyEntry *ThisGrid, float deltat) { 
  return SUCCESS; }

// Write module parameters to file
int NullProblem::WriteParameters(FILE *fptr) { 
  return SUCCESS; }
