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

#ifndef NULL_IMPLICIT_PROBLEM_DEFINED__
#define NULL_IMPLICIT_PROBLEM_DEFINED__

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "EnzoVector.h"
#include "ImplicitProblemABC.h"


class NullProblem : public virtual ImplicitProblemABC {

 public:

  // Constructor
  NullProblem();
  
  // Destructor
  ~NullProblem();

  // Problem Initializer
  int Initialize(HierarchyEntry &TopGrid, TopGridData &MetaData);
  
  // Problem setup/solver
  int Evolve(HierarchyEntry *ThisGrid, float deltat);

  // Write module parameters to file
  int WriteParameters(FILE *fptr);

};


#endif
