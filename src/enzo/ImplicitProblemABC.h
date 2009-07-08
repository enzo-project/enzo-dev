/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Implicit Problem Abstract Base Class, used in conjunction with 
/  nonlinear implicit solver
/
/  written by: Daniel Reynolds
/  date:       March, 2006
/  modified1:  
/
/  PURPOSE: This class defines problem-specific functions required for 
/  any implicit nonlinear solve (see InexactNewton.h for additional 
/  solver information).  For a given problem, the user must create a 
/  derived class that instantiates all of these operations on the 
/  appropriate data structures.
/
************************************************************************/

#ifndef IMPLICIT_PROBLEM_ABSTRACT_BASE_CLASS_DEFINED__
#define IMPLICIT_PROBLEM_ABSTRACT_BASE_CLASS_DEFINED__

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

class ImplicitProblemABC 
{

 public:

  // Destructor
  //  virtual ~ImplicitProblemABC() = 0;
  
  // Problem initializer
  virtual int Initialize(HierarchyEntry &TopGrid, 
			 TopGridData &MetaData) = 0;
  
  // Problem-specific Linear system setup function, sets up the 
  //   linear Newton system matrix J(u) given an updated state u.
  virtual int WriteParameters(FILE *fptr) = 0;
  
  // Problem-specific Linear solver function 
  //   solves ||J(u)*s - b|| to tolerance delta
  virtual int Evolve(HierarchyEntry *ThisGrid, float deltat) = 0;
  
};
  
#endif
