/****************************************************************************
 *                                                                          *
 * Copyright 2009 Laboratory for Computational Astrophysics                 *
 * Copyright 2009 Regents of the University of California                   *
 *                                                                          *
 * This software is released under the terms of the "Enzo Public License"   *
 * in the accompanying LICENSE file.                                        *
 *                                                                          *
 ****************************************************************************/
/***********************************************************************
/  Contains headers that must be included *prior* to inclusion
/  of macros_and_parameters.h, since that so beautifully re-defines 
/  'float', messing up any external library that actually uses 'float' 
/  in its header files.  
/
/  written by: Daniel R. Reynolds
/  date:       June, 2009
/  modified1:  
/
************************************************************************/

#ifndef IMPLICIT_PROBLEM_PREINCLUDES_DEFINED__
#define IMPLICIT_PROBLEM_PREINCLUDES_DEFINED__

#ifdef USE_MPI
#include "mpi.h"
#endif
#ifdef USE_GRACKLE
extern "C" {
#include <grackle.h>
}
#endif
/* #include <stdlib.h> */
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifdef USE_HYPRE
#include "HYPRE_sstruct_ls.h"
#endif
#include "performance.h"
#include "ErrorExceptions.h"

#endif
