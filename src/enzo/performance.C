/*****************************************************************************
 *                                                                           *
 * Copyright 2006 James Bordner
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Board of Trustees of the University of Illinois            *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
//======================================================================
//
// File:        performance.C
//
// Description: Performance-related code
//
//----------------------------------------------------------------------
//
// Namespaces:  jb
//
//----------------------------------------------------------------------
//
// James Bordner
// UCSD
//
//======================================================================

#include "performance.h"
#include "macros_and_parameters.h"

#ifdef USE_LCAPERF

void lcaperfInitialize (int max_level)
{

  // Initialize lcaperf

  // Define lcaperf attributes

  lcaperf.new_attribute ("timestep", LCAPERF_INT);
  lcaperf.new_attribute ("level",    LCAPERF_INT);

  // Define lcaperf counters

  lcaperf.new_counter ("count-zones",      LCAPERF_COUNTER_TYPE_USER_ABS);
  lcaperf.new_counter ("count-ghosts",     LCAPERF_COUNTER_TYPE_USER_ABS);
  lcaperf.new_counter ("count-grids",      LCAPERF_COUNTER_TYPE_USER_ABS);
  lcaperf.new_counter ("count-particles",  LCAPERF_COUNTER_TYPE_USER_ABS);
  lcaperf.new_counter ("time-sim",         LCAPERF_COUNTER_TYPE_USER_ABS);

#ifdef USE_PAPI
  lcaperf.new_counter("fp-ops",LCAPERF_COUNTER_TYPE_PAPI);
#endif

  for (int level=0; level <= max_level; level++) {

    char lcaperf_counter_name[30];

    sprintf (lcaperf_counter_name,"count-zones-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,LCAPERF_COUNTER_TYPE_USER_ABS);

    sprintf (lcaperf_counter_name,"count-ghosts-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,LCAPERF_COUNTER_TYPE_USER_ABS);

    sprintf (lcaperf_counter_name,"count-grids-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,LCAPERF_COUNTER_TYPE_USER_ABS);

    sprintf (lcaperf_counter_name,"count-particles-local-%"ISYM,level);
    lcaperf.new_counter(lcaperf_counter_name,LCAPERF_COUNTER_TYPE_USER_ABS);

  }

#ifdef USE_LCAPERF_HDF5
  lcaperf.new_counter ("hdf5-read-calls",  LCAPERF_COUNTER_TYPE_USER_REL);
  lcaperf.new_counter ("hdf5-read-bytes",  LCAPERF_COUNTER_TYPE_USER_REL);
  lcaperf.new_counter ("hdf5-write-calls", LCAPERF_COUNTER_TYPE_USER_REL);
  lcaperf.new_counter ("hdf5-write-bytes", LCAPERF_COUNTER_TYPE_USER_REL);
#endif
}

#endif /* USE_LCAPERF */


