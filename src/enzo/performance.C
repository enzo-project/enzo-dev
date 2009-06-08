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

#ifdef USE_JBPERF

void jbPerfInitialize (int max_level)
{

  // Initialize jbPerf

  // Define jbPerf attributes

  jbPerf.new_attribute ("timestep", JB_INT);
#ifdef JB_PERF_LEVELS
  jbPerf.new_attribute ("level",    JB_INT);
#endif

  // Define jbPerf counters

  jbPerf.new_counter ("count-zones",      JB_COUNTER_TYPE_USER_ABS);
  jbPerf.new_counter ("count-ghosts",     JB_COUNTER_TYPE_USER_ABS);
  jbPerf.new_counter ("count-grids",      JB_COUNTER_TYPE_USER_ABS);
  jbPerf.new_counter ("count-particles",  JB_COUNTER_TYPE_USER_ABS);
  jbPerf.new_counter ("time-sim",         JB_COUNTER_TYPE_USER_ABS);

#ifdef USE_PAPI
  jbPerf.new_counter("fp-ops",JB_COUNTER_TYPE_PAPI);
#endif

  for (int level=0; level <= max_level; level++) {

    char jb_counter_name[30];

    sprintf (jb_counter_name,"count-zones-local-%"ISYM,level);
    jbPerf.new_counter(jb_counter_name,JB_COUNTER_TYPE_USER_ABS);

    sprintf (jb_counter_name,"count-ghosts-local-%"ISYM,level);
    jbPerf.new_counter(jb_counter_name,JB_COUNTER_TYPE_USER_ABS);

    sprintf (jb_counter_name,"count-grids-local-%"ISYM,level);
    jbPerf.new_counter(jb_counter_name,JB_COUNTER_TYPE_USER_ABS);

    sprintf (jb_counter_name,"count-particles-local-%"ISYM,level);
    jbPerf.new_counter(jb_counter_name,JB_COUNTER_TYPE_USER_ABS);

  }

#ifdef USE_JBPERF_HDF5
  jbPerf.new_counter ("hdf5-read-calls",  JB_COUNTER_TYPE_USER_REL);
  jbPerf.new_counter ("hdf5-read-bytes",  JB_COUNTER_TYPE_USER_REL);
  jbPerf.new_counter ("hdf5-write-calls", JB_COUNTER_TYPE_USER_REL);
  jbPerf.new_counter ("hdf5-write-bytes", JB_COUNTER_TYPE_USER_REL);
#endif
}

#endif /* USE_JBPERF */


