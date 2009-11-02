#ifdef __macros_and_parameters_h_
ERROR: need performance.h to be included before macros_and_parameters.h
#endif

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
#ifndef PERFORMANCE_H
#define PERFORMANCE_H

//======================================================================
//
// File:        performance.h
//
// Description: Interface layer between Enzo and lcaperf
//
// To use lcaperf, compile with -DUSE_LCAPERF and link with -llcaperf
// To use PAPI with lcaperf, compile also with -DUSE_PAPI and link with -lpapi
// To instrument HDF5 read/writes with lcaperf, compile with -DUSE_LCAPERF_HDF5
// To instrument MPI send/recvs, link with -lpmpi
//
//----------------------------------------------------------------------
//
// James Bordner (jobordner@ucsd.edu)
// 2003-06-20
//
//======================================================================

//----------------------------------------------------------------------
// lcaperf
//----------------------------------------------------------------------
#ifdef USE_LCAPERF

#   include "lcaperf.h" // THIS MUST COME BEFORE int IS REDEFINED

#   define LCAPERF_DUMP_FREQUENCY 1  /* How frequently to dump data to files. */

#endif /* USE_LCAPERF */

#ifdef USE_LCAPERF

#  define LCAPERF_BEGIN(segment)    lcaperf.begin (segment)
#  define LCAPERF_END(segment)      lcaperf.end (segment)
#  define LCAPERF_START(region)     lcaperf.start (region)
#  define LCAPERF_STOP(region)      lcaperf.stop (region)

#else

#  define LCAPERF_BEGIN(segment)    /* This space intentionally left blank */ ;
#  define LCAPERF_END(segment)      /* This space intentionally left blank */ ;
#  define LCAPERF_START(region)     /* This space intentionally left blank */ ;
#  define LCAPERF_STOP(region)      /* This space intentionally left blank */ ;

#endif

#endif /* PERFORMANCE_H */

