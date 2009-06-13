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
// Description: Interface layer between Enzo and jbPerf
//
// To use jbPerf, compile with -DUSE_JBPERF and link with -ljbperf
// To use PAPI with jbPerf, compile also with -DUSE_PAPI and link with -lpapi
// To instrument HDF5 read/writes with jbPerf, compile with -DUSE_JBPERF_HDF5
// To instrument MPI send/recvs, link with -lpmpi
//
//----------------------------------------------------------------------
//
// James Bordner (jbordner@cosmos.ucsd.edu)
// 2003-06-20
//
//======================================================================

//----------------------------------------------------------------------
// jbPerf
//----------------------------------------------------------------------
#ifdef USE_JBPERF

#   include "jbPerf.h" // THIS MUST COME BEFORE int IS REDEFINED

#   define JB_ITER_PER_SEGMENT 1  /* How frequently to dump data to files. */

#endif /* USE_JBPERF */

#ifdef USE_JBPERF

#  define JBPERF_BEGIN(segment)    jbPerf.begin (segment)
#  define JBPERF_END(segment)      jbPerf.end (segment)
#  define JBPERF_START(region)     jbPerf.start (region)
#  define JBPERF_STOP(region)      jbPerf.stop (region)

#else

#  define JBPERF_BEGIN(segment)    /* This space intentionally left blank */ ;
#  define JBPERF_END(segment)      /* This space intentionally left blank */ ;
#  define JBPERF_START(region)     /* This space intentionally left blank */ ;
#  define JBPERF_STOP(region)      /* This space intentionally left blank */ ;

#endif

#endif /* PERFORMANCE_H */

