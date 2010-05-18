/***********************************************************************
/
/  ENZO Version and Options in Effect
/
/  written by: Robert Harkness
/  date:       May 2008
/
************************************************************************/
 
#ifdef USE_MPI
#include <mpi.h>
#endif
 
#include <stdio.h>
#include <strings.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"
#include "svn_version.def"
#include "fortran.def"
 
#ifdef MEM_TRACE
Eint64 mused(void);
#endif
 
 
int ENZO_OptionsinEffect(void) 
{

  FILE *opf;

  if (MyProcessorNumber == 0) {

    opf = fopen("Enzo_Options", "w");

    fprintf(opf, "ENZO Options in Effect\n");

    if (ENZO_SVN_REVISION != 0) {
      fprintf(opf,"=========================\n");
      fprintf(opf,"Enzo SVN Branch   %s\n",ENZO_SVN_BRANCH);
      fprintf(opf,"Enzo SVN Revision %s\n",ENZO_SVN_REVISION);
      fprintf(opf,"=========================\n");
    }

#ifdef SMALL_INTS
    fprintf(opf, " 32 bit Integer version\n");
#endif

#ifdef LARGE_INTS
    fprintf(opf, " 64 bit Integer version\n");
#endif

#ifdef CONFIG_PINT_4
    fprintf(opf, " 32 bit Integers for particle indices\n");
#endif
#ifdef CONFIG_PINT_8
    fprintf(opf, " 64 bit Integers for particle indices\n");
#endif

#ifdef INITS32
    fprintf(opf, " 32 bit Integer initial conditions\n");
#endif

#ifdef INITS64
    fprintf(opf, " 64 bit Integer initial conditions\n");
#endif

#ifdef CONFIG_BFLOAT_4
    fprintf(opf, " Float precision is 32 bits\n");
#endif

#ifdef CONFIG_BFLOAT_8
    fprintf(opf, " Float precision is 64 bits\n");
#endif

#ifdef CONFIG_PFLOAT_4
    fprintf(opf, " Position and time precision is 32 bits - NOT SUPPORTED!\n");
#endif

#ifdef CONFIG_PFLOAT_8
    fprintf(opf, " Position and time precision is 64 bits\n");
#endif

#ifdef CONFIG_PFLOAT_16
    fprintf(opf, " Position and time precision is 128 bits\n");
#endif

    fprintf(opf, "\n");
    fprintf(opf, "Optimizations in Effect\n");

#ifdef OOC_BOUNDARY
    fprintf(opf, "  Out-of-core Top Grid boundary conditions\n");
#endif

#ifdef FAST_SIB
    fprintf(opf, "  Fast Sibling Locator 1\n");
#endif

#ifdef FAST_SIB
    fprintf(opf, "  Fast Sibling Locator 2\n");
#endif

#ifdef FAST_SIB
    fprintf(opf, "  Fast Sibling Locator 3\n");
#endif

#ifdef FAST_SIB
    fprintf(opf, "  Fast Sibling Locator 4\n");
#endif

#ifdef FAST_SIB
    fprintf(opf, "  Fast Sibling Locator 5\n");
#endif

#ifdef STATIC_SIBLING_LIST
    fprintf(opf, "  Static allocation of Level Zero Sibling List\n");
#endif

#ifdef FLUX_FIX
    fprintf(opf, "  New Flux Correction scheme by Collins & Wagner\n");
#endif

#ifdef SAB
    fprintf(opf, "  AccelerationHack by Collins\n");
#endif

#ifdef USE_DT_LIMIT
    fprintf(opf, "  Use dt limit in AMR\n");
#endif

#ifdef UNIGRID
    fprintf(opf, "  Minimum memory start-up => non-apative mesh only\n");
#endif

#ifdef HDF5_USE_HDF5_GROUPS
    fprintf(opf, "  HDF5 groups for packed AMR\n");
#endif

#ifdef USE_HDF5_OUTPUT_BUFFERING
    fprintf(opf, "  HDF5 in-core buffering for packed AMR output\n");
#endif

#ifdef USE_HDF5_INPUT_BUFFERING
    fprintf(opf, "  HDF5 in-core buffering for packed AMR input\n");
#endif

#ifdef SINGLE_HDF5_OPEN_ON_INPUT
    fprintf(opf, "  HDF5 single open on input - no explicit task map\n");
#endif

#ifdef USE_NEW_RESTART_FORMAT
    fprintf(opf, "  New hierarchy format\n")
#endif

#ifdef FORCE_BUFFER_PURGE
    fprintf(opf, "  Force purge of communication buffers\n");
#endif

#ifdef FORCE_MSG_PROGRESS
    fprintf(opf, "  Force message progress with MPI_Barrier calls\n");
#endif

#ifdef TRANSFER
    fprintf(opf, "  Adaptive ray tracing enabled\n");
#else
    fprintf(opf, "  Adaptive ray tracing disabled\n");
#endif

#ifdef USE_PYTHON
    fprintf(opf, "  Inline python enabled\n");
#else
    fprintf(opf, "  Inline python disabled\n");
#endif

#ifdef FAST_SIB
    fprintf(opf, "  Fast sibiling search enabled\n");
#else
    fprintf(opf, "  Fast sibiling search disabled\n");
#endif

#ifdef USE_HDF4
    fprintf(opf, "  HDF4 reading enabled\n");
#else
    fprintf(opf, "  HDF4 reading disabled\n");
#endif

#ifdef FLUX_FIX
    fprintf(opf, "  Flux fix for subgrid siblings enabled\n");
#else
    fprintf(opf, "  Flux fix for subgrid siblings disabled\n");
#endif

#ifdef NEW_GRID_IO
    fprintf(opf, "  New Grid I/O enabled\n");
#else
    fprintf(opf, "  New Grid I/O disabled\n");
#endif

#ifdef BITWISE_IDENTICALITY
    fprintf(opf, "  Bitwise-identicality enabled\n");
#else
    fprintf(opf, "  Bitwise-identicality disabled\n");
#endif


    fprintf(opf, "\n");
    fprintf(opf, "Macro and Parameter Definitions\n");

    fprintf(opf, "  MAX_NUMBER_OF_TASKS                 %8d\n", MAX_NUMBER_OF_TASKS);
    fprintf(opf, "  MAX_NUMBER_OF_NODES                 %8d\n", MAX_NUMBER_OF_NODES);
#ifdef ENABLE_TASKMAP
    fprintf(opf, "  MAX_TASKS_PER_NODE                  %8d\n", MAX_TASKS_PER_NODE);
#endif
    fprintf(opf, "  MAX_NUMBER_OF_BARYON_FIELDS (>=6)   %8d\n", MAX_NUMBER_OF_BARYON_FIELDS);
    fprintf(opf, "  MAX_NUMBER_OF_SUBGRIDS              %8d\n", MAX_NUMBER_OF_SUBGRIDS);
    fprintf(opf, "  MAX_DEPTH_OF_HIERARCHY              %8d\n", MAX_DEPTH_OF_HIERARCHY);
    fprintf(opf, "  MAX_LINE_LENGTH                     %8d\n", MAX_LINE_LENGTH);

    fprintf(opf, "  MAX_NAME_LENGTH                     %8d\n", MAX_NAME_LENGTH);
    fprintf(opf, "  MAX_GRID_TAG_SIZE                   %8d\n", MAX_GRID_TAG_SIZE);
    fprintf(opf, "  MAX_TASK_TAG_SIZE                   %8d\n", MAX_TASK_TAG_SIZE);
    fprintf(opf, "  MAX_GROUP_TAG_SIZE                  %8d\n", MAX_GROUP_TAG_SIZE);
    fprintf(opf, "  MAX_CYCLE_TAG_SIZE                  %8d\n", MAX_CYCLE_TAG_SIZE);

    fprintf(opf, "  GRID FORMAT                              "); fprintf(opf, GRID_TAG_FORMAT);  fprintf(opf, "\n");
    fprintf(opf, "  TASK FORMAT                              "); fprintf(opf, TASK_TAG_FORMAT);  fprintf(opf, "\n");
    fprintf(opf, "  GROUP FORMAT                             "); fprintf(opf, GROUP_TAG_FORMAT); fprintf(opf, "\n");
    fprintf(opf, "  CYCLE FORMAT                             "); fprintf(opf, CYCLE_TAG_FORMAT); fprintf(opf, "\n");

    fprintf(opf, "  MAX_COUNTERS                        %8d\n", MAX_COUNTERS);
    fprintf(opf, "  DEFAULT_GHOST_ZONES (>=3)           %8d\n", DEFAULT_GHOST_ZONES);
    fprintf(opf, "  MAX_NUMBER_OF_OUTPUT_REDSHIFTS      %8d\n", MAX_NUMBER_OF_OUTPUT_REDSHIFTS);
    fprintf(opf, "  GRAVITY_BUFFER_SIZE                 %8d\n", GRAVITY_BUFFER_SIZE);
    fprintf(opf, "  MAX_FLAGGING_METHODS                %8d\n", MAX_FLAGGING_METHODS);
    fprintf(opf, "  MAX_STATIC_REGIONS                  %8d\n", MAX_STATIC_REGIONS);
    fprintf(opf, "  MAX_NUMBER_OF_PARTICLE_ATTRIBUTES   %8d\n", MAX_NUMBER_OF_PARTICLE_ATTRIBUTES);
    fprintf(opf, "  MAX_TIME_ACTIONS                    %8d\n", MAX_TIME_ACTIONS);
    fprintf(opf, "  MAX_CUBE_DUMPS                      %8d\n", MAX_CUBE_DUMPS);
    fprintf(opf, "  MAX_POTENTIAL_ITERATIONS            %8d\n", MAX_POTENTIAL_ITERATIONS);
    fprintf(opf, "  MAX_COLOR                           %8d\n", MAX_COLOR);
    fprintf(opf, "  MAX_ANY_SINGLE_DIRECTION            %8d\n", MAX_ANY_SINGLE_DIRECTION);

    fprintf(opf, "\n");

#ifdef SUN
    fprintf(opf, "Processor type is SUN\n");
#endif

#ifdef SUN_OLD
    fprintf(opf, "Processor type is SUN_OLD\n");
#endif

#ifdef SPP
    fprintf(opf, "Processor type is SPP\n");
#endif

#ifdef SP2
    fprintf(opf, "Processor type is SP2\n");
#endif

#ifdef BGL
    fprintf(opf, "Processor type is BGL\n");
#endif

#ifdef IRIS4
    fprintf(opf, "Processor type is IRIS4\n");
#endif

#ifdef COMPAQ
    fprintf(opf, "Processor type is COMPAQ\n");
#endif

#ifdef CONVEX
    fprintf(opf, "Processor type is CONVEX\n");
#endif

#ifdef LINUX
    fprintf(opf, "Processor type is LINUX\n");
#endif

#ifdef IA64
    fprintf(opf, "Processor type is IA64\n");
#endif

#ifdef CRAYX1
    fprintf(opf, "Processor type is CRAYX1\n");
#endif

#ifdef XT3
    fprintf(opf, "Processor type is CRAY XT3 or XT4\n");
#endif

#ifdef MEM_TRACE
    fprintf(opf, "Memory tracing enabled\n");
#endif

#ifdef TP_VELOCITY
    fprintf(opf, "Tracer particle velocity output is enabled\n");
#else
    fprintf(opf, "Tracer particle velocity output is disabled\n");
#endif

    fclose(opf);

  } // processor zero only

  return SUCCESS;

}
