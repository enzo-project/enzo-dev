/***********************************************************************
/
/  COMPUTE AND OUTPUT SOME SUMMARY INFORMATION ABOUT HIERARCHY
/
/  written by: Greg Bryan
/  date:       September, 1996
/  modified1:  Robert Harkness
/  date:       May, 2008
/
/  PURPOSE:
/
************************************************************************/
 
#include "performance.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
 
#if defined(MALLOC_IRIS4)
#include <sys/types.h>
#include <malloc.h>
#endif /* MALLOC_IRIS4 */
 
 
int OutputLevelInformation(FILE *fptr, TopGridData &MetaData,
			   LevelHierarchyEntry *LevelArray[])
{
 
  // NOTE: DON'T immediately return if this isn't the root
  //       processor, since lcaperf computes grid/zone-related
  //       counts here

  bool isRoot = (MyProcessorNumber == ROOT_PROCESSOR);

  /* Declarations */
 
  int level, maxdepth = 0, GridMemory, Cells, TotalCells, NumParticles;
  float AxialRatio, GridVolume;
  int Grids[MAX_DEPTH_OF_HIERARCHY];
  int LocalGrids[MAX_DEPTH_OF_HIERARCHY];

  long long Particles[MAX_DEPTH_OF_HIERARCHY];
  long long LocalParticles[MAX_DEPTH_OF_HIERARCHY];
  long long CellsActive[MAX_DEPTH_OF_HIERARCHY];
  long long CellsTotal[MAX_DEPTH_OF_HIERARCHY]; 
  long long LocalCellsActive[MAX_DEPTH_OF_HIERARCHY];
  long long LocalCellsTotal[MAX_DEPTH_OF_HIERARCHY]; 
  long long CellsFlagged[MAX_DEPTH_OF_HIERARCHY];

  float Coverage[MAX_DEPTH_OF_HIERARCHY], Memory[MAX_DEPTH_OF_HIERARCHY],
        MeanAxialRatio[MAX_DEPTH_OF_HIERARCHY],
        FractionFlagged[MAX_DEPTH_OF_HIERARCHY];

  LevelHierarchyEntry *Temp;
 
  /* Zero Hierarchy sums. */
 
  int   HierarchyGrids = 0;
  float HierarchyAxialRatio    = 0, HierarchyMemory = 0;
 
  /* Loop over levels, counting grids & approx memory used. */
 
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    /* Zero this level's summary info. */
 
    Grids[level]            = 0;
    Particles[level]        = 0;
    Memory[level]           = 0;
    Coverage[level]         = 0;
    MeanAxialRatio[level]   = 0;
    CellsActive[level]      = 0;
    CellsTotal[level]       = 0;
    CellsFlagged[level]     = 0;
    LocalCellsActive[level] = 0;
    LocalCellsTotal[level]  = 0;
    LocalGrids[level]       = 0;
    LocalParticles[level]   = 0;
 
    /* Loop over grids on this level. */
 
    Temp = LevelArray[level];
    while (Temp != NULL) {
 
      /* Incorporate grid information. */
 
      maxdepth = level;
 
      Temp->GridData->CollectGridInformation(GridMemory, GridVolume,
					     Cells, AxialRatio,
					     TotalCells,NumParticles);
      bool isLocal = Temp->GridData->isLocal();

      Grids[level]++;
      Particles[level]      += NumParticles;
      Memory[level]         += float(GridMemory);
      Coverage[level]       += GridVolume;
      CellsActive[level]    += Cells;
      CellsTotal[level]     += TotalCells;
      MeanAxialRatio[level] += AxialRatio;
      if (isLocal) {
	LocalGrids[level]++;
	LocalParticles[level]      += NumParticles;
	LocalCellsActive[level]    += Cells;
	LocalCellsTotal[level]     += TotalCells;
      }
 
      /* Get the actual number of flagged cells. */
 
#ifdef UNUSED
      int dummy;
      if (MetaData->StaticHierarchy == FALSE) {
	Temp->GridData->ClearFlaggingField();
	if (Temp->GridData->SetFlaggingField(dummy, level-1) == FAIL) {
	  ENZO_FAIL("Error in grid->SetFlaggingField!\n");
	}
	CellsFlagged[level] += Temp->GridData->FlagBufferZones();
	Temp->GridData->DeleteFlaggingField();
      }
#endif /* UNUSED */
 
      /* Next grid */
 
      Temp = Temp->NextGridThisLevel;
    }
 
    /* Add to overall sums. */
 
    HierarchyGrids         += Grids[level];
    HierarchyMemory        += Memory[level];
    HierarchyAxialRatio    += MeanAxialRatio[level];

    if (Grids[level] == 0) {
      MeanAxialRatio[level] = 0;
      FractionFlagged[level] = 0;
    } else {
      MeanAxialRatio[level] /= float(Grids[level]);
      FractionFlagged[level] = float(CellsFlagged[level])/
	float(max(CellsActive[level], 1));
    }

  }
 
  /* Get the total memory, if possible. */
 
  float TotalMemoryDeclared = 0, TotalMemoryUsed = 0;
#if defined(MALLOC_IRIS4)
  struct mallinfo proc;
  proc = mallinfo();
  TotalMemoryDeclared = float(proc.arena)/1.049e6;
  TotalMemoryUsed     = float(proc.usmblks+proc.uordblks)/1.049e6;
#endif /* MALLOC_IRIS4 */
 
  /* Write output (memory in MB). */

  if (isRoot) {
    fprintf(fptr, "Cycle %"ISYM"  Time %"GOUTSYM"  MaxDepth %"ISYM"  Grids %"ISYM"  Memory(MB) %"GSYM"  Ratio %"GSYM"\n",
      MetaData.CycleNumber,
      MetaData.Time, maxdepth, HierarchyGrids,
      float(HierarchyMemory)/1.049e6,
      HierarchyAxialRatio/float(HierarchyGrids) );
  }
 
#ifdef USE_LCAPERF
  long long jb_zones      = 0;
  long long jb_ghosts     = 0;
  long long jb_grids      = 0;
  long long jb_particles  = 0;
#endif

  for (level = 0; level <= MaximumRefinementLevel; level++) {

    if (isRoot) {
      fprintf(fptr, "  Level %"ISYM"  Grids %"ISYM"  Memory(MB) %"GSYM"  Coverage %"GSYM"  Ratio %"GSYM"  Flagged %"GSYM"  Active %lld\n",
      level, Grids[level], float(Memory[level])/1.049e6,
      Coverage[level], MeanAxialRatio[level], FractionFlagged[level],
      CellsActive[level]);
    }

#ifdef USE_LCAPERF

    // Assign to lcaperf counters "count-[zones|ghosts|grids]-<level>"

    if (level <= MaximumRefinementLevel) {

      char jb_counter_name[30];
      long long value;

      // lcaperf count-zones-local-<level>

      sprintf (jb_counter_name,"count-zones-local-%"ISYM,level);
      value = LocalCellsActive[level];
      lcaperf.assign(jb_counter_name, value);

      // lcaperf count-ghosts-local-<level>

      sprintf (jb_counter_name,"count-ghosts-local-%"ISYM,level);
      value = LocalCellsTotal[level]-LocalCellsActive[level];
      lcaperf.assign(jb_counter_name, value);

      // lcaperf count-grids-local-<level>

      sprintf (jb_counter_name,"count-grids-local-%"ISYM,level);
      value = LocalGrids[level];
      lcaperf.assign(jb_counter_name, value);

      // lcaperf count-particles-local-<level>

      sprintf (jb_counter_name,"count-particles-local-%"ISYM,level);
      value = LocalParticles[level];
      lcaperf.assign(jb_counter_name, value);
    }

    // Accumulate lcaperf counter values for "count-*"

    jb_zones     += CellsActive[level];
    jb_ghosts    += CellsTotal[level] - CellsActive[level];;
    jb_grids     += Grids[level];
    jb_particles += Particles[level];

#endif

  }

#ifdef USE_LCAPERF

  // Assign to lcaperf counters "count-[zones|ghosts|grids]"

  lcaperf.assign("count-zones",     jb_zones);
  lcaperf.assign("count-ghosts",    jb_ghosts);
  lcaperf.assign("count-grids",     jb_grids);
  lcaperf.assign("count-particles", jb_particles);
  
#endif

  if (isRoot) fprintf(fptr, "\n");

 
  return SUCCESS;
}
