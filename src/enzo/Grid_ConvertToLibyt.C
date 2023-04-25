#ifdef USE_LIBYT
/***********************************************************************
/
/  GRID CLASS (GIVE OUR INFO TO LIBYT)
/
/  written by: Matthew Turk
/  date:       April, 2023
/
/  PURPOSE:
/
/  RETURNS:
/    ENZO_SUCCESS or FAIL
/
************************************************************************/

#include <string>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "libyt/libyt.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

/* Note that previously we allowed WriteTime to be supplied and we do not now */
void grid::ConvertToLibyt(int LocalGridID, int GlobalGridID, int ParentID, int level, yt_grid *GridInfoArray)
{

    this->DebugCheck("Converting to libyt arrays");

    /* Declarations */

    /* libyt only wants us to set any array values if the grid is local.  Note
     * that this is opposite what we do in other places, where operations on
     * the hierarchy (that don't involve fields) are often replicated on all
     * processes. */
  
    if (ProcessorNumber == MyProcessorNumber) return;

    for (int i = 0; i < MAX_DIMENSION; i++) {
        GridInfoArray[LocalGridID].grid_dimensions[i] = this->GridDimension[i];
        GridInfoArray[LocalGridID].left_edge[i] = this->GridLeftEdge[i];
        GridInfoArray[LocalGridID].right_edge[i] = this->GridRightEdge[i];
    }
    GridInfoArray[LocalGridID].id = GlobalGridID;
    GridInfoArray[LocalGridID].parent_id = ParentID;
    GridInfoArray[LocalGridID].level = level;
    /* par_count_list can take multiple particle types */
    GridInfoArray[LocalGridID].par_count_list[0] = this->ReturnNumberOfParticles();

    for (int field = 0; field < NumberOfBaryonFields; field++) {
        /* These are pointers, and *not* copies. Ownership is retained here. *
         * Note that we are only setting with the *current* BaryonField values,
         * and ignoring OldBaryonField.
         *
         * The other thing to note here is that we need to be extra careful
         * about fields being in different order in different grids.  So what we
         * do is use our global lookup table to map the field_type enum to the
         * libyt field entry.
         * 
         * This means we do not need to do any species field identification.
         *
         * (This may be far too cautious in almost all cases.  I mean,
         * grid::WriteGrid assumes that DataLabel corresponds to the field
         * index, and that is also what we assume internally for libyt
         * initialization.)
         *
         * */
        int libyt_field = libyt_field_lookup[field];
        if (libyt_field == -1) continue;
        GridInfoArray[LocalGridID].field_data[libyt_field].data_ptr = BaryonField[field];
    }

}
#endif
