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
void grid::ConvertToLibyt(int LocalGridID, int GlobalGridID, int ParentID, int level, yt_grid *GridInfo)
{

    this->DebugCheck("Converting to libyt arrays");

    /* Declarations */

    /* libyt only wants us to set any array values if the grid is local.  Note
     * that this is opposite what we do in other places, where operations on
     * the hierarchy (that don't involve fields) are often replicated on all
     * processes. 
     *
     * *But* there may be times when we actually do want to do this.  So,
     * we don't do the check in here -- instead, we rely on it being checked
     * externally.
     *
     * In here, we'll happily do whatever we're told.
     *
     * */
  
    for (int i = 0; i < MAX_DIMENSION; i++) {
        GridInfo->grid_dimensions[i] = this->GridDimension[i];
        GridInfo->left_edge[i] = this->GridLeftEdge[i];
        GridInfo->right_edge[i] = this->GridRightEdge[i];
    }
    GridInfo->id = GlobalGridID;
    GridInfo->parent_id = ParentID;
    GridInfo->level = level;
    /* par_count_list can take multiple particle types */
    GridInfo->par_count_list[0] = this->ReturnNumberOfParticles();

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
        GridInfo->field_data[libyt_field].data_ptr = BaryonField[field];
    }

}
#endif
