/***********************************************************************
/
/  GRID CLASS (WRITE OUT HDF5 DATASPACES FOR PARTICLE TYPE)
/
/  written by: Matthew Turk
/  date:       August, 2009
/
/  PURPOSE:
/
************************************************************************/
 
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int grid::CreateParticleTypeGrouping(hid_t ptype_dset, 
                                     hid_t ptype_dspace,
                                     hid_t parent_group,
                                     hid_t file_id)
{
    int CurrentParticleType = ParticleType[0];
    int i;
    hsize_t start[1], count[1];
    hid_t new_dset, new_dspace;
    hdset_reg_ref_t reference;
    hsize_t TempOne[1] = {1};
    
    herr_t err = 0;

    start[0] = 0;
    char dset_name[22];
    char ptype_dset_fullname[255];
    H5Iget_name(ptype_dset, ptype_dset_fullname, 255);

    /* We iterate up to NoP because we check that inside the loop.
       If we are here, particles have been sorted by ParticleNumber. */
    for (i = 0; i <= NumberOfParticles; i++) {
        if ((i == NumberOfParticles)
         || (ParticleType[i] != CurrentParticleType)) {
            count[0] = i - start[0];
            /*fprintf(stderr, "Creating reference for %d of size %d starting at %d\n", 
                CurrentParticleType, count[0], start[0]);*/
            err = H5Sselect_hyperslab(
                ptype_dspace, H5S_SELECT_SET, start, NULL, count, NULL);
            if(err<0) ENZO_FAIL("Couldn't select hyperslab.");
            err = H5Rcreate(&reference, file_id, ptype_dset_fullname, H5R_DATASET_REGION, 
                            ptype_dspace);
            if(err<0){fprintf(stderr, "dsetname: %s\n", ptype_dset_fullname);
ENZO_FAIL("Couldn't create reference.");}
            new_dspace = H5Screate_simple(1, TempOne, NULL);
            if(err<0) ENZO_FAIL("Couldn't create new dataspace.");
            snprintf(dset_name, 22, "AddressParticleType%02d", CurrentParticleType);
            new_dset = H5Dcreate(parent_group, dset_name, H5T_STD_REF_DSETREG, 
                            new_dspace, H5P_DEFAULT);
            if(new_dset < 0) ENZO_FAIL("Couldn't create new dataset.");
            err = H5Dwrite(new_dset, H5T_STD_REF_DSETREG, H5S_ALL, H5S_ALL,
                     H5P_DEFAULT, reference);
            if(err < 0) ENZO_FAIL("Couldn't write new dataset.");
            err = H5Sclose(new_dspace);
            if(err < 0) ENZO_FAIL("Couldn't close new dataspace.");
            err = H5Dclose(new_dset);
            if(err < 0) ENZO_FAIL("Couldn't close new dataset.");

            if(i<NumberOfParticles)CurrentParticleType = ParticleType[i];
            start[0] = i;
        }
    }
    return SUCCESS;
}
