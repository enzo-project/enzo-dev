/***********************************************************************
/
/  READ HDF5 DATASET ATTRIBUTE
/
/  written by: Robert Harkness
/  date:       July 2002
/  modified1:  Robert Harkness
/  date:       May, 2008
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAILS HARD
/
************************************************************************/
 
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
void my_exit(int status);
 
// HDF5 function prototypes
 

 
 
 
 
int ReadAttr(char *Fname, int *Rank, int Dims[], int *NSeg, int *LSeg, FILE *log_fptr)
{
 
  hid_t       file_id, dset_id, attr_id;
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int         component_rank_attr;
  int         component_size_attr;
  int         field_rank_attr;
  int         field_dims_attr[3];
 
  int         dim;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
  if (io_log) fprintf(log_fptr, "H5Fopen with Name = %s\n", Fname);
 
  file_id = H5Fopen(Fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (io_log) fprintf(log_fptr, "H5Fopen id: %"ISYM"\n", file_id);
    if( file_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  if (io_log) fprintf(log_fptr, "H5Dopen with Name = %s\n", Fname);
 
  dset_id =  H5Dopen(file_id, Fname);
    if (io_log) fprintf(log_fptr, "H5Dopen id: %"ISYM"\n", dset_id);
    if( dset_id == h5_error ){my_exit(EXIT_FAILURE);}
 
 
  if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Rank\n");
 
  attr_id = H5Aopen_name(dset_id, "Component_Rank");
    if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aread(attr_id, HDF5_INT, &component_rank_attr);
    if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  if (io_log) fprintf(log_fptr, "COMPONENT RANK %"ISYM"\n", component_rank_attr);
 
 
  if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Component_Size\n");
 
  attr_id = H5Aopen_name(dset_id, "Component_Size");
    if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aread(attr_id, HDF5_INT, &component_size_attr);
    if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  if (io_log) fprintf(log_fptr, "COMPONENT SIZE %"ISYM"\n", component_size_attr);
 
 
  if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Rank\n");
 
  attr_id = H5Aopen_name(dset_id, "Rank");
    if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aread(attr_id, HDF5_INT, &field_rank_attr);
    if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  if (io_log) fprintf(log_fptr, "RANK %"ISYM"\n", field_rank_attr);
 
 
  if (io_log) fprintf(log_fptr, "H5Aopen_name with Name = Dimensions\n");
 
  attr_id = H5Aopen_name(dset_id, "Dimensions");
    if (io_log) fprintf(log_fptr, "H5Aopen_name id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aread(attr_id, HDF5_INT, field_dims_attr);
    if (io_log) fprintf(log_fptr, "H5Aread: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  for ( dim=0; dim < field_rank_attr; dim++ )
  {
    if (io_log) fprintf(log_fptr, "DIMS %"ISYM":  %"ISYM"\n", dim, field_dims_attr[dim]);
  }
 
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log_fptr, "H5Dclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log_fptr, "H5Fclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  *Rank = field_rank_attr;
  *NSeg = component_rank_attr;
  *LSeg = component_size_attr;
 
  for(dim = 0; dim < field_rank_attr; dim++)
  {
    Dims[dim] = field_dims_attr[dim];
  }
 
  return SUCCESS;
 
}
