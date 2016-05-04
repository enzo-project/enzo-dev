/***********************************************************************
/
/  WRITE HDF5 STRING ATTRIBUTE
/
/  written by: Robert Harkness
/  date:       July 2002
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAILS HARD
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include <hdf5.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
void my_exit(int status);
 
// HDF5 function prototypes
 

 
 
 
 
int WriteStringAttr(hid_t dset_id, char *Alabel, char *String, FILE *log_fptr)
{
 
  hid_t       attr_id;
  hid_t       attr_dsp_id;
  hid_t       attr_type_id;
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  const char  *NoString = "none";
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
  if (io_log) fprintf(log_fptr, "Enter WSA\n");
 
//  if (io_log) fprintf(log_fptr, "  Alabel: %"ISYM" %s\n", strlen(Alabel), Alabel);
//  if (io_log) fprintf(log_fptr, "  String: %"ISYM" %s\n", strlen(String), String);
 
  attr_dsp_id = H5Screate(H5S_SCALAR);
    if (io_log) fprintf(log_fptr, "  H5Screate attr_dsp_id: %"ISYM"\n", attr_dsp_id);
    if( attr_dsp_id == h5_error ){my_exit(EXIT_FAILURE);}
 
  attr_type_id = H5Tcopy(H5T_C_S1);
                 H5Tset_size(attr_type_id, 80);
 
  attr_id = H5Acreate(dset_id, Alabel, attr_type_id,  attr_dsp_id, H5P_DEFAULT);
    if (io_log) fprintf(log_fptr, "  H5Acreate attr_id: %"ISYM"\n", attr_id);
    if( attr_id == h5_error ){my_exit(EXIT_FAILURE);}

    if (String != NULL) {
      if( strlen(String) > 0 )
      {
        h5_status = H5Awrite(attr_id, attr_type_id, (void *) String);
      }
      else
      {
        h5_status = H5Awrite(attr_id, attr_type_id, (void *) NoString);
      }
    }
 
    if (io_log) fprintf(log_fptr, "  H5Awrite: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "  H5Aclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  h5_status = H5Sclose(attr_dsp_id);
    if (io_log) fprintf(log_fptr, "  H5Sclose: %"ISYM"\n", h5_status);
    if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
 
  return SUCCESS;
 
}
