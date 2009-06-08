/***********************************************************************
/
/  OUTPUT INT FIELD TO A HDF5 FILE - CHUNKS FOR LARGE ICs
/
/  written by: Robert Harkness
/  date:       August 2005
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
 
#include "macros_and_parameters.h"
 
// HDF5 function prototypes
 
#include "extern_hdf5.h"
 
 
 
 
int WriteIntField(int Rank, int Dims[3], int *Field, char *Name, int Part, int Npart,
               int GridRank, int Starts[3], int Ends[3], int Tops[3])
{
 
 
  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id;
  hid_t       int_file_type_id, int_mem_type_id;
 
  hsize_t     out_dims[3];
  hsize_t     dimm;
  hsize_t     slab_dims[4];
  hsize_t     slab_rank;
  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[4], file_count[4], file_block[4];
  hsize_t     attr_count;

  hsize_t     bsize;
 
  hssize_t    mem_offset;
  hssize_t    file_offset[4];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int         dim;
 
  int         component_rank_attr;
  int         component_size_attr;
  int         field_rank_attr;
  int         field_dims_attr[3];

  int         chunk;
  int         numchunks;
 
  FILE        *dumpfile;
  FILE        *log;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
#ifdef DUMP_OK
  int         dump_ok = 1;
#else
  int         dump_ok = 0;
#endif
 
  if (dump_ok) dumpfile = fopen("DumpWF","a");
  if (io_log) log = fopen("IO_Log","a");
 
  if (io_log) fprintf(log, "On entry to WriteIntField\n");
  if (io_log) fprintf(log, "  Rank %"ISYM"\n", Rank);
  if (io_log) fprintf(log, "  Dims %"ISYM"  %"ISYM"  %"ISYM"\n", Dims[0], Dims[1], Dims[2]);
  if (io_log) fprintf(log, "  Name %s\n", Name);
  if (io_log) fprintf(log, "  Part %"ISYM" of %"ISYM"\n", Part, Npart);
 
//  GB: Reverse dim ordering since we are using fortran array ordering
//      (actually, we don't have to do this here).
//
//  out_dims[Rank-dim-1] = Dims[dim];
 
  for ( dim = 0; dim < Rank; dim++ )
  {
    out_dims[dim] = Dims[dim];
  }
 
//  reverse these also?
 
  for ( dim =0; dim < GridRank; dim++ )
  {
    if ( Starts[dim] == INT_UNDEFINED )
      Starts[dim] = 0;
    if ( Ends[dim] == INT_UNDEFINED )
      Ends[dim] = Tops[dim] - 1;
  }
 
  slab_rank = Rank+1;
 
  slab_dims[0] = Npart;
 
  for ( dim = 1; dim < slab_rank; dim++ )
  {
    slab_dims[dim] = out_dims[dim-1];
  }
 
  if (io_log) fprintf(log, "  Extended Rank %"ISYM"\n", (int) slab_rank);
 
  for ( dim = 0; dim < slab_rank; dim++ )
  {
    if (io_log) fprintf(log, "    %"ISYM":  %"ISYM"\n", dim, (int) slab_dims[dim]);
  }
 
  dimm = 1;
 
  for ( dim = 0; dim < Rank; dim++ )
  {
    dimm = dimm * Dims[dim];
  }
 
  if (io_log) fprintf(log, "  Grid Elements %"ISYM"\n", (int) dimm);
 
  component_rank_attr = Npart;
  component_size_attr = dimm;
 
  field_rank_attr = Rank;
 
  for ( dim = 0; dim < Rank; dim++ )
  {
    field_dims_attr[dim] = Dims[dim];
  }
 
  /* HDF5
 
     The HDF4 DFSD interface does not have a direct analogue in HDF5.
     In particular, there is no "append" on open and there are no
     specific features to store rank, dimensions, scales etc.
     These are replaced by dataset attributes, which can be named.
     To eliminate redundancy and to keep the data structure close
     to the original, each ENZO file will contain a single dataset
     of the same name (as opposed to separate datasets for each
     component).
 
     The dataspace uses the number of components for each of the
     the dimensions of the dataset as an additional "dimension".
     For example, a 3D scalar field of 4x4x4 points will have a
     dataspace of {1,4,4,4}, while a 3D vector field of, say,
     {Vx,Vy,Vz} and 4x4x4 points will have a dataspace of {3,4,4,4}.
 
     Slab I/O is used in preparation for parallel I/O.
 
     It is ASSUMED that this routine is called with
     Part = {0,1,2,...,Npart-1} to write Npart components of equal length
     (the product of the field dimensions {Dims[i] for i=0,Rank-1}.
     Each component is offset in the file by {0,1,2,...,Npart-1} * the field size.
 
     The rank and dimensions of the field and the number of field
     components are HDF5 attributes.
 
  */
 
  int ll = sizeof(Eint);
 
  switch(ll)
  {
 
    case 4:
      int_mem_type_id = HDF5_I4;
      int_file_type_id = HDF5_FILE_I4;
      break;
    case 8:
      int_mem_type_id = HDF5_I8;
      int_file_type_id = HDF5_FILE_I8;
      break;
    default:
      int_mem_type_id = HDF5_I8;
      int_file_type_id = HDF5_FILE_I8;
  }
 
 
  int ii = sizeof(int);
 
  switch(ii)
  {
 
    case 4:
      mem_type_id = HDF5_I8;
      file_type_id = HDF5_FILE_I8;
      break;
 
    case 8:
      mem_type_id = HDF5_I8;
      file_type_id = HDF5_FILE_I8;
      break;
 
    default:
      mem_type_id = HDF5_I8;
      file_type_id = HDF5_FILE_I8;
 
  }

// If Rank = 3, chunk in planes of Y*Z

  if( Rank == 3 ) {
    numchunks = out_dims[0];
    bsize = out_dims[1]*out_dims[2];
    fprintf(stderr, "3D Chunk by planes of Y*Z\n");
  }

// If Rank = 2, what?

  if( Rank == 2 ) {
    fprintf(stderr, "Fix this for 2D!\n");
  }

// If Rank = 1, chunk according to size

  if( Rank == 1 ) {
    numchunks = nint(POW( ((double) dimm) , (1.0/3.0) ));
    bsize = dimm / numchunks;
    if ( dimm % bsize > 0 ) numchunks = numchunks+1;
    fprintf(stderr, "1D Chunk by N^(1/3)\n");
  }

  fprintf(stderr, "Numchunks %d\n", numchunks);
  fprintf(stderr, "Bsize %d\n", bsize);

  mem_dsp_id = H5Screate_simple(1, &bsize, NULL);
    if (io_log) fprintf(log, "H5Screate mem_dsp_id: %"ISYM"\n", mem_dsp_id);
    assert( mem_dsp_id != h5_error );

  file_dsp_id = H5Screate_simple(((Eint32) slab_rank), slab_dims, NULL);
    if (io_log) fprintf(log, "H5Screate file_dsp_id: %"ISYM"\n", file_dsp_id);
    assert( file_dsp_id != h5_error );

//  If Part is zero, create an HDF5 file with a single dataset of
//  the same name and attach the dataset attributes, otherwise
//  open an existing file and dataset.
 
  if ( Part == 0 )
  {
    if (io_log) fprintf(log, "Calling H5Fcreate with Name = %s\n", Name);
 
    file_id = H5Fcreate(Name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fcreate id: %"ISYM"\n", file_id);
      assert( file_id != h5_error );
 
    if (io_log) fprintf(log, "Calling H5Dcreate with Name = %s\n", Name);
 
    dset_id =  H5Dcreate(file_id, Name, int_file_type_id, file_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Dcreate id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );
 
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "Component_Rank",  int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id,  int_mem_type_id, &component_rank_attr);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "Component_Size", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id,  int_mem_type_id, &component_size_attr);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
 
    attr_count = 1;
 
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "Rank", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id,  int_mem_type_id, &field_rank_attr);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
 
    attr_count = Rank;
 
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "Dimensions", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id,  int_mem_type_id, field_dims_attr);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
 
    attr_count = GridRank;
 
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "TopGridStart", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id,  int_mem_type_id, Starts);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_count = GridRank;
 
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "TopGridEnd", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id,  int_mem_type_id, Ends);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    attr_count = GridRank;
 
    attr_dsp_id = H5Screate_simple(1, &attr_count, NULL);
      if (io_log) fprintf(log, "H5Screate_simple: %"ISYM"\n", attr_dsp_id);
      assert( attr_dsp_id != h5_error );
 
    attr_id = H5Acreate(dset_id, "TopGridDims", int_file_type_id, attr_dsp_id, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Acreate: %"ISYM"\n", attr_id);
      assert( attr_id != h5_error );
 
    h5_status = H5Awrite(attr_id,  int_mem_type_id, Tops);
      if (io_log) fprintf(log, "H5Awrite: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Aclose(attr_id);
      if (io_log) fprintf(log, "H5Aclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
 
    h5_status = H5Sclose(attr_dsp_id);
      if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
      assert( h5_status != h5_error );
  }
 
  else
 
  {
    if (io_log) fprintf(log, "Calling H5Fopen with Name = %s\n", Name);
 
    file_id = H5Fopen(Name, H5F_ACC_RDWR, H5P_DEFAULT);
      if (io_log) fprintf(log, "H5Fopen id: %"ISYM"\n", file_id);
      assert( file_id != h5_error );
 
    if (io_log) fprintf(log, "Calling H5Dopen with Name = %s\n", Name);
 
    dset_id =  H5Dopen(file_id, Name);
      if (io_log) fprintf(log, "H5Dopen id: %"ISYM"\n", dset_id);
      assert( dset_id != h5_error );
  }

// FOR-- 

// numchunks = out_dims[0];

for( chunk = 0; chunk < numchunks; chunk++) {

//  fprintf(stderr, "CHUNK %d\n", chunk);

// Data in memory is considered 1D, stride 1, with zero offset

  mem_stride = 1;      // contiguous elements
  mem_count = bsize;   // number of elements in chunk
  mem_offset = 0;      // zero offset in buffer
  mem_block = 1;       // single element blocks

  if( (dimm % bsize) != 0 ) bsize = (dimm % bsize); 

// 1D memory model
 
  h5_status =  H5Sselect_hyperslab(mem_dsp_id,  H5S_SELECT_SET, &mem_offset, &mem_stride, &mem_count, NULL);
    if (io_log) fprintf(log, "H5Sselect mem slab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 

//  Data in the file is (1+Rank)D with Npart components per grid point.
//  Offset[0] is the component Part of Npart components.  Data for each
//  Part are contiguous in the file, so stride = 1.

  if( Rank == 1) {
  file_stride[0] = 1;      // contiguous elements
  file_count[0] = 1;       // one component per call
  file_offset[0] = Part;   // component Part of Npart
  file_block[0] = 1;       // single element blocks

  file_stride[1] = 1;                   // contiguous elements
  file_count[1] = bsize;                // field dimensions
  file_offset[1] = chunk*bsize;               // complete field, no offset
  file_block[1] = 1;                    // single element blocks
  }

  if( Rank > 1 ) {
  file_stride[0] = 1;      // contiguous elements
  file_count[0] = 1;       // one component per call
  file_offset[0] = Part;   // component Part of Npart
  file_block[0] = 1;       // single element blocks

  file_stride[1] = 1;                   // contiguous elements
  file_count[1] = 1;                    // field dimensions
  file_offset[1] = chunk;               // complete field, no offset
  file_block[1] = 1;                    // single element blocks

  file_stride[2] = 1;                   // contiguous elements
  file_count[2] = out_dims[1];          // field dimensions
  file_offset[2] = 0;                   // complete field, no offset
  file_block[2] = 1;                    // single element blocks

  file_stride[3] = 1;                   // contiguous elements
  file_count[3] = out_dims[2];          // field dimensions
  file_offset[3] = 0;                   // complete field, no offset
  file_block[3] = 1;                    // single element blocks
  }

  h5_status = H5Sselect_hyperslab(file_dsp_id, H5S_SELECT_SET, file_offset, file_stride, file_count, NULL);
    if (io_log) fprintf(log, "H5Sselect file slab: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
  h5_status = H5Dwrite(dset_id, int_mem_type_id, mem_dsp_id, file_dsp_id,  H5P_DEFAULT, &Field[chunk*bsize]);
    if (io_log) fprintf(log, "H5Dwrite: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );

  }
// FOR end
 
  h5_status = H5Dclose(dset_id);
    if (io_log) fprintf(log, "H5Dclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(mem_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Sclose(file_dsp_id);
    if (io_log) fprintf(log, "H5Sclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  h5_status = H5Fclose(file_id);
    if (io_log) fprintf(log, "H5Fclose: %"ISYM"\n", h5_status);
    assert( h5_status != h5_error );
 
  if (io_log) fprintf(log, "Exit WriteIntField\n");
 
  if (dump_ok) fclose(dumpfile);
  if (io_log) fclose(log);
 
  return SUCCESS;
}
