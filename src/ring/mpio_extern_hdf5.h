
// HDF5 prototypes

extern "C" herr_t H5open(void);
extern "C" herr_t H5close(void);

extern "C" hid_t H5Fcreate(const char *name, unsigned flags, hid_t create_id, hid_t access_id);
extern "C" hid_t H5Fopen(const char *name, unsigned flags, hid_t access_id);
extern "C" herr_t H5Fclose(hid_t file_id);

extern "C" hid_t H5Gcreate(hid_t loc_id, const char *name, size_t size_hint);
extern "C" hid_t H5Gopen(hid_t loc_id, const char *name);
extern "C" herr_t H5Gclose(hid_t group_id);

extern "C" hid_t H5Dcreate(hid_t loc_id, const char *name, hid_t type_id, hid_t space_id, hid_t create_plist_id);
extern "C" hid_t H5Dopen(hid_t loc_id, const char *name);
extern "C" herr_t H5Dwrite(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void *buffer);
extern "C" herr_t H5Dread(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, void *buffer);
extern "C" hid_t H5Dget_space(hid_t dataset_id);
extern "C" hid_t H5Dget_type(hid_t dataset_id);
extern "C" herr_t H5Dclose(hid_t dataset_id);

extern "C" hid_t H5Screate_simple(Eint32 rank, const hsize_t *dims, const hsize_t *maxdims);
extern "C" Eint32 H5Sget_simple_extent_ndims(hid_t space_id);
extern "C" Eint32 H5Sget_simple_extent_dims(hid_t space_id, hsize_t *dims, hsize_t *maxdims);
extern "C" herr_t H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op, const hssize_t *start, const hsize_t *stride, const hsize_t *count, const hsize_t *block);
extern "C" herr_t H5Sclose(hid_t space_id);

extern "C" hid_t H5Acreate(hid_t loc_id, const char *name, hid_t type_id, hid_t space_id, hid_t create_plist);
extern "C" hid_t H5Aopen_name(hid_t loc_id, const char *name);
extern "C" herr_t H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buffer);
extern "C" herr_t H5Aread(hid_t attr_id, hid_t mem_type_id, void *buffer);
extern "C" herr_t H5Aclose(hid_t attr_id);

extern "C" hid_t H5Tcopy(hid_t type_id);
extern "C" herr_t H5Tset_size(hid_t type_id, size_t size);

extern "C" hid_t H5Pcreate (hid_t type);
extern "C" herr_t H5Pset_fapl_mpio(hid_t fapl_id, MPI_Comm comm, MPI_Info info);
extern "C" herr_t H5Pget_fapl_mpio(hid_t fapl_id, MPI_Comm *comm, MPI_Info *info);
extern "C" herr_t H5Pset_dxpl_mpio(hid_t dxpl_id, H5FD_mpio_xfer_t xfer_mode);
extern "C" herr_t H5Pget_dxpl_mpio(hid_t dxpl_id, H5FD_mpio_xfer_t *xfer_mode);
extern "C" herr_t H5Pclose(hid_t plist);

