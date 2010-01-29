module hdf_io

#ifdef ENZO
  use hdf5
  use grafic_types
  implicit none

  integer(HID_T), private,save :: plist_id, file_id, dset_id
  integer, private, dimension(3) :: topgrid_dims

  interface hdf_write_attr
     module procedure hdf_write_attr_int, hdf_write_attr_real, &
          hdf_write_attr_int_scalar, hdf_write_attr_real_scalar
  end interface

  interface hdf_read_attr
     module procedure hdf_read_attr_int, hdf_read_attr_real, &
          hdf_read_attr_int_scalar, hdf_read_attr_real_scalar
  end interface

  interface write_phdf5_data
     module procedure write_phdf5_data_single, write_phdf5_data_double
  end interface

  interface read_phdf5_data
     module procedure read_phdf5_data_single, read_phdf5_data_double
  end interface
#endif

contains

#ifdef ENZO
  subroutine init_phdf5(filename, overwrite_in, serial_in)
    character(128), intent(in) :: filename
    logical, optional :: overwrite_in  ! read (1) or write (0)
    logical, optional :: serial_in
    logical :: overwrite, serial
    integer :: mdc_nelmts
    integer(SIZE_T) :: rdcc_nelmts, rdcc_nbytes
    integer(SIZE_T) :: new_cache_size, new_buf_size
    real :: rdcc_w0
    
    include 'mpif.h'
    integer ierr

    overwrite = .false.
    serial = .false.
    if (present(overwrite_in)) overwrite = overwrite_in
    if (present(serial_in)) serial = serial_in

    new_cache_size = 16*1024*1024
    new_buf_size = 1*1024*1024

    ! Open HDF5 file for parallel I/O
    CALL h5open_f(ierr)
    if (serial) then
       plist_id = H5P_DEFAULT_F
    else
       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
       CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)
       CALL h5pset_sieve_buf_size_f(plist_id, new_buf_size, ierr)
       CALL h5pget_cache_f(plist_id, mdc_nelmts, rdcc_nelmts, rdcc_nbytes, &
            rdcc_w0, ierr)
       rdcc_nbytes = new_cache_size
       CALL h5pset_cache_f(plist_id, mdc_nelmts, rdcc_nelmts, rdcc_nbytes, &
            rdcc_w0, ierr)
    endif
    
    ! Create/open the file
    if (overwrite) then
       CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr, &
            access_prp = plist_id)
    else
       CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, ierr, &
            access_prp = plist_id)
    endif
    CALL h5pclose_f(plist_id, ierr)

    ! Create properties list for collective dataset access
    if (.not.serial) then
       CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
       CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, ierr)
    endif
  end subroutine init_phdf5

  subroutine close_phdf5_file
    integer ierr
    CALL h5dclose_f(dset_id, ierr)
    CALL h5pclose_f(plist_id, ierr)
    CALL h5fclose_f(file_id, ierr)
    CALL h5close_f(ierr)
  end subroutine close_phdf5_file

  !--------------------------------------------------------------
  !---------------------- HEADER ROUTINES -----------------------
  !--------------------------------------------------------------
  
  subroutine grafic_hdf5_read_header_white(filename, nx, ny, nz, iseed)
    character(128), intent(in) :: filename
    integer, intent(out) :: nx, ny, nz, iseed
    integer :: ierr, dims(3)
    
    CALL h5open_f(ierr)
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)
    CALL hdf_read_attr(file_id, "Dimensions", dims, 3)
    CALL hdf_read_attr(file_id, "Random Seed", iseed)
    CALL h5fclose_f(file_id, ierr)
    CALL h5close_f(ierr)
    nx = dims(1)
    ny = dims(2)
    nz = dims(3)

  end subroutine grafic_hdf5_read_header_white

  subroutine grafic_hdf5_write_header_white(filename, nx, ny, nz, iseed)
    character(128), intent(in) :: filename
    integer, intent(in) :: nx, ny, nz, iseed
    integer :: ierr, dims(3)

    dims = (/nx, ny, nz/)

    CALL h5open_f(ierr)
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)
    CALL hdf_write_attr(file_id, "Dimensions", dims, 3)
    CALL hdf_write_attr(file_id, "Random Seed", iseed)
    CALL h5fclose_f(file_id, ierr)
    CALL h5close_f(ierr)
  end subroutine grafic_hdf5_write_header_white

  subroutine grafic_hdf5_read_header(filename, params, hcosmo)
    character(128), intent(in) :: filename
    type(taille), intent(out) :: params
    type(cosmo), intent(out) :: hcosmo
    integer i, ierr, dims(3)
    real offset(3)

    CALL h5open_f(ierr)
    CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, ierr)

    CALL hdf_read_attr(file_id, "Dimensions", dims, 3)
    CALL hdf_read_attr(file_id, "dx", params%dx)
    CALL hdf_read_attr(file_id, "Offset", offset, 3)
    CALL hdf_read_attr(file_id, "a_start", hcosmo%astart)
    CALL hdf_read_attr(file_id, "omega_b", hcosmo%omegab)
    CALL hdf_read_attr(file_id, "omega_m", hcosmo%omegam)
    CALL hdf_read_attr(file_id, "omega_v", hcosmo%omegav)
    CALL hdf_read_attr(file_id, "h0", hcosmo%h0)
    CALL hdf_read_attr(file_id, "vfact", hcosmo%vfact)

    CALL h5fclose_f(file_id, ierr)
    CALL h5close_f(ierr)

    params%nx = dims(1)
    params%ny = dims(2)
    params%nz = dims(3)
    params%lx = offset(1)
    params%ly = offset(2)
    params%lz = offset(3)

  end subroutine grafic_hdf5_read_header

  subroutine grafic_hdf5_write_header(filename, params, hcosmo)
    character(128), intent(in) :: filename
    type(taille), intent(in) :: params
    type(cosmo), intent(in) :: hcosmo
    integer ierr, dims(3)
    real offset(3)

    dims = (/params%nx, params%ny, params%nz/)
    offset = (/params%lx, params%ly, params%lz/)

    CALL h5open_f(ierr)
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, ierr)

    CALL hdf_write_attr(file_id, "Dimensions", dims, 3)
    CALL hdf_write_attr(file_id, "dx", params%dx)
    CALL hdf_write_attr(file_id, "Offset", offset, 3)
    CALL hdf_write_attr(file_id, "a_start", hcosmo%astart)
    CALL hdf_write_attr(file_id, "omega_b", hcosmo%omegab)
    CALL hdf_write_attr(file_id, "omega_m", hcosmo%omegam)
    CALL hdf_write_attr(file_id, "omega_v", hcosmo%omegav)
    CALL hdf_write_attr(file_id, "h0", hcosmo%h0)
    CALL hdf_write_attr(file_id, "vfact", hcosmo%vfact)

    CALL h5fclose_f(file_id, ierr)
    CALL h5close_f(ierr)

  end subroutine grafic_hdf5_write_header

  subroutine set_topgrid_dims(dims)
    integer, intent(in), dimension(3) :: dims
    topgrid_dims = dims
  end subroutine set_topgrid_dims

  !---------------------------------------------------------------
  !---------------------- DATASET ROUTINES -----------------------
  !---------------------------------------------------------------

  subroutine create_hdf5_dset(name, datatype, rank, dimsf, filepart, n_filepart, &
       level)
    character(128), intent(in) :: name
    integer(HID_T), intent(in) :: datatype
    integer, intent(in) :: rank, filepart, n_filepart
    integer, intent(in) :: dimsf(rank)
    integer, optional :: level
    integer(HID_T) :: filespace
    integer(HSIZE_T), allocatable :: file_dims(:)
    integer :: file_rank, ierr, datasize, this_level
    integer, dimension(3) :: top_start, top_end, top_dims
    
    ! Enzo expects all velocity components in one file.  So the
    ! filespace dimensions are (4,nx,ny,nz) for 3-rank vectors.
    file_rank = rank+1
    allocate(file_dims(file_rank))
    file_dims(file_rank) = n_filepart
    file_dims(1:rank) = dimsf(1:rank)

    datasize = product(dimsf)
    CALL h5screate_simple_f(file_rank, file_dims, filespace, ierr)
    CALL h5dcreate_f(file_id, name, datatype, filespace, dset_id, ierr)
    CALL h5sclose_f(filespace, ierr)

    CALL hdf_write_attr(dset_id, "Component_Rank", n_filepart)
    CALL hdf_write_attr(dset_id, "Component_Size", datasize)
    CALL hdf_write_attr(dset_id, "Dimensions", dimsf, 3)
    CALL hdf_write_attr(dset_id, "Rank", rank)

    ! In the public version of Enzo, the nested grid initalization
    ! routine expects the attributes TopGridDims, TopGridStart,
    ! TopGridEnd to correctly place the subgrids (and to check against
    ! parameters).  Having the subgrids centered at 0.5 makes life
    ! easier because we only need the grid level and dimensions to
    ! compute the start and end coordinates of the grid.
    if (present(level)) then
       top_dims = dimsf / 2**level
       top_start = topgrid_dims/2 - top_dims/2
       top_end = topgrid_dims/2 + top_dims/2
       CALL hdf_write_attr(dset_id, "TopGridStart", top_start, 3)
       CALL hdf_write_attr(dset_id, "TopGridEnd", top_end, 3)
       CALL hdf_write_attr(dset_id, "TopGridDims", topgrid_dims, 3)
    endif

    deallocate(file_dims)

  end subroutine create_hdf5_dset

  subroutine open_hdf5_dset(name)
    character(128), intent(in) :: name
    integer ierr
    CALL h5dopen_f(file_id, name, dset_id, ierr)
  end subroutine open_hdf5_dset

  !------------------------- WRITING -------------------------

  subroutine write_phdf5_data_single(buffer, zstart, local_nz, nz, ny, nx, &
       filepart, n_filepart, masked_in)
    real(sp), dimension(:), intent(in) :: buffer
    integer, intent(in) :: filepart, n_filepart
    integer, intent(in) :: nx, ny, nz
    integer, intent(in) :: zstart, local_nz
    logical, optional, intent(in) :: masked_in
        
    integer(HID_T) memspace, filespace
    integer(4) :: mem_rank, file_rank
    integer(HSIZE_T), dimension(:), allocatable :: mem_dims, mem_offset, &
         mem_odims, dims, odims, offset
    integer :: ierr, n2x
    logical :: masked

    masked = .false.
    if (present(masked_in)) masked = masked_in

    mem_rank = 3
    if (nx.le.1 .and. ny.le.1) mem_rank = 1
    file_rank = mem_rank + 1

    allocate(mem_dims(mem_rank))
    allocate(mem_offset(mem_rank))
    allocate(mem_odims(mem_rank))
    allocate(dims(file_rank))
    allocate(odims(file_rank))
    allocate(offset(file_rank))

    if (mem_rank.eq.3) then
       if (masked) then 
          n2x = nx
       else
          n2x = 2*(nx/2+1)
       endif
       dims = (/nx, ny, local_nz, 1/)
       odims = (/nx, ny, nz, n_filepart/)
       offset = (/0, 0, zstart, filepart/)
       mem_odims = (/n2x, ny, max(local_nz,1)/)  ! accounts for padding
       mem_dims = (/nx, ny, local_nz/)
       mem_offset = (/0, 0, 0/)
    elseif (mem_rank.eq.1) then
       dims = (/local_nz, 1/)
       odims = (/nz, n_filepart/)
       offset = (/zstart, filepart/)
       mem_odims = (/local_nz/)
       mem_dims = (/local_nz/)
       mem_offset = (/0/)
    endif

    ! Create memory- (1D) and file-space
    CALL h5screate_simple_f(mem_rank, mem_odims, memspace, ierr)
    CALL h5dget_space_f(dset_id, filespace, ierr)

    ! Select the hyperslab we're going to read/write to
    if (local_nz > 0) then
       CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, mem_offset, &
            mem_dims, ierr)
       CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, ierr)
    else
       CALL h5sselect_none_f(memspace, ierr)
       CALL h5sselect_none_f(filespace, ierr)
    end if

    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, buffer, odims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)

    deallocate(mem_dims, mem_offset, mem_odims)
    deallocate(dims, odims, offset)

  end subroutine write_phdf5_data_single

  subroutine write_phdf5_data_double(buffer, zstart, local_nz, nz, ny, nx, &
       filepart, n_filepart, masked_in)
    real(dp), dimension(:), intent(in) :: buffer
    integer, intent(in) :: filepart, n_filepart
    integer, intent(in) :: nx, ny, nz
    integer, intent(in) :: zstart, local_nz
    logical, optional, intent(in) :: masked_in
        
    integer(HID_T) memspace, filespace
    integer(4) :: mem_rank, file_rank
    integer(HSIZE_T), dimension(:), allocatable :: mem_dims, mem_offset, &
         mem_odims, dims, odims, offset
    integer :: ierr, n2x
    logical :: masked

    masked = .false.
    if (present(masked_in)) masked = masked_in

    mem_rank = 3
    if (nx.le.1 .and. ny.le.1) mem_rank = 1
    file_rank = mem_rank + 1

    allocate(mem_dims(mem_rank))
    allocate(mem_offset(mem_rank))
    allocate(mem_odims(mem_rank))
    allocate(dims(file_rank))
    allocate(odims(file_rank))
    allocate(offset(file_rank))

    if (mem_rank.eq.3) then
       if (masked) then 
          n2x = nx
       else
          n2x = 2*(nx/2+1)
       endif
       dims = (/nx, ny, local_nz, 1/)
       odims = (/nx, ny, nz, n_filepart/)
       offset = (/0, 0, zstart, filepart/)
       mem_odims = (/n2x, ny, max(local_nz,1)/)  ! accounts for padding
       mem_dims = (/nx, ny, local_nz/)
       mem_offset = (/0, 0, 0/)
    elseif (mem_rank.eq.1) then
       dims = (/local_nz, 1/)
       odims = (/nz, n_filepart/)
       offset = (/zstart, filepart/)
       mem_odims = (/local_nz/)
       mem_dims = (/local_nz/)
       mem_offset = (/0/)
    endif

    ! Create memory- (1D) and file-space
    CALL h5screate_simple_f(mem_rank, mem_odims, memspace, ierr)
    CALL h5dget_space_f(dset_id, filespace, ierr)

    ! Select the hyperslab we're going to read/write to
    if (local_nz > 0) then
       CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, mem_offset, &
            mem_dims, ierr)
       CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, ierr)
    else
       CALL h5sselect_none_f(memspace, ierr)
       CALL h5sselect_none_f(filespace, ierr)
    end if

    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buffer, odims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)

    deallocate(mem_dims, mem_offset, mem_odims)
    deallocate(dims, odims, offset)

  end subroutine write_phdf5_data_double

  !------------------------- READING -------------------------

  subroutine read_phdf5_data_single(buffer, zstart, local_nz, ny, nx, filepart, &
       padding_in)
    real(sp), dimension(:), intent(inout) :: buffer
    integer, intent(in) :: zstart, local_nz, ny, nx, filepart
    logical, optional, intent(in) :: padding_in

    integer(HID_T) memspace, filespace
    integer(4) :: mem_rank, file_rank
    integer(HSIZE_T), dimension(:), allocatable :: mem_dims, mem_offset, &
         mem_odims, dims, odims, offset
    integer :: ierr, n2x
    logical :: padding

    padding = .false.
    if (present(padding_in)) padding = padding_in

    mem_rank = 3
    if (nx.le.1 .and. ny.le.1) mem_rank = 1
    file_rank = mem_rank + 1

    allocate(mem_dims(mem_rank))
    allocate(mem_offset(mem_rank))
    allocate(mem_odims(mem_rank))
    allocate(dims(file_rank))
    allocate(odims(file_rank))
    allocate(offset(file_rank))

    if (mem_rank.eq.3) then
       if (padding) then
          n2x = nx
       else
          n2x = 2*(nx/2+1)
       endif
       dims = (/nx, ny, local_nz, 1/)
       offset = (/0, 0, zstart, filepart/)
       mem_odims = (/n2x, ny, local_nz/)  ! accounts for padding
       mem_dims = (/nx, ny, local_nz/)
       mem_offset = (/0, 0, 0/)
    elseif (mem_rank.eq.1) then
       dims = (/local_nz, 1/)
       offset = (/zstart, filepart/)
       mem_odims = (/local_nz/)
       mem_dims = (/local_nz/)
       mem_offset = (/0/)
    endif

    ! Create a file dataspace independently, then select the hyperslab
    CALL h5dget_space_f(dset_id, filespace, ierr)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, ierr)

    ! Now create a memory dataspace
    CALL h5screate_simple_f(mem_rank, mem_odims, memspace, ierr)
    CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, mem_offset, mem_dims, ierr)

    ! Finally read the data collectively
    CALL h5dread_f(dset_id, H5T_NATIVE_REAL, buffer, dims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)

    deallocate(mem_dims, mem_offset, mem_odims)
    deallocate(dims, odims, offset)
    
  end subroutine read_phdf5_data_single

  subroutine read_phdf5_data_double(buffer, zstart, local_nz, ny, nx, filepart, &
       padding_in)
    real(dp), dimension(:), intent(inout) :: buffer
    integer, intent(in) :: zstart, local_nz, ny, nx, filepart
    logical, optional, intent(in) :: padding_in
    
    integer(HID_T) memspace, filespace
    integer(4) :: mem_rank, file_rank
    integer(HSIZE_T), dimension(:), allocatable :: mem_dims, mem_offset, &
         mem_odims, dims, odims, offset
    integer :: ierr, n2x
    logical :: padding

    padding = .false.
    if (present(padding_in)) padding = padding_in

    mem_rank = 3
    if (nx.le.1 .and. ny.le.1) mem_rank = 1
    file_rank = mem_rank + 1

    allocate(mem_dims(mem_rank))
    allocate(mem_offset(mem_rank))
    allocate(mem_odims(mem_rank))
    allocate(dims(file_rank))
    allocate(odims(file_rank))
    allocate(offset(file_rank))

    if (mem_rank.eq.3) then
       if (padding) then
          n2x = nx
       else
          n2x = 2*(nx/2+1)
       endif
       dims = (/nx, ny, local_nz, 1/)
       offset = (/0, 0, zstart, filepart/)
       mem_odims = (/n2x, ny, local_nz/)  ! accounts for padding
       mem_dims = (/nx, ny, local_nz/)
       mem_offset = (/0, 0, 0/)
    elseif (mem_rank.eq.1) then
       dims = (/local_nz, 1/)
       offset = (/zstart, filepart/)
       mem_odims = (/local_nz/)
       mem_dims = (/local_nz/)
       mem_offset = (/0/)
    endif

    ! Create a file dataspace independently, then select the hyperslab
    CALL h5dget_space_f(dset_id, filespace, ierr)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, ierr)

    ! Now create a memory dataspace
    CALL h5screate_simple_f(mem_rank, mem_odims, memspace, ierr)
    CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, mem_offset, mem_dims, ierr)

    ! Finally read the data collectively
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buffer, dims, ierr, &
         file_space_id = filespace, mem_space_id = memspace, &
         xfer_prp = plist_id)

    CALL h5sclose_f(filespace, ierr)
    CALL h5sclose_f(memspace, ierr)

    deallocate(mem_dims, mem_offset, mem_odims)
    deallocate(dims, odims, offset)
    
  end subroutine read_phdf5_data_double

  !----------------------------------------------------------------
  !--------------------- ATTRIBUTE ROUTINES -----------------------
  !----------------------------------------------------------------
  
  ! integer type
  subroutine hdf_write_attr_int(group_id, name, value, n)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in) :: name
    integer, intent(in) :: n
    integer, intent(in) :: value(n)

    integer(HID_T) :: aspace_id, attr_id
    integer(HSIZE_T) :: n8(1)
    integer ierr
    integer(4) :: rank

    n8(1) = n
    rank = 1
    CALL h5screate_simple_f(rank, n8, aspace_id, ierr)
    CALL h5acreate_f(group_id, name, H5T_NATIVE_INTEGER, aspace_id, attr_id, ierr)
    CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, n8, ierr)
    CALL h5aclose_f(attr_id, ierr)
    CALL h5sclose_f(aspace_id, ierr)
  end subroutine hdf_write_attr_int

  subroutine hdf_write_attr_int_scalar(group_id, name, value)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in) :: name
    integer, intent(in) :: value

    integer(HID_T) :: aspace_id, attr_id
    integer(HSIZE_T) :: n8(1) = (/1/)
    integer ierr
    integer(4) :: rank = 1

    CALL h5screate_simple_f(rank, n8, aspace_id, ierr)
    CALL h5acreate_f(group_id, name, H5T_NATIVE_INTEGER, aspace_id, attr_id, ierr)
    CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, (/value/), n8, ierr)
    CALL h5aclose_f(attr_id, ierr)
    CALL h5sclose_f(aspace_id, ierr)
  end subroutine hdf_write_attr_int_scalar

  subroutine hdf_read_attr_int(group_id, name, value, n)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in) :: name
    integer, intent(in) :: n
    integer, intent(out) :: value(n)

    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: nn
    integer ierr

    nn(1) = n
    CALL h5aopen_name_f(group_id, name, attr_id, ierr)
    CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, value, nn, ierr)
    CALL h5aclose_f(attr_id, ierr)
  end subroutine hdf_read_attr_int

  subroutine hdf_read_attr_int_scalar(group_id, name, value)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in) :: name
    integer, intent(out) :: value

    integer(HID_T) :: attr_id
    integer(HSIZE_T) :: n(1)
    integer ierr

    n(1) = 1
    CALL h5aopen_name_f(group_id, name, attr_id, ierr)
    CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, value, n, ierr)
    CALL h5aclose_f(attr_id, ierr)
  end subroutine hdf_read_attr_int_scalar

  ! real type
  subroutine hdf_write_attr_real(group_id, name, value, n)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in) :: name
    integer, intent(in) :: n
    real, intent(in) :: value(n)

    integer(HID_T) :: aspace_id, attr_id
    integer(HSIZE_T) :: n8(1)
    integer ierr
    integer(4) :: rank

    n8(1) = n
    rank = 1
    CALL h5screate_simple_f(rank, n8, aspace_id, ierr)
    CALL h5acreate_f(group_id, name, H5T_NATIVE_REAL, aspace_id, attr_id, ierr)
    CALL h5awrite_f(attr_id, H5T_NATIVE_REAL, value, n8, ierr)
    CALL h5aclose_f(attr_id, ierr)
    CALL h5sclose_f(aspace_id, ierr)
  end subroutine hdf_write_attr_real

  subroutine hdf_write_attr_real_scalar(group_id, name, value)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in) :: name
    real, intent(in) :: value

    integer(HID_T) :: aspace_id, attr_id
    integer(HSIZE_T) :: n8(1) = (/1/)
    integer ierr
    integer(4) :: rank = 1

    CALL h5screate_simple_f(rank, n8, aspace_id, ierr)
    CALL h5acreate_f(group_id, name, H5T_NATIVE_REAL, aspace_id, attr_id, ierr)
    CALL h5awrite_f(attr_id, H5T_NATIVE_REAL, (/value/), n8, ierr)
    CALL h5aclose_f(attr_id, ierr)
    CALL h5sclose_f(aspace_id, ierr)
  end subroutine hdf_write_attr_real_scalar

  subroutine hdf_read_attr_real(group_id, name, value, n)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in) :: name
    integer, intent(in) :: n
    real, intent(out) :: value(n)

    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: nn
    integer ierr

    nn(1) = n
    CALL h5aopen_name_f(group_id, name, attr_id, ierr)
    CALL h5aread_f(attr_id, H5T_NATIVE_REAL, value, nn, ierr)
    CALL h5aclose_f(attr_id, ierr)
  end subroutine hdf_read_attr_real

  subroutine hdf_read_attr_real_scalar(group_id, name, value)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in) :: name
    real, intent(out) :: value

    integer(HID_T) :: attr_id
    integer(HSIZE_T) :: n(1)
    integer ierr
    
    n(1) = 1
    CALL h5aopen_name_f(group_id, name, attr_id, ierr)
    CALL h5aread_f(attr_id, H5T_NATIVE_REAL, value, n, ierr)
    CALL h5aclose_f(attr_id, ierr)
  end subroutine hdf_read_attr_real_scalar
#endif

end module hdf_io
