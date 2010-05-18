module enzo_io
#ifdef ENZO

  use grafic_types
  use hdf5
  use hdf_io

  integer, public, save :: enzo_handle = 77

  interface crop_buffer
     module procedure crop_buffer_single, crop_buffer_double
  end interface

#endif
contains
#ifdef ENZO

  subroutine enzo_filestat(filename, newname, filepart, n_filepart)
    implicit none
    character(128), intent(in)  :: filename
    character(*), intent(out) :: newname
    integer, intent(out) :: filepart, n_filepart
    character :: dimchar

    dimchar = filename(8:8)
    if (filename == "ic_deltab") then
       newname = "GridDensity"
       n_filepart = 1
       filepart = 0
    else if (filename(:7) == "ic_velb") then
#ifdef ONE_DIMENSION
       newname = "GridVelocities" // "_" // dimchar
       n_filepart = 1
       filepart = 0
#else
       newname = "GridVelocities"
       n_filepart = 3
       filepart = ichar(dimchar)-120
#endif
    else if (filename(:7) == "ic_velc") then
#ifdef ONE_DIMENSION
       newname = "ParticleVelocities" // "_" // dimchar
       n_filepart = 1
       filepart = 0
#else
       newname = "ParticleVelocities"
       n_filepart = 3
       filepart = ichar(dimchar)-120
#endif
    else ! white noise
       newname = trim(filename)
       n_filepart = 1
       filepart = 0
    endif
    newname = trim(newname)

  end subroutine enzo_filestat

  function get_ndims(filename)
    implicit none
    character(128), intent(in) :: filename
    integer :: get_ndims

    if (filename == "GridDensity") then
       get_ndims = 1
    else if (filename(:14) == "GridVelocities") then
#ifdef ONE_DIMENSION
       get_ndims = 1
#else
       get_ndims = 3
#endif
    else if (filename(:18) == "ParticleVelocities") then
#ifdef ONE_DIMENSION
       get_ndims = 1
#else
       get_ndims = 3
#endif
    else
       print*, "Unknown filename for enzo :: ", filename
       stop
    endif
  end function get_ndims

  function get_level_ext(filename)
    implicit none
    character(128), intent(in) :: filename
    character(2) :: get_level_ext
    integer :: dot_pos
    dot_pos = index(filename, ".")
    if (dot_pos.gt.0) then
       get_level_ext = filename(dot_pos:dot_pos+1)
    else
       get_level_ext = ""
    endif
  end function get_level_ext

  function get_dsetname(filename)
    character(128), intent(in) :: filename
    character(128) :: get_dsetname
    integer :: dot_pos
    dot_pos = index(filename, ".")
    if (dot_pos.gt.0) then
       get_dsetname = filename(1:dot_pos-1)
    else
       get_dsetname = filename
    endif
  end function get_dsetname

!-----------------------------------------------------------------------

  subroutine create_mask(mask, htin, htout, degramax, level, width, buffercells, &
       local_nz, local_z_start, write_paramfile)

    implicit none

    ! Arguments
    logical, allocatable, dimension(:,:,:), intent(inout) :: mask
    type(taille), intent(in) :: htin
    type(taille), intent(inout) :: htout
    real(sp), intent(in) :: width
    integer(i4b), intent(in) :: buffercells, level, degramax
    integer(i4b), intent(inout) :: local_nz, local_z_start
    logical, intent(in) :: write_paramfile

    ! Locals
    integer(i4b) :: maxlevel, this_width, i, start_idx, end_idx
    integer(i4b) :: nx, ny, nz, nx_trim, ny_trim, nz_trim, ntrim
    integer(i4b) :: local_x0, local_x1, local_y0, local_y1, local_z0, local_z1
    integer(i4b) :: border_up, border_bt, clip_bt, clip_up
    real :: left_edge, right_edge
#ifdef DOUB
    character(*), parameter :: format_str1 = '(A,"[",I0,"] = ",F15.12, F15.12, F15.12)'
    character(*), parameter :: format_str2 = '(A," = ",F15.12, F15.12, F15.12)'
#else
    character(*), parameter :: format_str1 = '(A,"[",I0,"] = ",F15.6, F15.6, F15.6)'
    character(*), parameter :: format_str2 = '(A," = ",F15.6, F15.6, F15.6)'
#endif

    maxlevel = nint(log(real(degramax)) / log(2.0))
    nx = htout%nx
    ny = htout%ny
    nz = htout%nz

    ! First create the mask with all .false. and then include all of
    ! the cells we want to write.
    if (allocated(mask)) then
       print*, 'mask: mask already allocated?! Break.'
       stop
    else
       allocate(mask(nx,ny,local_nz))
       mask = .false.
    endif

    ! Adjust the width of the coarsest subgrid to this level (in cell
    ! widths on this level)
    if (level.gt.0) then
       this_width = nint(htin%nx * width)  ! in cell widths of the finest grid
       this_width = ishft(this_width, level-maxlevel) ! cell widths of this level
       do i = 2, level  ! remove buffer cells from width
          this_width = this_width - 2*ishft(buffercells, level-i+1)
       enddo
    else
       this_width = ishft(htin%nx, level-maxlevel)
    endif

    start_idx = (nx - this_width)/2 + 1
    end_idx = (nx + this_width)/2

    ! Exit if this processor doesn't have any data inside the subgrid
!    if (start_idx.gt.local_z_start+local_nz .or. &
!         end_idx.lt.local_z_start) then
!       local_nz = 0
!       return
!    endif

    ! Crop to width of this subgrid
    local_x0 = start_idx
    local_x1 = end_idx
    local_y0 = start_idx
    local_y1 = end_idx
    local_z0 = max(min(start_idx - local_z_start, local_nz), 1)
    local_z1 = max(min(end_idx - local_z_start, local_nz), 1)
    mask(local_x0:local_x1, local_y0:local_y1, local_z0:local_z1) = .true.
    if (local_z_start+local_nz < start_idx .or. local_z_start >= end_idx) then
       local_nz = 0
       local_z_start = 0
    else
       local_nz = local_z1 - local_z0 + 1
       local_z_start = max(local_z_start - start_idx + 1, 0)
    end if

    ! Adjust header to reflect cropping only at the max level.  This
    ! is done in the main program for the coarse levels.
    htout%nx = this_width
    htout%ny = this_width
    htout%nz = this_width
    !htout%dx = htin%dx / 2**level

    ! Write part of the enzo parameter file if requested
    if (write_paramfile .and. level.gt.0) then
       left_edge = real(start_idx-1) / real(nx)
       right_edge = real(end_idx) / real(nx)
       if (level.eq.maxlevel) then
          write(enzo_handle, format_str2) &
               "RefineRegionLeftEdge ", left_edge, left_edge, left_edge
          write(enzo_handle, format_str2) &
               "RefineRegionRightEdge", right_edge, right_edge, right_edge
          write(enzo_handle, *) ''
          write(enzo_handle, '(A,I0)') &
               "CosmologySimulationNumberOfInitialGrids = ", maxlevel+1
       endif
       write(enzo_handle, '(A,"[",I0,"] = ",I6,I6,I6)') &
            "CosmologySimulationGridDimension", &
            level, this_width, this_width, this_width
       write(enzo_handle, format_str1) &
            "CosmologySimulationGridLeftEdge", &
            level, left_edge, left_edge, left_edge
       write(enzo_handle, format_str1) &
            "CosmologySimulationGridRightEdge", &
            level, right_edge, right_edge, right_edge
       write(enzo_handle, '(A,"[",I0,"]     = ",I1)') &
            "CosmologySimulationGridLevel", level, level
    endif

  end subroutine create_mask


  subroutine crop_buffer_single(dbuffer, workout, dims, mask)

    ! Arguments
    implicit none
    real(sp), dimension(:), allocatable, intent(in) :: dbuffer
    real(sp), dimension(:), allocatable, intent(out) :: workout
    integer, intent(inout) :: dims(3)
    logical, dimension(:,:,:), intent(in) :: mask
    
    ! Locals
    logical, dimension(:,:,:), allocatable :: workmask
    integer, dimension(2) :: xrange, yrange, zrange
    integer, dimension(3) :: old_dims
    integer :: i, j, k, index1, index2, n2x

    xrange(1) = 99999; xrange(2) = -99999
    yrange(1) = 99999; yrange(2) = -99999
    zrange(1) = 99999; zrange(2) = -99999

    if (dims(1).eq.dims(2)) then
       n2x = 2*(dims(1)/2+1)
    else
       n2x = dims(1)
    endif
    old_dims = (/n2x, dims(2), dims(3)/)

    ! Must reshape mask to include padding
    allocate(workmask(n2x, dims(2), dims(3)))
    workmask = .false.
    workmask(1:dims(1), 1:dims(2), 1:dims(3)) = mask

    do k = 1, dims(3)
       do j = 1, dims(2)
          do i = 1, dims(1)
             if (workmask(i,j,k)) then
                xrange(1) = min(xrange(1), i)
                xrange(2) = max(xrange(2), i)
                yrange(1) = min(yrange(1), j)
                yrange(2) = max(yrange(2), j)
                zrange(1) = min(zrange(1), k)
                zrange(2) = max(zrange(2), k)
             end if
          end do
       end do
    end do
    dims = (/xrange(2)-xrange(1), yrange(2)-yrange(1), zrange(2)-zrange(1)/)
    where (dims > 0)
       dims = dims+1
    end where

    allocate(workout(int(dims(1)*dims(2),8)*dims(3)))
    if (count(workmask) == size(workmask)) then
       workout = dbuffer
    else
       do k = 1, dims(3)
          do j = 1, dims(2)
             index1 = ((k-1)*dims(2) + j-1)*dims(1) + 1
             index2 = ((zrange(1)+k-2)*old_dims(2) + yrange(1)+j-2) * old_dims(1) + &
                  xrange(1)
             do i = 1, dims(1)
                workout(index1) = dbuffer(index2)
                index1 = index1 + 1
                index2 = index2 + 1
             enddo
          enddo
       enddo
    end if
    deallocate(workmask)

  end subroutine crop_buffer_single

  subroutine crop_buffer_double(dbuffer, workout, dims, mask)

    ! Arguments
    implicit none
    real(dp), dimension(:), allocatable, intent(in) :: dbuffer
    real(dp), dimension(:), allocatable, intent(out) :: workout
    integer, intent(inout) :: dims(3)
    logical, dimension(:,:,:), intent(in) :: mask
    
    ! Locals
    logical, dimension(:,:,:), allocatable :: workmask
    integer, dimension(2) :: xrange, yrange, zrange
    integer, dimension(3) :: old_dims
    integer :: i, j, k, index1, index2, n2x

    xrange(1) = 99999; xrange(2) = -99999
    yrange(1) = 99999; yrange(2) = -99999
    zrange(1) = 99999; zrange(2) = -99999

    if (dims(1).eq.dims(2)) then
       n2x = 2*(dims(1)/2+1)
    else
       n2x = dims(1)
    endif
    old_dims = (/n2x, dims(2), dims(3)/)

    ! Must reshape mask to include padding
    allocate(workmask(n2x, dims(2), dims(3)))
    workmask = .false.
    workmask(1:dims(1), 1:dims(2), 1:dims(3)) = mask

    do k = 1, dims(3)
       do j = 1, dims(2)
          do i = 1, dims(1)
             if (workmask(i,j,k)) then
                xrange(1) = min(xrange(1), i)
                xrange(2) = max(xrange(2), i)
                yrange(1) = min(yrange(1), j)
                yrange(2) = max(yrange(2), j)
                zrange(1) = min(zrange(1), k)
                zrange(2) = max(zrange(2), k)
             end if
          end do
       end do
    end do
    dims = (/xrange(2)-xrange(1), yrange(2)-yrange(1), zrange(2)-zrange(1)/)
    where (dims > 0)
       dims = dims+1
    end where

    allocate(workout(int(dims(1)*dims(2),8)*dims(3)))
    if (count(workmask) == size(workmask)) then
       workout = dbuffer
    else
       do k = 1, dims(3)
          do j = 1, dims(2)
             index1 = ((k-1)*dims(2) + j-1)*dims(1) + 1
             index2 = ((zrange(1)+k-2)*old_dims(2) + yrange(1)+j-2) * old_dims(1) + &
                  xrange(1)
             do i = 1, dims(1)
                workout(index1) = dbuffer(index2)
                index1 = index1 + 1
                index2 = index2 + 1
             enddo
          enddo
       enddo
    end if
    deallocate(workmask)

  end subroutine crop_buffer_double

#endif
end module enzo_io
