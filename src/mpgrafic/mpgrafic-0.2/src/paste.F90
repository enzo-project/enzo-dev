#include <../config.h>
module paste

  use grafic_types
  use grafic_io
  !use derivative, only: init_kgrid,kx,ky,kz
  implicit none

contains

  subroutine grid_paste(local_z_start,local_nz,head_taille,head_cosmo,buffer,small_file_name)

    !Arguments
    integer, intent(in) :: local_z_start,local_nz
    type(cosmo), intent(in) :: head_cosmo
    type(taille), intent(in) :: head_taille
#ifdef DOUB
    real(dp), dimension(:), intent(inout) :: buffer
#else
    real(sp), dimension(:), intent(inout) :: buffer
#endif
    character(len=128), intent(in) :: small_file_name

    ! Local variables
    type(cosmo) :: small_head_cosmo
    type(taille) :: small_head_taille
    integer :: nx, nxs, ny, nys, nz, nzs, n2x, n2xs
    real(dp) :: dx, dxs
    integer :: myid, ierr
    logical :: debut=.false., fin=.false.
    integer :: local_nzs, local_z_starts, local_z_start_end
    integer :: i,j,k,kk
    integer(i8b) :: index_big, index_small
    
#ifdef DOUB
    real(dp), dimension(:), allocatable :: sbuffer
#else
    real(sp), dimension(:), allocatable :: sbuffer
#endif

    call init_grafic_io()

    ! Read small file header, and verify its compliance with big cube header
    call grafic_read_header(small_file_name,small_head_taille,small_head_cosmo)
    nx=head_taille%nx
    ny=head_taille%ny
    nz=head_taille%nz
    dx=head_taille%dx
    nxs=small_head_taille%nx
    nys=small_head_taille%ny
    nzs=small_head_taille%nz
    dxs=small_head_taille%dx

    n2x=2*(nx/2+1)
    n2xs=2*(nxs/2+1)

    call mpi_barrier(MPI_COMM_WORLD,ierr)
    
!!$    call mpi_comm_rank(MPI_COMM_WORLD,myid,ierr)
!!$    if (nxs > nx .or. nys > ny .or. nzs > nz) then
!!$       if (myid==0) print*,'Wrong box sizes, aborting'
!!$       call mpi_finalize(ierr)
!!$       stop
!!$    endif
!!$
!!$    if (nx*dx /= nxs*dxs) then
!!$       if (myid==0) print*,'Incompatible physical box sizes, aborting'
!!$       call mpi_finalize(ierr)
!!$       stop
!!$    endif
!!$
!!$    if (head_cosmo%astart/=small_head_cosmo%astart .or. &
!!$         & head_cosmo%omegam/=small_head_cosmo%omegam .or. &
!!$         & head_cosmo%omegav/=small_head_cosmo%omegav .or. &
!!$         & head_cosmo%h0/=small_head_cosmo%h0) then
!!$       if (myid==0) print*,'Incompatible cosmological parameters, aborting'
!!$       call mpi_finalize(ierr)
!!$       stop
!!$    endif

    ! verify that modes need to be pasted
    ! debut is .true. if positive kz modes need padding
    ! fin is .true. if negative kz modes need padding
    ! Beware that all values of the pair (debut,fin) can happen ...
    ! The pasting policy is to get rid of all Nyquist planes in the small
    ! box, and retain the corresponding (non-Nyquist) frequencies (positive
    ! and negative) of the large box.

    ! First element is a positive (small) kx frequency (Nyquist excluded)   
    if (local_z_start <= nzs/2-1) debut=.true.
    ! Last element is a negative (small) kx frequency
    if (local_z_start+local_nz-1 >= nz-nzs/2+1) fin=.true.

    if (debut) then
       !local_nzs = min(local_z_start+local_nz,nzs/2)-local_z_start
       local_nzs = min(local_z_start+local_nz,nzs)-local_z_start
       local_z_starts = local_z_start
       allocate(sbuffer(local_nzs*nys*n2xs))
       ! Read small k grid. This avoids complicated communication ...
       call grafic_read(sbuffer,local_nzs,local_z_starts,nzs,nys,nxs, &
            & small_file_name,padding_in=.true.,serial_in=.true.)
       ! Fill "debut": positive kx frequency modes
       do k=1,local_nzs
          ! Fill positive ky frequency modes, including ky=0, skipping nyq_y
          do j=1,nys/2
             do i=1,n2xs-2 ! Skip nyq_z
                index_big = int((k-1)*ny+j-1,8)*n2x+i
                index_small = ((k-1)*nys+j-1)*n2xs+i
                buffer(index_big)=sbuffer(index_small)
             enddo
          enddo
          ! Now fill negative ky frequency modes
          do j=nys/2+2,nys
             do i=1,n2xs-2
                index_big = int((k-1)*ny+j+ny-nys-1,8)*n2x+i
                index_small = ((k-1)*nys+j-1)*n2xs+i
                buffer(index_big)=sbuffer(index_small)
             enddo
          enddo
       enddo
       deallocate(sbuffer)
    endif

    if (fin) then
       local_z_start_end = max(local_z_start,nz-nzs/2+1)
       local_z_starts = local_z_start_end-nz+nzs ! Small cube index
       local_nzs = local_nz - (local_z_start_end - local_z_start)

       allocate(sbuffer(local_nzs*nys*n2xs))
       call grafic_read(sbuffer,local_nzs,local_z_starts,nzs,nys,nxs, &
            & small_file_name,padding_in=.true.,serial_in=.true.)
       ! Fill "fin": negative kx frequency modes
       do k=1,local_nzs
          kk = k+local_z_start_end-local_z_start
          ! Fill positive ky frequency modes, including ky=0, skipping nyq_y
          do j=1,nys/2
             do i=1,n2xs-2 ! Skip nyq_z
                index_big = int((kk-1)*ny+j-1,8)*n2x+i
                index_small = ((k-1)*nys+j-1)*n2xs+i
                buffer(index_big)=sbuffer(index_small)
             enddo
          enddo
          ! Fill negative ky frequency modes
          do j=nys/2+2,nys
             do i=1,n2xs-2
                index_big = int((kk-1)*ny+j+ny-nys-1,8)*n2x+i
                index_small = ((k-1)*nys+j-1)*n2xs+i
                buffer(index_big)=sbuffer(index_small)
             enddo
          enddo
       enddo
       deallocate(sbuffer)
    endif


  end subroutine grid_paste

end module paste
