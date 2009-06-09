#include <../config.h>
module normalize

  use grafic_types

contains

  subroutine mpnorm(local_nz,nz,ny,nx,total_local_size,buffer)

    ! Arguments
    integer, intent(in) :: local_nz,nz,ny,nx,total_local_size
#ifdef DOUB
    real(dp), dimension(total_local_size), intent(inout) :: buffer
#else
    real(sp), dimension(total_local_size), intent(inout) :: buffer
#endif

    ! Local variables
    integer :: myid, nproc, ierr
    real(dp) :: lmean, lstdev
    integer :: i1,i2,i3
    integer(i8b) :: index
    integer :: n2x
    real(dp) :: mean, stdev
    
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    n2x = 2*(nx/2+1)

    lmean=0
    do i3=1,local_nz
       do i2=1,ny
          do i1=1,nx ! Real space, ignore padding zone
             index=int((i3-1)*ny+i2-1,8)*n2x+i1
#ifdef DOUB
             lmean = lmean+buffer(index)
#else
             lmean = lmean+real(buffer(index),kind=dp)
#endif
          enddo
       enddo
    enddo

    call mpi_allreduce(lmean,mean,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    mean = mean/real(int(nx*ny,8)*nz,kind=dp)

    if (myid==0) then
       print*,'Mean value of the white noise box = ',mean
       print*,'Will be removed'
    endif
    
    lstdev=0
    do i3=1,local_nz
       do i2=1,ny
          do i1=1,nx ! Real space, skip padding zone
             index=int((i3-1)*ny+i2-1,8)*n2x+i1
#ifdef DOUB
             buffer(index) = buffer(index)-mean
             lstdev = lstdev+buffer(index)**2
#else
             buffer(index) = buffer(index)-real(mean,kind=sp)
             lstdev = lstdev + real(buffer(index),kind=8)**2
#endif
          enddo
       enddo
    enddo
    
    call mpi_allreduce(lstdev,stdev,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    stdev = stdev/real(int(nx*ny,8)*nz,kind=dp)
    stdev = sqrt(stdev)
    if (myid==0) then
       print*,'Standard deviation of the white noise box = ',stdev
       print*,'Will be factored out so that input stdev=1'
    endif
    
#ifdef DOUB
    buffer = buffer/stdev
#else
    buffer = buffer/real(stdev,kind=sp)
#endif
    
    
  end subroutine mpnorm

end module normalize
