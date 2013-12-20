#include <../config.h>
module transform

  use grafic_types
  IMPLICIT NONE

  logical, parameter, private :: use_workspace = .true.
#ifdef DOUB
  real(dp), allocatable, dimension(:) :: work
#else
  real(sp), allocatable, dimension(:) :: work
#endif
  

#ifdef ADD1US
#define  rfftw3d_f77_create_plan  rfftw3d_f77_create_plan_
#define  rfftwnd_f77_destroy_plan rfftwnd_f77_destroy_plan_
#define  rfftwnd_f77_one_real_to_complex rfftwnd_f77_one_real_to_complex_
#define  rfftw3d_f77_mpi_create_plan  rfftw3d_f77_mpi_create_plan_
#define  rfftwnd_f77_mpi_destroy_plan rfftwnd_f77_mpi_destroy_plan_
#define  rfftwnd_f77_mpi rfftwnd_f77_mpi_
#define  rfftwnd_f77_mpi_local_sizes rfftwnd_f77_mpi_local_sizes_
#endif

  interface fft_mpi
     module procedure fft_mpi_single, fft_mpi_double
  end interface

contains

  subroutine fft_forward(l,m,n,input,output)

    integer :: l,m,n
    real*8, dimension(l,m,n), intent(in) :: input
    complex*16, dimension(l/2+1,m,n), intent(out) :: output
    integer*8 :: plan

    call rfftw3d_create_plan(plan,l,m,n,FFTW_REAL_TO_COMPLEX, &
         & FFTW_ESTIMATE)
    call rfftwnd_one_real_to_complex(plan,input,output)
    call rfftwnd_destroy_plan(plan)

  end subroutine fft_forward

  ! This covers both forward and backward ffts, 
  ! depending on plan

  subroutine fft_mpi_double(plan,input)

    ! Arguments
    integer(i8b), intent(in) :: plan
    real(dp), dimension(:), intent(inout) :: input

    if (use_workspace) then
       call rfftwnd_mpi(plan, 1, input, work, 1, FFTW_NORMAL_ORDER) 
    else
       call rfftwnd_mpi(plan, 1, input, 0, 0, FFTW_NORMAL_ORDER) 
    endif

  end subroutine fft_mpi_double

  subroutine fft_mpi_single(plan,input)

    ! Arguments
    integer(i8b), intent(in) :: plan
    real(sp), dimension(:), intent(inout) :: input

    if (use_workspace) then
       call rfftwnd_mpi(plan, 1, input, work, 1, FFTW_NORMAL_ORDER) 
    else
       call rfftwnd_mpi(plan, 1, input, 0, 0, FFTW_NORMAL_ORDER) 
    endif

  end subroutine fft_mpi_single

  subroutine init_fftw(size)
    integer, intent(in) :: size
    if (.not.allocated(work) .and. use_workspace) allocate(work(size))
  end subroutine init_fftw

  subroutine cleanup_fftw
    if (allocated(work) .and. use_workspace) deallocate(work)
  end subroutine cleanup_fftw

end module transform
