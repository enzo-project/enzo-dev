#include <../config.h>
module grafic_types

  integer, parameter :: dp=selected_real_kind(12,200)
  integer, parameter :: sp=selected_real_kind(5,30)
  integer, parameter :: i8b=selected_int_kind(16)
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp))
  integer, parameter :: spc=kind((1.0_sp,1.0_sp))
  
  REAL(kind=dp), PARAMETER, public :: PI=3.141592653589793238462643383279502884197_dp

  type taille
     integer :: nx
     integer :: ny
     integer :: nz
     real(kind=sp) :: dx
     real(kind=sp) :: lx
     real(kind=sp) :: ly
     real(kind=sp) :: lz
  end type taille

  type cosmo
     real(kind=sp) :: astart
     real(kind=sp) :: omegab
     real(kind=sp) :: omegam
     real(kind=sp) :: omegav
     real(kind=sp) :: h0
     real(kind=sp) :: vfact
  end type cosmo

  ! FFTW stuff
  INTEGER,PARAMETER:: FFTW_FORWARD=-1, FFTW_BACKWARD=1

  INTEGER,PARAMETER:: FFTW_REAL_TO_COMPLEX=-1, FFTW_COMPLEX_TO_REAL=1

  INTEGER,PARAMETER:: FFTW_ESTIMATE=0, FFTW_MEASURE=1
  
  INTEGER,PARAMETER:: FFTW_OUT_OF_PLACE=0, FFTW_IN_PLACE=8

  INTEGER,PARAMETER:: FFTW_USE_WISDOM=16

  INTEGER,PARAMETER:: FFTW_THREADSAFE=128

  INTEGER,PARAMETER:: FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0
  INTEGER,PARAMETER:: FFMPI_COMM_WORLDTW_SCRAMBLED_INPUT=8192
  INTEGER,PARAMETER:: FFTW_SCRAMBLED_OUTPUT=16384

#ifdef ADD1US
#define  rfftw3d_f77_create_plan  rfftw3d_f77_create_plan_
#define  rfftwnd_f77_destroy_plan rfftwnd_f77_destroy_plan_
#define  rfftwnd_f77_one_real_to_complex rfftwnd_f77_one_real_to_complex_
#define  rfftw3d_f77_mpi_create_plan  rfftw3d_f77_mpi_create_plan_
#define  rfftwnd_f77_mpi_destroy_plan rfftwnd_f77_mpi_destroy_plan_
#define  rfftwnd_f77_mpi rfftwnd_f77_mpi_
#define  rfftwnd_f77_mpi_local_sizes rfftwnd_f77_mpi_local_sizes_ 
#endif

  INCLUDE 'mpif.h'


end module grafic_types
