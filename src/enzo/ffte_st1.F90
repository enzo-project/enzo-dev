      subroutine ffte_st1(x, n1, idir)

      implicit none

      integer :: n1, idir
      complex :: x(n1)

      real :: factor
      real :: scale
      complex, allocatable :: work(:)

      integer :: nwork, jdir
      integer :: m1

      m1 = n1
      nwork = n1*2
      jdir = idir

      allocate( work(nwork) )

      call ffte_zfft1d(x, m1,    0, work)
      call ffte_zfft1d(x, m1, jdir, work)

      deallocate( work )

!     factor = 1.0/real(n1)

!     if( jdir == 1 ) then
!       scale = 1.0
!     else
!       scale = factor
!     end if

      return
      end
