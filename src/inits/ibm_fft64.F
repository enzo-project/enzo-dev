#include "../enzo/fortran.def"

#ifdef SP2

#ifdef CONFIG_BFLOAT_4

      subroutine ibm_1d(x, rank, n1, n2, n3, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1)

      REAL*4 :: factor
      REAL*4 :: scale
      REAL*8, allocatable :: aux1(:), aux2(:)

      integer*4 :: inc1x, inc2x, inc1y, inc2y
      integer*4 :: naux1, naux2, jdir, mult
      integer*4 :: m1, m2, m3, i0, i1

      m1 = n1
      m2 = n2
      m3 = n3
      jdir = -idir
      mult = 1
      inc1x = 1
      inc2x = 0
      inc1y = 1
      inc2y = 0
      i0 = 0
      i1 = 1

! SCFT

      if( n1 <= 8192 ) then
        naux1 = 20000
      else
        naux1 = 20000+(1.14_RKIND)*n1
      end if

      if( n1 <= 8192 ) then
        naux2 = 20000
      else
        naux2 = 20000+(1.14_RKIND)*n1
      end if 

      allocate ( aux1(naux1) )
      allocate ( aux2(naux2) )

      factor = 1._RKIND/REAL(n1,RKIND)

      if( jdir == 1 ) then
        scale = 1._RKIND
      else
        scale = factor
      end if

!     call scft(1,x,1,0,x,1,0,n1,1,jdir,scale,aux1,naux1,aux2,naux2)
!     call scft(0,x,1,0,x,1,0,n1,1,jdir,scale,aux1,naux1,aux2,naux2)

      call scft(i1,x,i1,i0,x,i1,i0,m1,i1,jdir,scale,
     &          aux1,naux1,aux2,naux2)
      call scft(i0,x,i1,i0,x,i1,i0,m1,i1,jdir,scale,
     &          aux1,naux1,aux2,naux2)

      deallocate ( aux1 )
      deallocate ( aux2 )

      return
      end




      subroutine ibm_2d(x, rank, n1, n2, n3, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1,n2)

      REAL*4 :: factor
      REAL*4 :: scale
      REAL*8, allocatable :: aux1(:), aux2(:)

      integer*4 :: inc1x, inc2x, inc1y, inc2y
      integer*4 :: naux1, naux2, jdir
      integer*4 :: m1, m2, m3, i0, i1

      m1 = n1
      m2 = n2
      m3 = n3
      jdir = -idir
      inc1x = 1
      inc1y = 1
      inc2x = n1
      inc2y = n1
      i0 = 0
      i1 = 1

! SCFT2

      if( max(n1,n2) <= 8192 ) then
        naux1 = 40000
      else
        naux1 = 40000+(1.14_RKIND)*(n1+n2)
      end if

      if( max(n1,n2) < 252 ) then
        naux2 = 20000
      else
        naux2 = 20000+(max(n1,n2)+256)*(min(64_IKIND,n1,n2)+1.14_RKIND)
      end if

      factor = 1._RKIND/(REAL(n1,RKIND)*REAL(n2,RKIND))

      if( jdir == 1 ) then
        scale = 1._RKIND
      else
        scale = factor
      end if

      allocate ( aux1(naux1) )
      allocate ( aux2(naux2) )

      call scft2(i1,x,inc1x,inc2x,x,inc1y,inc2y,m1,m2,
     &           jdir,scale,aux1,naux1,aux2,naux2)
      call scft2(i0,x,inc1x,inc2x,x,inc1y,inc2y,m1,m2,
     &           jdir,scale,aux1,naux1,aux2,naux2)

      deallocate ( aux1 )
      deallocate ( aux2 )

      return
      end




      subroutine ibm_3d( x, rank, n1, n2, n3, idir )

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1,n2,n3)

      REAL*4 :: factor
      REAL*4:: scale
      REAL*8, allocatable :: aux(:)

      integer*4 :: inc2x, inc3x, naux, jdir
      integer*4 :: m1, m2, m3, i0, i1

      m1 = n1
      m2 = n2
      m3 = n3
      jdir = -idir
      inc2x=n1
      inc3x=n1*n2
      i0 = 0
      i1 = 1

! SCFT3

      if( max(n2,n3) < 252 ) then
        if( n1 <= 8192 ) then
          naux = 60000
        else
          naux = 60000+(2.28_RKIND)*n1
        end if
      end if

      if(( n2 >= 252 ) .and. ( n3 < 252 )) then
        if( n1 <= 8192 ) then
          naux = 60000                +(n2+256)*(min(64,n1)+2.28_RKIND)
        else
          naux = 60000+(2.28_RKIND)*n1+(n2+256)*(min(64,n1)+2.28_RKIND)
        end if
      end if

      if(( n2 < 252 ) .and. ( n3 >= 252 )) then
        if( n1 <= 8192 ) then
          naux=60000                +(n3+256)*(min(64,n1*n2)+2.28_RKIND)
        else
          naux=60000+(2.28_RKIND)*n1+(n3+256)*(min(64,n1*n2)+2.28_RKIND)
        end if
      end if

      if(( n2 >= 252 ) .and. ( n3 >= 252 )) then

        naux = max(
     &        60000+(2.28_RKIND)*n1+(n2+256)*(min(64,n1)+2.28_RKIND),
     &        60000+(2.28_RKIND)*n1+(n3+256)*(min(64,n1*n2)+2.28_RKIND))

      end if

      write(*,*) naux

      factor = 1._RKIND/(REAL(n1,RKIND)*REAL(n2,RKIND)*REAL(n3,RKIND))

      if( jdir == 1 ) then
        scale = 1._RKIND
      else
        scale = factor
      end if

      allocate ( aux(naux) )
      call scft3(x, inc2x, inc3x, x, inc2x, inc3x, m1, m2, m3,
     &           jdir, scale, aux, naux)
      deallocate ( aux )

      end

#endif





#ifdef CONFIG_BFLOAT_8

      subroutine ibm_1d(x, rank, n1, n2, n3, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1)

      REAL*8 :: factor
      REAL*8 :: scale
      REAL*8, allocatable :: aux1(:), aux2(:)

      integer*4 :: inc1x, inc2x, inc1y, inc2y
      integer*4 :: naux1, naux2, jdir, mult
      integer*4 :: m1, m2, m3, i0, i1

      m1 = n1
      m2 = n2
      m3 = n3
      jdir = -idir
      mult = 1
      inc1x = 1
      inc2x = 0
      inc1y = 1
      inc2y = 0
      i0 = 0
      i1 = 1

! DCFT

      if( n1 <= 2048 ) then
        naux1 = 20000
      else
        naux1 = 20000+(2.28_RKIND)*n1
      end if

      if( n1 <= 2048 ) then
        naux2 = 20000
      else
        naux2 = 20000+(2.28_RKIND)*n1
      end if

      allocate ( aux1(naux1) )
      allocate ( aux2(naux2) )

      factor = 1._RKIND/REAL(n1,RKIND)

      if( jdir == 1 ) then
        scale = 1._RKIND
      else
        scale = factor
      end if

      call dcft(i1,x,i1,i0,x,i1,i0,m1,i1,jdir,scale,
     &          aux1,naux1,aux2,naux2)
      call dcft(i0,x,i1,i0,x,i1,i0,m1,i1,jdir,scale,
     &          aux1,naux1,aux2,naux2)

      deallocate ( aux1 )
      deallocate ( aux2 )

      return
      end




      subroutine ibm_2d(x, rank, n1, n2, n3, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1,n2)

      REAL*8 :: factor
      REAL*8 :: scale
      REAL*8, allocatable :: aux1(:), aux2(:)

      integer*4 :: inc1x, inc2x, inc1y, inc2y
      integer*4 :: naux1, naux2, jdir
      integer*4 :: m1, m2, m3, i0, i1

      m1 = n1
      m2 = n2
      m3 = n3
      jdir = -idir
      inc1x = 1
      inc1y = 1
      inc2x = n1
      inc2y = n1
      i0 = 0
      i1 = 1

! DCFT2

      if( max(n1,n2) <= 2048 ) then
        naux1 = 40000
      else
        naux1 = 40000+(2.28_RKIND)*(n1+n2)
      end if

      if( max(n1,n2) < 252 ) then
        naux2 = 20000
      else
        naux2 = 20000+(2*max(n1,n2)+256)*(min(64,n1,n2)+2.28_RKIND)
      end if

      factor = 1._RKIND/(REAL(n1,RKIND)*REAL(n2,RKIND))

      if( jdir == 1 ) then
        scale = 1._RKIND
      else
        scale = factor
      end if

      allocate ( aux1(naux1) )
      allocate ( aux2(naux2) )

      call dcft2(i1,x,inc1x,inc2x,x,inc1y,inc2y,m1,m2,
     &           jdir,scale,aux1,naux1,aux2,naux2)
      call dcft2(i0,x,inc1x,inc2x,x,inc1y,inc2y,m1,m2,
     &           jdir,scale,aux1,naux1,aux2,naux2)

      deallocate ( aux1 )
      deallocate ( aux2 )

      return
      end




      subroutine ibm_3d( x, rank, n1, n2, n3, idir )

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC :: rank, n1, n2, n3, idir
      CMPLX_PREC :: x(n1,n2,n3)

      REAL*8 :: factor
      REAL*8 :: scale
      REAL*8, allocatable :: aux(:)

      integer*4 :: inc2x, inc3x, naux, jdir
      integer*4 :: m1, m2, m3, i0, i1

      m1 = n1
      m2 = n2
      m3 = n3
      jdir = -idir
      inc2x=m1
      inc3x=m1*m2
      i0 = 0
      i1 = 1

! DCFT3

      if( max(n2,n3) < 252 ) then
        if( n1 <= 2048 ) then
          naux = 60000
        else
          naux = 60000+(4.56_RKIND)*n1
        end if
      end if

      if(( n2 >= 252 ) .and. ( n3 < 252 )) then
        if( n1 <= 2048 ) then
          naux = 60000+(2*n2+256)*(min(64,n1)+4.56_RKIND)
        else
          naux = 60000+(4.56_RKIND)*n1
     &          +(2*n2+256)*(min(64,n1)+4.56_RKIND)
        end if
      end if

      if(( n2 < 252 ) .and. ( n3 >= 252 )) then
        if( n1 <= 2048 ) then
          naux = 60000+(2*n3+256)*(min(64,n1*n2)+4.56_RKIND)
        else
          naux = 60000+(2.28_RKIND)*n1
     &          +(2*n3+256)*(min(64,n1*n2)+4.56_RKIND)
        end if
      end if

      if(( n2 >= 252 ) .and. ( n3 >= 252 )) then

        naux = max(
     &        60000+(4.56_RKIND)*n1+(2*n2+256)*(min(64,n1)+4.56_RKIND),
     &        60000+(4.56_RKIND)*n1
     &             +(2*n3+256)*(min(64,n1*n2)+4.56_RKIND) )

      end if

      write(*,*) naux

      factor = 1._RKIND/(REAL(n1,RKIND)*REAL(n2,RKIND)*REAL(n3,RKIND))

      if( jdir == 1 ) then
        scale = 1._RKIND
      else
        scale = factor
      end if

      allocate ( aux(naux) )
      call dcft3(x, inc2x, inc3x, x, inc2x, inc3x, m1, m2, m3,
     &           jdir, scale, aux, naux)
      deallocate ( aux )

      return
      end

#endif

#else

      subroutine ibm_1d(x, rank, n1, n2, n3, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC rank, n1, n2, n3, idir
      CMPLX_PREC x(n1)

      write(0,'("Dummy IBM 1D FFT - error")')
      call stop_all_cpus

      return
      end

      subroutine ibm_2d(x, rank, n1, n2, n3, idir)

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC rank, n1, n2, n3, idir
      CMPLX_PREC x(n1,n2)

      write(0,'("Dummy IBM 2D FFT - error")')
      call stop_all_cpus

      return
      end

      subroutine ibm_3d( x, rank, n1, n2, n3, idir )

      implicit none
#include "../enzo/fortran_types.def"

      INTG_PREC rank, n1, n2, n3, idir
      CMPLX_PREC x(n1,n2,n3)

      write(0,'("Dummy IBM 3D FFT - error")')
      call stop_all_cpus

      return
      end

#endif
