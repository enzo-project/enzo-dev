!     ffte: a fast fourier transform package
!
!     (c) copyright software, 2000-2004, all rights reserved
!                by
!         daisuke takahashi
!         graduate school of systems and information engineering
!         university of tsukuba
!         1-1-1 tennodai, tsukuba, ibaraki 305-8573, japan
!         e-mail: daisuke@cs.tsukuba.ac.jp

!     Modified by Robert Harkness, SDSC, MArch 2006




      module ffte_param

      implicit none
#include "fortran_types.def"

      INTEGER, parameter :: fftkind = RKIND
      INTEGER, parameter :: fftintk = IKIND

! the maximum supported 2-d transform length is 65536.
      INTEGER(fftintk), parameter :: nda2=65536

! the maximum supported 3-d transform length is 4096.
      INTEGER(fftintk), parameter :: nda3=4096
      INTEGER(fftintk), parameter :: nda4=256

! the parameter nblk is a blocking parameter.
!     INTEGER(fftintk), parameter :: nblk=8  !(for pentiumiii and athlon)
      INTEGER(fftintk), parameter :: nblk=16 !(for pentium4, athlon xp, opteron,
                                    ! itanium and itanium2)

! the parameter np is a padding parameter to avoid cache
! conflicts in the fft routines.
!     INTEGER(fftintk), parameter :: np=2 !(for pentiumiii)
!     INTEGER(fftintk), parameter :: np=4 !(for athlon, athlon xp, opteron and itanium)
      INTEGER(fftintk), parameter :: np=8 !(for pentium4 and itanium2)

! size of l2 cache
      INTEGER(fftintk), parameter :: l2size=1048576

      end module ffte_param




!     ffte: a fast fourier transform package
!
!     (c) copyright software, 2000-2004, all rights reserved
!                by
!         daisuke takahashi
!         graduate school of systems and information engineering
!         university of tsukuba
!         1-1-1 tennodai, tsukuba, ibaraki 305-8573, japan
!         e-mail: daisuke@cs.tsukuba.ac.jp
!
!
!     1-d CMPLX_PREC fft routine
!
!     fortran77 source program
!
!     call zfft1d(a,n,iopt,b)
!
!     a(n) is CMPLX_PREC input/output vector (complex*16)
!     b(n) is work vector (complex*16)
!     n is the length of the transforms (integer*4)
!       -----------------------------------
!         n  =  (2**ip) * (3**iq) * (5**ir)
!       -----------------------------------
!     iopt  =  0 for initializing the coefficients (integer*4)
!           =  -1 for forward transform
!           =  +1 for inverse transform
!
!     written by daisuke takahashi

      subroutine ffte_zfft1d(a,n,iopt,b)

      use ffte_param

      implicit none

      INTG_PREC :: n, iopt
      CMPLX_PREC :: a(*),b(*)

      INTG_PREC :: i
      INTG_PREC :: nd
      INTG_PREC :: n1, n2, m1, m2
      INTG_PREC :: nw2, nw3, nw4
      R_PREC :: dn

      CMPLX_PREC :: c((nda2+np)*(nblk+1)+np)
      CMPLX_PREC :: w1(nda2/2+np),w2(nda2/2+np)
      CMPLX_PREC :: ww((nda2+np)*4+np)
      INTG_PREC :: ip(3),ip1(3),ip2(3)

      save w1,w2,ww

      call factor(n,ip)

      if (iopt .eq. 1) then
        do i = 1,n
          a(i) = conjg(a(i))
        end do
      end if

      if (n .le. min(l2size/16/3,nda2)) then
        if (iopt .eq. 0) then
          call settbl(w1,n)
          return
        end if
        call fft235(a,b,w1,n,ip)
      else
        do i = 1,3
          ip1(i) = (ip(i)+1)/2
          ip2(i) = ip(i)-ip1(i)
        end do
        n1 = (2**ip1(1))*(3**ip1(2))*(5**ip1(3))
        n2 = (2**ip2(1))*(3**ip2(2))*(5**ip2(3))
        if (2**ip1(1) .lt. nblk .or. 2**ip2(1) .lt. nblk) then
          m1 = min(n1,(2**(ip1(1)/2))*(3**(ip1(2)/2))*(5**(ip1(3)/2)))
          m2 = min(n2,(2**(ip2(1)/2))*(3**(ip2(2)/2))*(5**(ip2(3)/2)))
        else
          m1 = min(n1,max(nblk,2**(ip1(1)/2)))
          m2 = min(n2,max(nblk,2**(ip2(1)/2)))
        end if
        nw2 = m1*m2+np
        nw3 = nw2+m1*(n2/m2)+np
        nw4 = nw3+m2*(n1/m1)+np

        if (iopt .eq. 0) then
          call settbl(w1,n1)
          call settbl(w2,n2)
          call settbls(ww,ww(nw2+1),ww(nw3+1),ww(nw4+1),n1,n2,m1,m2)
          return
        end if

        nd = (n2+np)*nblk+np

!$omp parallel private(c)
        call zfft1d0(a,a,b,c,c(nd+1),w1,w2,ww,ww(nw2+1),ww(nw3+1),&
     &               ww(nw4+1),n1,n2,m1,m2,ip1,ip2)
!$omp end parallel

      end if

      if (iopt .eq. 1) then
        dn = 1.0_fftkind/REAL(n,fftkind)
        do i = 1,n
          a(i) = conjg(a(i))*dn
        end do
      end if

      return
      end subroutine ffte_zfft1d


      subroutine zfft1d0(a1,a2,b,c,d,w1,w2,ww1,ww2,ww3,ww4,n1,n2,m1,m2,&
     &                   ip1,ip2)

      use ffte_param

      implicit none

      INTG_PREC :: n1, n2, m1, m2
      CMPLX_PREC :: a1(n1,*),a2(n2,*),b(n1,*),c(n2+np,*),d(*)
      CMPLX_PREC :: w1(*),w2(*)
      CMPLX_PREC :: ww1(m1,*),ww2(m1,*),ww3(m2,*),ww4(n1/m1,*)
      CMPLX_PREC :: temp
      INTG_PREC :: ip1(*),ip2(*)

      INTG_PREC :: i, j
      INTG_PREC :: ii, jj
      INTG_PREC :: ij, ik, ir, is, ij0

!$omp do private(ij,ij0,ir,j,temp)
      do ii = 1,n1,nblk

        do jj = 1,n2,nblk
          do i = ii,min(ii+nblk-1,n1)
            do j = jj,min(jj+nblk-1,n2)
              c(j,i-ii+1) = a1(i,j)
            end do
          end do
        end do

        do i = ii,min(ii+nblk-1,n1)
          call fft235(c(1,i-ii+1),d,w2,n2,ip2)
        end do

        if (2**ip1(1) .lt. nblk .or. 2**ip2(1) .lt. nblk) then
          do is = 1,n2/m2
            do ik = 1,m2
              j = ik+(is-1)*m2
              do i = ii,min(ii+nblk-1,n1)
                ir = (i-1)/m1+1
                ij = mod(i-1,m1)+1
                b(i,j) = c(j,i-ii+1)*(ww1(ij,ik)*ww2(ij,is)&
     &                *ww3(ik,ir)*ww4(ir,is))
              end do
            end do
          end do
        else
          ir = (ii-1)/m1+1
          ij0 = mod(ii-1,m1)+1
          do is = 1,n2/m2
            do ik = 1,m2
              temp = ww3(ik,ir)*ww4(ir,is)
              j = ik+(is-1)*m2
              ij = ij0
              do i = ii,min(ii+nblk-1,n1)
                b(i,j) = c(j,i-ii+1)*(ww1(ij,ik)*ww2(ij,is)*temp)
                ij = ij+1
              end do
            end do
          end do
        end if
      end do

!$omp do
      do jj = 1,n2,nblk
        do j = jj,min(jj+nblk-1,n2)
          call fft235(b(1,j),c,w1,n1,ip1)
        end do
        do i = 1,n1
          do j = jj,min(jj+nblk-1,n2)
            a2(j,i) = b(i,j)
          end do
        end do
      end do

      return
      end subroutine zfft1d0


      subroutine settbls(w1,w2,w3,w4,n1,n2,m1,m2)

      use ffte_param

      implicit none

      INTG_PREC :: n1, n2, m1, m2
      R_PREC :: w1(2,m1,*),w2(2,m1,*),w3(2,m2,*),w4(2,n1/m1,*)

      INTG_PREC :: j, k
      INTG_PREC :: ir, is
      R_PREC :: pi2, px

      pi2 = 8.0_fftkind*atan(1.0_fftkind)
      px = -pi2/(REAL(n1,fftkind)*REAL(n2,fftkind))

!$omp parallel
!$omp do
      do k = 1,m2
        do j = 1,m1
          w1(1,j,k) = cos(px*REAL(j-1,fftkind)*REAL(k-1,fftkind))
          w1(2,j,k) = sin(px*REAL(j-1,fftkind)*REAL(k-1,fftkind))
        end do
        do ir = 1,n1/m1
          w3(1,k,ir) = cos(px*REAL(k-1,fftkind)*REAL(ir-1,fftkind)*REAL(m1,fftkind))
          w3(2,k,ir) = sin(px*REAL(k-1,fftkind)*REAL(ir-1,fftkind)*REAL(m1,fftkind))
        end do
      end do
      do is = 1,n2/m2
        do j = 1,m1
          w2(1,j,is) = cos(px*REAL(j-1,fftkind)*REAL(is-1,fftkind)*REAL(m2,fftkind))
          w2(2,j,is) = sin(px*REAL(j-1,fftkind)*REAL(is-1,fftkind)*REAL(m2,fftkind))
        end do
        do ir = 1,n1/m1
          w4(1,ir,is) = cos(px*REAL(ir-1,fftkind)*REAL(m1,fftkind)*REAL(is-1,fftkind)*REAL(m2,fftkind))
          w4(2,ir,is) = sin(px*REAL(ir-1,fftkind)*REAL(m1,fftkind)*REAL(is-1,fftkind)*REAL(m2,fftkind))
        end do
      end do
!$omp end parallel

      return
      end subroutine settbls



!     ffte: a fast fourier transform package
!
!     (c) copyright software, 2000-2004, all rights reserved
!                by
!         daisuke takahashi
!         graduate school of systems and information engineering
!         university of tsukuba
!         1-1-1 tennodai, tsukuba, ibaraki 305-8573, japan
!         e-mail: daisuke@cs.tsukuba.ac.jp
!
!
!     2-d CMPLX_PREC fft routine
!
!     fortran77 source program
!
!     call zfft2d(a,nx,ny,iopt)
!
!     a(nx,ny) is CMPLX_PREC input/output vector (complex*16)
!     nx is the length of the transforms in the x-direction (integer*4)
!     ny is the length of the transforms in the y-direction (integer*4)
!       ------------------------------------
!         nx  =  (2**ip) * (3**iq) * (5**ir)
!         ny  =  (2**jp) * (3**jq) * (5**jr)
!       ------------------------------------
!     iopt  =  0 for initializing the coefficients (integer*4)
!           =  -1 for forward transform
!           =  +1 for inverse transform
!
!     written by daisuke takahashi

      subroutine ffte_zfft2d(a,nx,ny,iopt)

      use ffte_param

      implicit none

      INTG_PREC :: nx, ny, iopt
      CMPLX_PREC :: a(*)
      CMPLX_PREC :: b((nda2+np)*(nblk+1)+np)
      CMPLX_PREC :: wx(nda2/2+np),wy(nda2/2+np)
      INTG_PREC :: lnx(3),lny(3)

      INTG_PREC :: i, nc
      R_PREC :: dn

      save wx,wy

      call factor(nx,lnx)
      call factor(ny,lny)

      if (iopt .eq. 0) then
        call settbl(wx,nx)
        call settbl(wy,ny)
        return
      end if

      if (iopt .eq. 1) then
        do i = 1,nx*ny
          a(i) = conjg(a(i))
        end do
      end if

      nc = (ny+np)*nblk+np
!$omp parallel private(b)
      call zfft2d0(a,b,b(nc+1),wx,wy,nx,ny,lnx,lny)
!$omp end parallel

      if (iopt .eq. 1) then
        dn = 1.0_fftkind/(REAL(nx,fftkind)*REAL(ny,fftkind))
        do i = 1,nx*ny
          a(i) = conjg(a(i))*dn
        end do
      end if

      return
      end subroutine ffte_zfft2d


      subroutine zfft2d0(a,b,c,wx,wy,nx,ny,lnx,lny)

      use ffte_param

      implicit none

      INTG_PREC :: nx, ny
      CMPLX_PREC :: a(nx,*),b(ny+np,*),c(*)
      CMPLX_PREC :: wx(*),wy(*)
      INTG_PREC :: lnx(*),lny(*)

      INTG_PREC :: i, j
      INTG_PREC :: ii, jj

!$omp do
      do ii = 1,nx,nblk
        do jj = 1,ny,nblk
          do i = ii,min(ii+nblk-1,nx)
            do j = jj,min(jj+nblk-1,ny)
              b(j,i-ii+1) = a(i,j)
            end do
          end do
        end do
        do i = ii,min(ii+nblk-1,nx)
          call fft235(b(1,i-ii+1),c,wy,ny,lny)
        end do
        do j = 1,ny
          do i = ii,min(ii+nblk-1,nx)
            a(i,j) = b(j,i-ii+1)
          end do
        end do
      end do
!$omp do
      do j = 1,ny
        call fft235(a(1,j),b,wx,nx,lnx)
      end do

      return
      end subroutine zfft2d0



!     ffte: a fast fourier transform package
!
!     (c) copyright software, 2000-2004, all rights reserved
!                by
!         daisuke takahashi
!         graduate school of systems and information engineering
!         university of tsukuba
!         1-1-1 tennodai, tsukuba, ibaraki 305-8573, japan
!         e-mail: daisuke@cs.tsukuba.ac.jp
!
!
!     3-d CMPLX_PREC fft routine
!
!     fortran77 source program
!
!     call zfft3d(a,nx,ny,nz,iopt)
!
!     a(nx,ny,nz) is CMPLX_PREC input/output vector (complex*16)
!     nx is the length of the transforms in the x-direction (integer*4)
!     ny is the length of the transforms in the y-direction (integer*4)
!     nz is the length of the transforms in the z-direction (integer*4)
!       ------------------------------------
!         nx  =  (2**ip) * (3**iq) * (5**ir)
!         ny  =  (2**jp) * (3**jq) * (5**jr)
!         nz  =  (2**kp) * (3**kq) * (5**kr)
!       ------------------------------------
!     iopt  =  0 for initializing the coefficients (integer*4)
!           =  -1 for forward transform
!           =  +1 for inverse transform
!
!     written by daisuke takahashi

      subroutine ffte_zfft3d(a,nx,ny,nz,iopt)

      use ffte_param

      implicit none

      INTG_PREC :: nx, ny, nz, iopt
      CMPLX_PREC :: a(*)
      CMPLX_PREC :: b((nda3+np)*(nblk+1)+np)
      CMPLX_PREC :: wx(nda3/2+np),wy(nda3/2+np),wz(nda3/2+np)
      INTG_PREC :: lnx(3),lny(3),lnz(3)

      INTG_PREC :: i, nc
      R_PREC :: dn

      save wx,wy,wz

      call factor(nx,lnx)
      call factor(ny,lny)
      call factor(nz,lnz)

      if (iopt .eq. 0) then
        call settbl(wx,nx)
        call settbl(wy,ny)
        call settbl(wz,nz)
        return
      end if

      if (iopt .eq. 1) then
        do i = 1,nx*ny*nz
          a(i) = conjg(a(i))
        end do
      end if

      nc = (max(ny,nz)+np)*nblk+np
!$omp parallel private(b)
      call zfft3d0(a,b,b,b(nc+1),wx,wy,wz,nx,ny,nz,lnx,lny,lnz)
!$omp end parallel

      if (iopt .eq. 1) then
        dn = 1.0_fftkind/(REAL(nx,fftkind)*REAL(ny,fftkind)*REAL(nz,fftkind))
        do i = 1,nx*ny*nz
          a(i) = conjg(a(i))*dn
        end do
      end if

      return
      end subroutine ffte_zfft3d


      subroutine zfft3d0(a,by,bz,c,wx,wy,wz,nx,ny,nz,lnx,lny,lnz)

      use ffte_param

      implicit none

      INTG_PREC :: nx, ny, nz
      CMPLX_PREC :: a(nx,ny,*),by(ny+np,*),bz(nz+np,*),c(*)
      CMPLX_PREC :: wx(*),wy(*),wz(*)
      INTG_PREC :: lnx(*),lny(*),lnz(*)

      INTG_PREC :: i, j, k
      INTG_PREC :: ii, jj, kk

!$omp do
      do j = 1,ny
        do ii = 1,nx,nblk
          do kk = 1,nz,nblk
            do i = ii,min(ii+nblk-1,nx)
              do k = kk,min(kk+nblk-1,nz)
                bz(k,i-ii+1) = a(i,j,k)
              end do
            end do
          end do
          do i = ii,min(ii+nblk-1,nx)
            call fft235(bz(1,i-ii+1),c,wz,nz,lnz)
          end do
          do k = 1,nz
            do i = ii,min(ii+nblk-1,nx)
              a(i,j,k) = bz(k,i-ii+1)
            end do
          end do
        end do
      end do
!$omp do
      do k = 1,nz
        do ii = 1,nx,nblk
          do jj = 1,ny,nblk
            do i = ii,min(ii+nblk-1,nx)
              do j = jj,min(jj+nblk-1,ny)
                by(j,i-ii+1) = a(i,j,k)
              end do
            end do
          end do
          do i = ii,min(ii+nblk-1,nx)
            call fft235(by(1,i-ii+1),c,wy,ny,lny)
          end do
          do j = 1,ny
            do i = ii,min(ii+nblk-1,nx)
              a(i,j,k) = by(j,i-ii+1)
            end do
          end do
        end do
        do j = 1,ny
          call fft235(a(1,j,k),c,wx,nx,lnx)
        end do
      end do

      return
      end subroutine zfft3d0



!     ffte: a fast fourier transform package
!
!     (c) copyright software, 2000-2004, all rights reserved
!                by
!         daisuke takahashi
!         graduate school of systems and information engineering
!         university of tsukuba
!         1-1-1 tennodai, tsukuba, ibaraki 305-8573, japan
!         e-mail: daisuke@cs.tsukuba.ac.jp
!
!
!     radix-2, 3, 4, 5 and 8 fft routine
!
!     fortran77 source program
!
!     written by daisuke takahashi

      subroutine fft235(a,b,w,n,ip)

      use ffte_param

      implicit none

      INTG_PREC :: n
      INTG_PREC :: ip(*)
      CMPLX_PREC :: a(*),b(*),w(*)

      INTG_PREC :: j, k, l, m
      INTG_PREC :: key
      INTG_PREC :: kp4, kp8

      if (ip(1) .ne. 1) then
        kp4 = 2-mod(ip(1)+2,3_fftintk)
        kp8 = (ip(1)-kp4)/3
      else
        kp4 = 0
        kp8 = 0
      end if

      key = 1
      j = 1
      l = n
      m = 1
      do 10 k = 1,kp8
        l = l/8
        if (l .ge. 2) then
          if (key .ge. 0) then
            call fft8(a,b,w(j),m,l)
          else
            call fft8(b,a,w(j),m,l)
          end if
          key = -key
        else
          if (key .ge. 0) then
            call fft8(a,a,w(j),m,l)
          else
            call fft8(b,a,w(j),m,l)
          end if
        end if
        m = m*8
        j = j+l
   10 continue

      do 20 k = 1,ip(3)
        l = l/5
        if (l .ge. 2) then
          if (key .ge. 0) then
            call fft5(a,b,w(j),m,l)
          else
            call fft5(b,a,w(j),m,l)
          end if
          key = -key
        else
          if (key .ge. 0) then
            call fft5(a,a,w(j),m,l)
          else
            call fft5(b,a,w(j),m,l)
          end if
        end if
        m = m*5
        j = j+l
   20 continue

      do 30 k = 1,kp4
        l = l/4
        if (l .ge. 2) then
          if (key .ge. 0) then
            call fft4(a,b,w(j),m,l)
          else
            call fft4(b,a,w(j),m,l)
          end if
          key = -key
        else
          if (key .ge. 0) then
            call fft4(a,a,w(j),m,l)
          else
            call fft4(b,a,w(j),m,l)
          end if
        end if
        m = m*4
        j = j+l
   30 continue

      do 40 k = 1,ip(2)
        l = l/3
        if (l .ge. 2) then
          if (key .ge. 0) then
            call fft3(a,b,w(j),m,l)
          else
            call fft3(b,a,w(j),m,l)
          end if
          key = -key
        else
          if (key .ge. 0) then
            call fft3(a,a,w(j),m,l)
          else
            call fft3(b,a,w(j),m,l)
          end if
        end if
        m = m*3
        j = j+l
   40 continue

      if (ip(1) .eq. 1) then
        if (key .ge. 0) then
          call fft2(a,a,m)
        else
          call fft2(b,a,m)
        end if
      end if

      return
      end subroutine fft235


      subroutine fft3(a,b,w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      CMPLX_PREC :: a(*),b(*),w(*)

      if (m .eq. 1) then
        call fft3a(a,b,w,l)
      else
        call fft3b(a,b,w,m,l)
      end if

      return
      end subroutine fft3


      subroutine fft4(a,b,w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      CMPLX_PREC :: a(*),b(*),w(*)

      if (m .eq. 1) then
        call fft4a(a,b,w,l)
      else
        call fft4b(a,b,w,m,l)
      end if

      return
      end subroutine fft4


      subroutine fft5(a,b,w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      CMPLX_PREC :: a(*),b(*),w(*)

      if (m .eq. 1) then
        call fft5a(a,b,w,l)
      else
        call fft5b(a,b,w,m,l)
      end if

      return
      end subroutine fft5


      subroutine fft8(a,b,w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      CMPLX_PREC :: a(*),b(*),w(*)

      if (m .eq. 1) then
        call fft8a(a,b,w,l)
      else
        call fft8b(a,b,w,m,l)
      end if

      return
      end subroutine fft8


      subroutine settbl(w,n)

      use ffte_param

      implicit none

      INTG_PREC :: n
      CMPLX_PREC :: w(*)

      INTG_PREC :: i, j, k, l
      INTG_PREC :: ip(3)
      INTG_PREC :: kp4, kp8

      call factor(n,ip)

      if (ip(1) .ne. 1) then
        kp4 = 2-mod(ip(1)+2,3_fftintk)
        kp8 = (ip(1)-kp4)/3
      else
        kp4 = 0
        kp8 = 0
      end if

      j = 1
      l = n

      do 10 k = 1,kp8
        l = l/8
        call settbl0(w(j),8_IKIND,l)
        j = j+l
   10 continue

      do 20 k = 1,ip(3)
        l = l/5
        call settbl0(w(j),5_IKIND,l)
        j = j+l
   20 continue

      do 30 k = 1,kp4
        l = l/4
        call settbl0(w(j),4_IKIND,l)
        j = j+l
   30 continue

      do 40 k = 1,ip(2)
        l = l/3
        call settbl0(w(j),3_IKIND,l)
        j = j+l
   40 continue

      return
      end subroutine settbl


      subroutine settbl0(w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      R_PREC :: w(2,*)
      INTG_PREC :: i
      R_PREC :: pi2, px

      pi2 = 8.0_fftkind*atan(1.0_fftkind)
      px = -pi2/(REAL(m,fftkind)*REAL(l,fftkind))

      do 10 i = 1,l
        w(1,i) = cos(px*REAL(i-1,fftkind))
        w(2,i) = sin(px*REAL(i-1,fftkind))
   10 continue

      return
      end subroutine settbl0


      subroutine settbl2(w,n1,n2)

      use ffte_param

      implicit none

      INTG_PREC :: n1, n2
      R_PREC :: w(2,n1,*)
      INTG_PREC :: j, k
      R_PREC :: pi2, px

      pi2 = 8.0_fftkind*atan(1.0_fftkind)
      px = -pi2/(REAL(n1,fftkind)*REAL(n2,fftkind))

      do 20 k = 1,n2
        do 10 j = 1,n1
          w(1,j,k) = cos(px*REAL(j-1,fftkind)*REAL(k-1,fftkind))
          w(2,j,k) = sin(px*REAL(j-1,fftkind)*REAL(k-1,fftkind))
   10   continue
   20 continue

      return
      end subroutine settbl2


      subroutine factor(n,ip)

      use ffte_param

      implicit none

      INTG_PREC :: n, n2
      INTG_PREC :: ip(*)

      ip(1) = 0
      ip(2) = 0
      ip(3) = 0
      n2 = n

      if (mod(n,2_fftintk) .ne. 0 .and. &
          mod(n,3_fftintk) .ne. 0 .and. &
     &    mod(n,5_fftintk) .ne. 0) return

   10 continue

      if (n2 .le. 1) return
      if (mod(n2,2_fftintk) .eq. 0) then
        ip(1) = ip(1)+1
        n2 = n2/2
        go to 10
      else if (mod(n2,3_fftintk) .eq. 0) then
        ip(2) = ip(2)+1
        n2 = n2/3
        go to 10
      else if (mod(n2,5_fftintk) .eq. 0) then
        ip(3) = ip(3)+1
        n2 = n2/5
        go to 10
      end if

      return
      end subroutine factor



!     ffte: a fast fourier transform package
!
!     (c) copyright software, 2000-2004, all rights reserved
!                by
!         daisuke takahashi
!         graduate school of systems and information engineering
!         university of tsukuba
!         1-1-1 tennodai, tsukuba, ibaraki 305-8573, japan
!         e-mail: daisuke@cs.tsukuba.ac.jp
!
!
!     radix-2, 3, 4, 5 and 8 fft kernel routine
!
!     fortran77 source program
!
!     written by daisuke takahashi

      subroutine fft2(a,b,m)

      use ffte_param

      implicit none

      INTG_PREC :: m
      R_PREC :: a(2,m,*),b(2,m,*)

      INTG_PREC :: i
      R_PREC :: x0, y0, x1, y1

      do i = 1,m
        x0 = a(1,i,1)
        y0 = a(2,i,1)
        x1 = a(1,i,2)
        y1 = a(2,i,2)
        b(1,i,1) = x0+x1
        b(2,i,1) = y0+y1
        b(1,i,2) = x0-x1
        b(2,i,2) = y0-y1
      end do

      return
      end subroutine fft2


      subroutine fft3a(a,b,w,l)

      use ffte_param

      implicit none

      INTG_PREC :: l
      R_PREC :: a(2,l,*),b(2,3,*),w(2,*)

      INTG_PREC :: i, j
      R_PREC :: wr1, wi1, wr2, wi2
      R_PREC :: x0, y0, x1, y1, x2, y2

      R_PREC :: c31, c32
      data c31/0.86602540378443865_fftkind/
      data c32/0.5_fftkind/

      do j = 1,l
        wr1 = w(1,j)
        wi1 = w(2,j)
        wr2 = wr1*wr1-wi1*wi1
        wi2 = wr1*wi1+wr1*wi1
        x0 = a(1,j,2)+a(1,j,3)
        y0 = a(2,j,2)+a(2,j,3)
        x1 = a(1,j,1)-c32*x0
        y1 = a(2,j,1)-c32*y0
        x2 = c31*(a(2,j,2)-a(2,j,3))
        y2 = c31*(a(1,j,3)-a(1,j,2))
        b(1,1,j) = a(1,j,1)+x0
        b(2,1,j) = a(2,j,1)+y0
        b(1,2,j) = wr1*(x1+x2)-wi1*(y1+y2)
        b(2,2,j) = wr1*(y1+y2)+wi1*(x1+x2)
        b(1,3,j) = wr2*(x1-x2)-wi2*(y1-y2)
        b(2,3,j) = wr2*(y1-y2)+wi2*(x1-x2)
      end do

      return
      end subroutine fft3a


      subroutine fft3b(a,b,w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      R_PREC :: a(2,m,l,*),b(2,m,3,*),w(2,*)

      INTG_PREC :: i, j
      R_PREC :: wr1, wi1, wr2, wi2
      R_PREC :: x0, y0, x1, y1, x2, y2

      R_PREC :: c31, c32
      data c31/0.86602540378443865_fftkind/
      data c32/0.5_fftkind/

      do i = 1,m
        x0 = a(1,i,1,2)+a(1,i,1,3)
        y0 = a(2,i,1,2)+a(2,i,1,3)
        x1 = a(1,i,1,1)-c32*x0
        y1 = a(2,i,1,1)-c32*y0
        x2 = c31*(a(2,i,1,2)-a(2,i,1,3))
        y2 = c31*(a(1,i,1,3)-a(1,i,1,2))
        b(1,i,1,1) = a(1,i,1,1)+x0
        b(2,i,1,1) = a(2,i,1,1)+y0
        b(1,i,2,1) = x1+x2
        b(2,i,2,1) = y1+y2
        b(1,i,3,1) = x1-x2
        b(2,i,3,1) = y1-y2
      end do

      do j = 2,l
        wr1 = w(1,j)
        wi1 = w(2,j)
        wr2 = wr1*wr1-wi1*wi1
        wi2 = wr1*wi1+wr1*wi1
        do i = 1,m
          x0 = a(1,i,j,2)+a(1,i,j,3)
          y0 = a(2,i,j,2)+a(2,i,j,3)
          x1 = a(1,i,j,1)-c32*x0
          y1 = a(2,i,j,1)-c32*y0
          x2 = c31*(a(2,i,j,2)-a(2,i,j,3))
          y2 = c31*(a(1,i,j,3)-a(1,i,j,2))
          b(1,i,1,j) = a(1,i,j,1)+x0
          b(2,i,1,j) = a(2,i,j,1)+y0
          b(1,i,2,j) = wr1*(x1+x2)-wi1*(y1+y2)
          b(2,i,2,j) = wr1*(y1+y2)+wi1*(x1+x2)
          b(1,i,3,j) = wr2*(x1-x2)-wi2*(y1-y2)
          b(2,i,3,j) = wr2*(y1-y2)+wi2*(x1-x2)
        end do
      end do

      return
      end subroutine fft3b


      subroutine fft4a(a,b,w,l)

      use ffte_param

      implicit none

      INTG_PREC :: l
      R_PREC :: a(2,l,*),b(2,4,*),w(2,*)

      INTG_PREC :: j
      R_PREC :: wr1, wi1, wr2, wi2, wr3, wi3
      R_PREC :: x0, y0, x1, y1, x2, y2, x3, y3

      do j = 1,l
        wr1 = w(1,j)
        wi1 = w(2,j)
        wr2 = wr1*wr1-wi1*wi1
        wi2 = wr1*wi1+wr1*wi1
        wr3 = wr1*wr2-wi1*wi2
        wi3 = wr1*wi2+wi1*wr2
        x0 = a(1,j,1)+a(1,j,3)
        y0 = a(2,j,1)+a(2,j,3)
        x1 = a(1,j,1)-a(1,j,3)
        y1 = a(2,j,1)-a(2,j,3)
        x2 = a(1,j,2)+a(1,j,4)
        y2 = a(2,j,2)+a(2,j,4)
        x3 = a(2,j,2)-a(2,j,4)
        y3 = a(1,j,4)-a(1,j,2)
        b(1,1,j) = x0+x2
        b(2,1,j) = y0+y2
        b(1,3,j) = wr2*(x0-x2)-wi2*(y0-y2)
        b(2,3,j) = wr2*(y0-y2)+wi2*(x0-x2)
        b(1,2,j) = wr1*(x1+x3)-wi1*(y1+y3)
        b(2,2,j) = wr1*(y1+y3)+wi1*(x1+x3)
        b(1,4,j) = wr3*(x1-x3)-wi3*(y1-y3)
        b(2,4,j) = wr3*(y1-y3)+wi3*(x1-x3)
      end do

      return
      end subroutine fft4a


!pgi$r opt=1
      subroutine fft4b(a,b,w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      R_PREC :: a(2,m,l,*),b(2,m,4,*),w(2,*)

      INTG_PREC :: i, j
      R_PREC :: wr1, wi1, wr2, wi2, wr3, wi3
      R_PREC :: x0, y0, x1, y1, x2, y2, x3, y3

      do i = 1,m
        x0 = a(1,i,1,1)+a(1,i,1,3)
        y0 = a(2,i,1,1)+a(2,i,1,3)
        x1 = a(1,i,1,1)-a(1,i,1,3)
        y1 = a(2,i,1,1)-a(2,i,1,3)
        x2 = a(1,i,1,2)+a(1,i,1,4)
        y2 = a(2,i,1,2)+a(2,i,1,4)
        x3 = a(2,i,1,2)-a(2,i,1,4)
        y3 = a(1,i,1,4)-a(1,i,1,2)
        b(1,i,1,1) = x0+x2
        b(2,i,1,1) = y0+y2
        b(1,i,3,1) = x0-x2
        b(2,i,3,1) = y0-y2
        b(1,i,2,1) = x1+x3
        b(2,i,2,1) = y1+y3
        b(1,i,4,1) = x1-x3
        b(2,i,4,1) = y1-y3
      end do

      do j = 2,l
        wr1 = w(1,j)
        wi1 = w(2,j)
        wr2 = wr1*wr1-wi1*wi1
        wi2 = wr1*wi1+wr1*wi1
        wr3 = wr1*wr2-wi1*wi2
        wi3 = wr1*wi2+wi1*wr2
        do i = 1,m
          x0 = a(1,i,j,1)+a(1,i,j,3)
          y0 = a(2,i,j,1)+a(2,i,j,3)
          x1 = a(1,i,j,1)-a(1,i,j,3)
          y1 = a(2,i,j,1)-a(2,i,j,3)
          x2 = a(1,i,j,2)+a(1,i,j,4)
          y2 = a(2,i,j,2)+a(2,i,j,4)
          x3 = a(2,i,j,2)-a(2,i,j,4)
          y3 = a(1,i,j,4)-a(1,i,j,2)
          b(1,i,1,j) = x0+x2
          b(2,i,1,j) = y0+y2
          b(1,i,3,j) = wr2*(x0-x2)-wi2*(y0-y2)
          b(2,i,3,j) = wr2*(y0-y2)+wi2*(x0-x2)
          b(1,i,2,j) = wr1*(x1+x3)-wi1*(y1+y3)
          b(2,i,2,j) = wr1*(y1+y3)+wi1*(x1+x3)
          b(1,i,4,j) = wr3*(x1-x3)-wi3*(y1-y3)
          b(2,i,4,j) = wr3*(y1-y3)+wi3*(x1-x3)
        end do
      end do

      return
      end subroutine fft4b


      subroutine fft5a(a,b,w,l)

      use ffte_param

      implicit none

      INTG_PREC :: l
      R_PREC :: a(2,l,*),b(2,5,*),w(2,*)

      INTG_PREC :: j
      R_PREC :: wr1, wi1, wr2, wi2, wr3, wi3, wr4, wi4
      R_PREC :: x0, y0, x1, y1, x2, y2, x3, y3
      R_PREC :: x4, y4, x5, y5, x6, y6, x7, y7
      R_PREC :: x8, y8, x9, y9, x10, y10

      R_PREC :: c51, c52, c53, c54
      data c51/0.95105651629515357_fftkind/
      data c52/0.61803398874989485_fftkind/
      data c53/0.55901699437494742_fftkind/
      data c54/0.25_fftkind/

      do j = 1,l
        wr1 = w(1,j)
        wi1 = w(2,j)
        wr2 = wr1*wr1-wi1*wi1
        wi2 = wr1*wi1+wr1*wi1
        wr3 = wr1*wr2-wi1*wi2
        wi3 = wr1*wi2+wi1*wr2
        wr4 = wr2*wr2-wi2*wi2
        wi4 = wr2*wi2+wr2*wi2
        x0 = a(1,j,2)+a(1,j,5)
        y0 = a(2,j,2)+a(2,j,5)
        x1 = a(1,j,3)+a(1,j,4)
        y1 = a(2,j,3)+a(2,j,4)
        x2 = c51*(a(1,j,2)-a(1,j,5))
        y2 = c51*(a(2,j,2)-a(2,j,5))
        x3 = c51*(a(1,j,3)-a(1,j,4))
        y3 = c51*(a(2,j,3)-a(2,j,4))
        x4 = x0+x1
        y4 = y0+y1
        x5 = c53*(x0-x1)
        y5 = c53*(y0-y1)
        x6 = a(1,j,1)-c54*x4
        y6 = a(2,j,1)-c54*y4
        x7 = x6+x5
        y7 = y6+y5
        x8 = x6-x5
        y8 = y6-y5
        x9 = y2+c52*y3
        y9 = -x2-c52*x3
        x10 = c52*y2-y3
        y10 = x3-c52*x2
        b(1,1,j) = a(1,j,1)+x4
        b(2,1,j) = a(2,j,1)+y4
        b(1,2,j) = wr1*(x7+x9)-wi1*(y7+y9)
        b(2,2,j) = wr1*(y7+y9)+wi1*(x7+x9)
        b(1,3,j) = wr2*(x8+x10)-wi2*(y8+y10)
        b(2,3,j) = wr2*(y8+y10)+wi2*(x8+x10)
        b(1,4,j) = wr3*(x8-x10)-wi3*(y8-y10)
        b(2,4,j) = wr3*(y8-y10)+wi3*(x8-x10)
        b(1,5,j) = wr4*(x7-x9)-wi4*(y7-y9)
        b(2,5,j) = wr4*(y7-y9)+wi4*(x7-x9)
      end do

      return
      end subroutine fft5a


      subroutine fft5b(a,b,w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      R_PREC :: a(2,m,l,*),b(2,m,5,*),w(2,*)

      INTG_PREC :: i, j
      R_PREC :: wr1, wi1, wr2, wi2, wr3, wi3, wr4, wi4
      R_PREC :: x0, y0, x1, y1, x2, y2, x3, y3
      R_PREC :: x4, y4, x5, y5, x6, y6, x7, y7
      R_PREC :: x8, y8, x9, y9, x10, y10

      R_PREC :: c51, c52, c53, c54
      data c51/0.95105651629515357_fftkind/
      data c52/0.61803398874989485_fftkind/
      data c53/0.55901699437494742_fftkind/
      data c54/0.25_fftkind/

      do i = 1,m
        x0 = a(1,i,1,2)+a(1,i,1,5)
        y0 = a(2,i,1,2)+a(2,i,1,5)
        x1 = a(1,i,1,3)+a(1,i,1,4)
        y1 = a(2,i,1,3)+a(2,i,1,4)
        x2 = c51*(a(1,i,1,2)-a(1,i,1,5))
        y2 = c51*(a(2,i,1,2)-a(2,i,1,5))
        x3 = c51*(a(1,i,1,3)-a(1,i,1,4))
        y3 = c51*(a(2,i,1,3)-a(2,i,1,4))
        x4 = x0+x1
        y4 = y0+y1
        x5 = c53*(x0-x1)
        y5 = c53*(y0-y1)
        x6 = a(1,i,1,1)-c54*x4
        y6 = a(2,i,1,1)-c54*y4
        x7 = x6+x5
        y7 = y6+y5
        x8 = x6-x5
        y8 = y6-y5
        x9 = y2+c52*y3
        y9 = -x2-c52*x3
        x10 = c52*y2-y3
        y10 = x3-c52*x2
        b(1,i,1,1) = a(1,i,1,1)+x4
        b(2,i,1,1) = a(2,i,1,1)+y4
        b(1,i,2,1) = x7+x9
        b(2,i,2,1) = y7+y9
        b(1,i,3,1) = x8+x10
        b(2,i,3,1) = y8+y10
        b(1,i,4,1) = x8-x10
        b(2,i,4,1) = y8-y10
        b(1,i,5,1) = x7-x9
        b(2,i,5,1) = y7-y9
      end do

      do j = 2,l
        wr1 = w(1,j)
        wi1 = w(2,j)
        wr2 = wr1*wr1-wi1*wi1
        wi2 = wr1*wi1+wr1*wi1
        wr3 = wr1*wr2-wi1*wi2
        wi3 = wr1*wi2+wi1*wr2
        wr4 = wr2*wr2-wi2*wi2
        wi4 = wr2*wi2+wr2*wi2
        do i = 1,m
          x0 = a(1,i,j,2)+a(1,i,j,5)
          y0 = a(2,i,j,2)+a(2,i,j,5)
          x1 = a(1,i,j,3)+a(1,i,j,4)
          y1 = a(2,i,j,3)+a(2,i,j,4)
          x2 = c51*(a(1,i,j,2)-a(1,i,j,5))
          y2 = c51*(a(2,i,j,2)-a(2,i,j,5))
          x3 = c51*(a(1,i,j,3)-a(1,i,j,4))
          y3 = c51*(a(2,i,j,3)-a(2,i,j,4))
          x4 = x0+x1
          y4 = y0+y1
          x5 = c53*(x0-x1)
          y5 = c53*(y0-y1)
          x6 = a(1,i,j,1)-c54*x4
          y6 = a(2,i,j,1)-c54*y4
          x7 = x6+x5
          y7 = y6+y5
          x8 = x6-x5
          y8 = y6-y5
          x9 = y2+c52*y3
          y9 = -x2-c52*x3
          x10 = c52*y2-y3
          y10 = x3-c52*x2
          b(1,i,1,j) = a(1,i,j,1)+x4
          b(2,i,1,j) = a(2,i,j,1)+y4
          b(1,i,2,j) = wr1*(x7+x9)-wi1*(y7+y9)
          b(2,i,2,j) = wr1*(y7+y9)+wi1*(x7+x9)
          b(1,i,3,j) = wr2*(x8+x10)-wi2*(y8+y10)
          b(2,i,3,j) = wr2*(y8+y10)+wi2*(x8+x10)
          b(1,i,4,j) = wr3*(x8-x10)-wi3*(y8-y10)
          b(2,i,4,j) = wr3*(y8-y10)+wi3*(x8-x10)
          b(1,i,5,j) = wr4*(x7-x9)-wi4*(y7-y9)
          b(2,i,5,j) = wr4*(y7-y9)+wi4*(x7-x9)
        end do
      end do

      return
      end subroutine fft5b


      subroutine fft8a(a,b,w,l)

      use ffte_param

      implicit none

      INTG_PREC :: l
      R_PREC :: a(2,l,*),b(2,8,*),w(2,*)

      INTG_PREC :: j
      R_PREC :: wr1, wi1, wr2, wi2, wr3, wi3, wr4, wi4
      R_PREC :: wr5, wi5, wr6, wi6, wr7, wi7
      R_PREC :: x0, y0, x1, y1, x2, y2, x3, y3
      R_PREC :: x4, y4, x5, y5, x6, y6, x7, y7
      R_PREC :: u0, v0
      R_PREC :: u1, v1
      R_PREC :: u2, v2
      R_PREC :: u3, v3

      R_PREC :: c81
      data c81/0.70710678118654752_fftkind/

      do j = 1,l
        wr1 = w(1,j)
        wi1 = w(2,j)
        wr2 = wr1*wr1-wi1*wi1
        wi2 = wr1*wi1+wr1*wi1
        wr3 = wr1*wr2-wi1*wi2
        wi3 = wr1*wi2+wi1*wr2
        wr4 = wr2*wr2-wi2*wi2
        wi4 = wr2*wi2+wr2*wi2
        wr5 = wr2*wr3-wi2*wi3
        wi5 = wr2*wi3+wi2*wr3
        wr6 = wr3*wr3-wi3*wi3
        wi6 = wr3*wi3+wr3*wi3
        wr7 = wr3*wr4-wi3*wi4
        wi7 = wr3*wi4+wi3*wr4
        x0 = a(1,j,1)+a(1,j,5)
        y0 = a(2,j,1)+a(2,j,5)
        x1 = a(1,j,1)-a(1,j,5)
        y1 = a(2,j,1)-a(2,j,5)
        x2 = a(1,j,3)+a(1,j,7)
        y2 = a(2,j,3)+a(2,j,7)
        x3 = a(2,j,3)-a(2,j,7)
        y3 = a(1,j,7)-a(1,j,3)
        u0 = x0+x2
        v0 = y0+y2
        u1 = x0-x2
        v1 = y0-y2
        x4 = a(1,j,2)+a(1,j,6)
        y4 = a(2,j,2)+a(2,j,6)
        x5 = a(1,j,2)-a(1,j,6)
        y5 = a(2,j,2)-a(2,j,6)
        x6 = a(1,j,4)+a(1,j,8)
        y6 = a(2,j,4)+a(2,j,8)
        x7 = a(1,j,4)-a(1,j,8)
        y7 = a(2,j,4)-a(2,j,8)
        u2 = x4+x6
        v2 = y4+y6
        u3 = y4-y6
        v3 = x6-x4
        b(1,1,j) = u0+u2
        b(2,1,j) = v0+v2
        b(1,5,j) = wr4*(u0-u2)-wi4*(v0-v2)
        b(2,5,j) = wr4*(v0-v2)+wi4*(u0-u2)
        b(1,3,j) = wr2*(u1+u3)-wi2*(v1+v3)
        b(2,3,j) = wr2*(v1+v3)+wi2*(u1+u3)
        b(1,7,j) = wr6*(u1-u3)-wi6*(v1-v3)
        b(2,7,j) = wr6*(v1-v3)+wi6*(u1-u3)
        u0 = x1+c81*(x5-x7)
        v0 = y1+c81*(y5-y7)
        u1 = x1-c81*(x5-x7)
        v1 = y1-c81*(y5-y7)
        u2 = x3+c81*(y5+y7)
        v2 = y3-c81*(x5+x7)
        u3 = x3-c81*(y5+y7)
        v3 = y3+c81*(x5+x7)
        b(1,2,j) = wr1*(u0+u2)-wi1*(v0+v2)
        b(2,2,j) = wr1*(v0+v2)+wi1*(u0+u2)
        b(1,6,j) = wr5*(u1+u3)-wi5*(v1+v3)
        b(2,6,j) = wr5*(v1+v3)+wi5*(u1+u3)
        b(1,4,j) = wr3*(u1-u3)-wi3*(v1-v3)
        b(2,4,j) = wr3*(v1-v3)+wi3*(u1-u3)
        b(1,8,j) = wr7*(u0-u2)-wi7*(v0-v2)
        b(2,8,j) = wr7*(v0-v2)+wi7*(u0-u2)
      end do

      return
      end subroutine fft8a


      subroutine fft8b(a,b,w,m,l)

      use ffte_param

      implicit none

      INTG_PREC :: m, l
      R_PREC :: a(2,m,l,*),b(2,m,8,*),w(2,*)

      INTG_PREC :: i, j
      R_PREC :: wr1, wi1, wr2, wi2, wr3, wi3, wr4, wi4
      R_PREC :: wr5, wi5, wr6, wi6, wr7, wi7
      R_PREC :: x0, y0, x1, y1, x2, y2, x3, y3
      R_PREC :: x4, y4, x5, y5, x6, y6, x7, y7
      R_PREC :: u0, v0
      R_PREC :: u1, v1
      R_PREC :: u2, v2
      R_PREC :: u3, v3

      R_PREC :: c81
      data c81/0.70710678118654752_fftkind/

      do i = 1,m
        x0 = a(1,i,1,1)+a(1,i,1,5)
        y0 = a(2,i,1,1)+a(2,i,1,5)
        x1 = a(1,i,1,1)-a(1,i,1,5)
        y1 = a(2,i,1,1)-a(2,i,1,5)
        x2 = a(1,i,1,3)+a(1,i,1,7)
        y2 = a(2,i,1,3)+a(2,i,1,7)
        x3 = a(2,i,1,3)-a(2,i,1,7)
        y3 = a(1,i,1,7)-a(1,i,1,3)
        u0 = x0+x2
        v0 = y0+y2
        u1 = x0-x2
        v1 = y0-y2
        x4 = a(1,i,1,2)+a(1,i,1,6)
        y4 = a(2,i,1,2)+a(2,i,1,6)
        x5 = a(1,i,1,2)-a(1,i,1,6)
        y5 = a(2,i,1,2)-a(2,i,1,6)
        x6 = a(1,i,1,4)+a(1,i,1,8)
        y6 = a(2,i,1,4)+a(2,i,1,8)
        x7 = a(1,i,1,4)-a(1,i,1,8)
        y7 = a(2,i,1,4)-a(2,i,1,8)
        u2 = x4+x6
        v2 = y4+y6
        u3 = y4-y6
        v3 = x6-x4
        b(1,i,1,1) = u0+u2
        b(2,i,1,1) = v0+v2
        b(1,i,5,1) = u0-u2
        b(2,i,5,1) = v0-v2
        b(1,i,3,1) = u1+u3
        b(2,i,3,1) = v1+v3
        b(1,i,7,1) = u1-u3
        b(2,i,7,1) = v1-v3
        u0 = x1+c81*(x5-x7)
        v0 = y1+c81*(y5-y7)
        u1 = x1-c81*(x5-x7)
        v1 = y1-c81*(y5-y7)
        u2 = x3+c81*(y5+y7)
        v2 = y3-c81*(x5+x7)
        u3 = x3-c81*(y5+y7)
        v3 = y3+c81*(x5+x7)
        b(1,i,2,1) = u0+u2
        b(2,i,2,1) = v0+v2
        b(1,i,6,1) = u1+u3
        b(2,i,6,1) = v1+v3
        b(1,i,4,1) = u1-u3
        b(2,i,4,1) = v1-v3
        b(1,i,8,1) = u0-u2
        b(2,i,8,1) = v0-v2
      end do

      do j = 2,l
        wr1 = w(1,j)
        wi1 = w(2,j)
        wr2 = wr1*wr1-wi1*wi1
        wi2 = wr1*wi1+wr1*wi1
        wr3 = wr1*wr2-wi1*wi2
        wi3 = wr1*wi2+wi1*wr2
        wr4 = wr2*wr2-wi2*wi2
        wi4 = wr2*wi2+wr2*wi2
        wr5 = wr2*wr3-wi2*wi3
        wi5 = wr2*wi3+wi2*wr3
        wr6 = wr3*wr3-wi3*wi3
        wi6 = wr3*wi3+wr3*wi3
        wr7 = wr3*wr4-wi3*wi4
        wi7 = wr3*wi4+wi3*wr4
        do i = 1,m
          x0 = a(1,i,j,1)+a(1,i,j,5)
          y0 = a(2,i,j,1)+a(2,i,j,5)
          x1 = a(1,i,j,1)-a(1,i,j,5)
          y1 = a(2,i,j,1)-a(2,i,j,5)
          x2 = a(1,i,j,3)+a(1,i,j,7)
          y2 = a(2,i,j,3)+a(2,i,j,7)
          x3 = a(2,i,j,3)-a(2,i,j,7)
          y3 = a(1,i,j,7)-a(1,i,j,3)
          u0 = x0+x2
          v0 = y0+y2
          u1 = x0-x2
          v1 = y0-y2
          x4 = a(1,i,j,2)+a(1,i,j,6)
          y4 = a(2,i,j,2)+a(2,i,j,6)
          x5 = a(1,i,j,2)-a(1,i,j,6)
          y5 = a(2,i,j,2)-a(2,i,j,6)
          x6 = a(1,i,j,4)+a(1,i,j,8)
          y6 = a(2,i,j,4)+a(2,i,j,8)
          x7 = a(1,i,j,4)-a(1,i,j,8)
          y7 = a(2,i,j,4)-a(2,i,j,8)
          u2 = x4+x6
          v2 = y4+y6
          u3 = y4-y6
          v3 = x6-x4
          b(1,i,1,j) = u0+u2
          b(2,i,1,j) = v0+v2
          b(1,i,5,j) = wr4*(u0-u2)-wi4*(v0-v2)
          b(2,i,5,j) = wr4*(v0-v2)+wi4*(u0-u2)
          b(1,i,3,j) = wr2*(u1+u3)-wi2*(v1+v3)
          b(2,i,3,j) = wr2*(v1+v3)+wi2*(u1+u3)
          b(1,i,7,j) = wr6*(u1-u3)-wi6*(v1-v3)
          b(2,i,7,j) = wr6*(v1-v3)+wi6*(u1-u3)
          u0 = x1+c81*(x5-x7)
          v0 = y1+c81*(y5-y7)
          u1 = x1-c81*(x5-x7)
          v1 = y1-c81*(y5-y7)
          u2 = x3+c81*(y5+y7)
          v2 = y3-c81*(x5+x7)
          u3 = x3-c81*(y5+y7)
          v3 = y3+c81*(x5+x7)
          b(1,i,2,j) = wr1*(u0+u2)-wi1*(v0+v2)
          b(2,i,2,j) = wr1*(v0+v2)+wi1*(u0+u2)
          b(1,i,6,j) = wr5*(u1+u3)-wi5*(v1+v3)
          b(2,i,6,j) = wr5*(v1+v3)+wi5*(u1+u3)
          b(1,i,4,j) = wr3*(u1-u3)-wi3*(v1-v3)
          b(2,i,4,j) = wr3*(v1-v3)+wi3*(u1-u3)
          b(1,i,8,j) = wr7*(u0-u2)-wi7*(v0-v2)
          b(2,i,8,j) = wr7*(v0-v2)+wi7*(u0-u2)
        end do
      end do

      return
      end subroutine fft8b
