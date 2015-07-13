#include "../enzo/fortran.def"
      subroutine rotate3d(x,n1,n2,n3,y)

!     Simple blocked 3D rotation
!     Block size may be too small for some systems
!     Not worth worrying about unless the grids are large

!     Written by: Robert Harkness
!     Date:       November, 2003

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      INTG_PREC n1,n2,n3
      CMPLX_PREC x(n1,n2,n3), y(n2,n3,n1)

!     Locals

      INTG_PREC i,j,k
      INTG_PREC ii,jj,kk
      INTG_PREC ib,jb,kb
      INTG_PREC bs

!     Blocksize

      bs = 64

      ib=n1/max(n1/bs,1_IKIND)
      jb=n2/max(n2/bs,1_IKIND)
      kb=n3/max(n3/bs,1_IKIND)

!     do i=1,n1
!     do k=1,n3
!     do j=1,n2
!     y(j,k,i) = x(i,j,k)
!     end do
!     end do
!     end do
!     go to 666

      do ii=1,n1,ib
        do kk=1,n3,kb
          do jj=1,n2,jb

            do i=ii,min(ii+ib-1,n1)
              do k=kk,min(kk+kb-1,n3)
                do j=jj,min(jj+jb-1,n2)
                  y(j,k,i) = x(i,j,k)
                end do
              end do 
            end do

          end do
        end do
      end do

  666 continue

      return
      end
