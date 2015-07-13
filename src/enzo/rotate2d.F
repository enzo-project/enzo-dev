#include "fortran.def"
      subroutine rotate2d(x,n1,n2,y)

      implicit none
#include "fortran_types.def"

      INTG_PREC n1,n2
      CMPLX_PREC x(n1,n2), y(n2,n1)

      INTG_PREC i,j
      INTG_PREC ii,jj
      INTG_PREC ib,jb
      INTG_PREC bs

      bs = 64
      ib = n1/max(n1/bs,1)
      jb = n2/max(n2/bs,1)

!     do i=1,n1
!     do j=1,n2
!     y(j,i) = x(i,j)
!     end do
!     end do

!     go to 666

      do ii=1,n1,ib
        do jj=1,n2,jb

          do i=ii,min(ii+ib-1,n1)
            do j=jj,min(jj+jb-1,n2)
              y(j,i) = x(i,j)
            end do
          end do 

        end do
      end do

  666 continue

      return
      end
