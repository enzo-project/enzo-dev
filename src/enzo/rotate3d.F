      subroutine rotate3d(x,n1,n2,n3,y)

      implicit none

      integer n1,n2,n3
      complex x(n1,n2,n3), y(n2,n3,n1)

      integer i,j,k
      integer ii,jj,kk
      integer ib,jb,kb
      integer bs

      bs = 64

      ib=n1/max(n1/bs,1)
      jb=n2/max(n2/bs,1)
      kb=n3/max(n3/bs,1)

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
