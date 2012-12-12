#include "../enzo/fortran.def"
      subroutine spline(x,y,n,yp1,ypn,y2)

      implicit none
#include "../enzo/fortran_types.def"

!     Parameters

      INTG_PREC, parameter :: nmax = 1000

!     Arguments

      INTG_PREC :: n
      R_PREC :: x(n), y(n), y2(n)
      R_PREC :: yp1, ypn

!     Locals

      R_PREC :: u(nmax)
      R_PREC :: qn, un, p, sig
      INTG_PREC :: i, k

      if (yp1 .gt. 0.99e30_RKIND) then
        y2(1) = 0.0_RKIND
        u(1) = 0.0_RKIND
      else
        y2(1) = -0.5_RKIND
        u(1) = (3.0_RKIND/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i = 2, n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2.0_RKIND
        y2(i) = (sig-1.0_RKIND)/p
        u(i) = (6.0_RKIND*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if (ypn .gt. 0.99e30_RKIND) then
        qn = 0.0_RKIND
        un = 0.0_RKIND
      else
        qn = 0.5_RKIND
        un = (3.0_RKIND/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0_RKIND)

      do k = n-1, 1, -1
        y2(k) = y2(k)*y2(k+1)+u(k)
      end do

      return
      end



      subroutine splint(xa,ya,y2a,n,x,y)

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      INTG_PREC :: n
      R_PREC :: xa(n), ya(n), y2a(n)
      R_PREC :: x, y

!     Locals

      INTG_PREC :: klo, khi, k
      R_PREC :: h, a, b

      klo = 1
      khi = n

    1 continue
      if (khi-klo > 1) then
        k = (khi+klo)/2
        if(xa(k) > x)then
          khi = k
        else
          klo = k
        endif
      goto 1
      endif

      h = xa(khi)-xa(klo)

      if (h == 0.0_RKIND) stop 'bad xa input.'

      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+
     &      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_RKIND

      return
      end
