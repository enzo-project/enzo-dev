#include "../enzo/fortran.def"
!---------------------------------------------------------------------------
!  Routines from Numerical Recipes, unless I am mistaken (RH).
!---------------------------------------------------------------------------

      subroutine qromo(func,a,b,ss,choose)

!     Romberg integration on an open interval. returns as ss
!     the integral of the function func from a to b. using
!     any specified integrating sunroutine "choose" and
!     Romberg's' method. normally "choose" will be an open
!     formula, not evaluating the function at the endpoints.
!     it is assumed that "choose" triples the number of steps
!     on each call, and that its error series contains only
!     even powers of the number of steps. the routines
!     "midpnt", "midinf", "midsql", "midsqu", are possible
!     choices for "choose".

      implicit none
#include "../enzo/fortran_types.def"

!     Parameters

      INTG_PREC, parameter :: jmax=18, jmaxp=jmax+1, km=4, k=km+1
      R_PREC, parameter :: eps = 1.e-6_RKIND

!     Arguments

      R_PREC :: a, b, ss

!     Externals

      external :: func
      external :: choose

!     Locals

      R_PREC :: s(jmaxp), h(jmaxp)
      R_PREC :: dss
      INTG_PREC :: j

      h(1) = 1._RKIND

      do j = 1, jmax

        call choose(func,a,b,s(j),j)

        if (j.ge.k) then

          call polint(h(j-km),s(j-km),k,0._RKIND,ss,dss)
!         write(*,*) j,dss,ss,s(j)

          if (abs(dss) .lt. (eps*abs(ss)) ) return

        endif

        s(j+1) = s(j)

        h(j+1) = h(j)/9._RKIND

!     this is where the assumption of step tripling and an
!     even error series is used.

      end do

      stop 'too many steps.'

      end


      subroutine polint(xa,ya,n,x,y,dy)

!     this is a routine for polynomial interpolation or 
!     extrapolation. geven arrays xa and ya, each of lenth
!     n, and value x, this routine will return a value y,
!     and an error estimator dy. if p(x) is the polynomial
!     of degree n-1 such that p(xa_j) = ya_j, j = 1,...,n, then
!     the returned value y = p(x).

      implicit none
#include "../enzo/fortran_types.def"

!     Parameter

!     change nmax as desired to be the largest anticipated
!     value of n.

      INTG_PREC, parameter :: nmax = 10

!     Arguments

      INTG_PREC :: n
      R_PREC :: xa(n), ya(n)
      R_PREC :: x, y, dy

!     Locals

      R_PREC :: c(nmax), d(nmax)
      R_PREC :: dif, dift, ho, hp, w, den
      INTG_PREC :: ns
      INTG_PREC :: i, m

      ns = 1
      dif = abs(x-xa(1))

      do i = 1, n

!     here we find the index ns of the closest table entry.

        dift = abs(x-xa(i))
        if (dift.lt.dif) then
          ns = i
          dif = dift
        endif
        c(i) = ya(i)

!     and initialize the tableau of c's' and d's'.

        d(i) = ya(i)
      end do

      y = ya(ns)

!     this is the initial appromation to y.

      ns = ns-1

      do m = 1, n-1

!     for each column of the tableau,

        do i = 1, n-m

!     we loop over the current c's' and d's' and update them.

          ho = xa(i)-x
          hp = xa(i+m)-x
          w = c(i+1)-d(i)
          den = ho-hp
          if(den == 0._RKIND) stop 'den  =  0'

!     this error can only occur if two input xa's' are identical

          den = w/den
          d(i) = hp*den
          c(i) = ho*den

        end do

        if ((2*ns) < (n-m)) then

!     after each column in the tableau is completed, we
!     decided which correction , c or d, we want to add 
!     to our accumulating value of y, i.e. which path to
!     take through the tableau --- forking up or down.
!     we do this in such a way as to take the most 
!     "straight line" route through the tableau to its apex,
!     updating ns accordingly to keep track of where we are.
!     this route keeps the partial approximationa centered
!     (insofar as possible) on the target x. the last dy
!     added is thus the error indication.

          dy = c(ns+1)
        else
          dy = d(ns)
          ns = ns-1
        endif

        y = y+dy

      end do

      return
      end


      subroutine midinf(funk,aa,bb,s,n)


!     this subroutine is an exact replacement for "midpnt",
!     i.e. returns as "s" the n-th stage of refinement of the
!     integral of "func" from aa to bb, except that the function
!     is evaluated at evenly spaced points in 1/x rather than 
!     in x. this allows the upper limit bb to be as large and
!     positive as the computer allows, or the lower limit aa
!     to be as large and negative, but not both. aa and bb must
!     have the same sign.

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      INTG_PREC :: n
      R_PREC :: aa, bb, s

      R_PREC :: funk
      external :: funk

!     Locals

      R_PREC :: a, b
      R_PREC :: del, ddel, x, sum, tnm
      INTG_PREC :: j
      INTG_PREC :: it

!     Statement function

      R_PREC :: func

      func(x) = funk(1._RKIND/x)/x**2

      save it

      b = 1._RKIND/aa
      a = 1._RKIND/bb

      if (n == 1) then
        s = (b-a)*func(0.5_RKIND*(a+b))
        it = 1
      else
        tnm = it
        del = (b-a)/(3._RKIND*tnm)
        ddel = del+del
        x = a+0.5_RKIND*del
        sum = 0._RKIND

        do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
        end do

        s = (s+(b-a)*sum/tnm)/3._RKIND
        it = 3*it
      endif

      return
      end


      subroutine midpnt(func,a,b,s,n)


!     this routine computes the n-th stage of refinement of an 
!     extended midpoint rule. "func" is input as the name of the 
!     function to be integrated between limits a and b, also input.
!     when called with n = 1, the routine returns as "s" the crudest
!     estimate of integral. subsequent calls with n = 2,3,... (in that
!     sequential order) will improve the accuracy of "s" by adding
!     (2/3)*3^{n-1} additional interior points. "s" should not be
!     modified between sequential calls.

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      R_PREC :: a, b, s
      INTG_PREC :: n

      R_PREC :: func
      external :: func

!     Locals

      R_PREC :: del, ddel, tnm, sum, x
      INTG_PREC :: j
      INTG_PREC :: it

      save it


      if (n.eq.1) then
        s = (b-a)*func(0.5_RKIND*(a+b))
        it = 1

!     2*it points will be added on the next refinement

      else
        tnm = it
        del = (b-a)/(3._RKIND*tnm)
        ddel = del+del
        x = a+0.5_RKIND*del
        sum = 0._RKIND

        do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
        end do

        s = (s+(b-a)*sum/tnm)/3._RKIND

!     the new sum id combined with the old integral to give
!     a refined integral.

        it = 3*it
      endif

      return
      end


      subroutine midsql(funk,aa,bb,s,n)


!     this routine is an exact replacement for midpnt except that it
!     allows for an inverse square-root singularity int the integrand
!     at the lower limit aa

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      R_PREC :: aa, bb, ss, s
      INTG_PREC :: n

      R_PREC :: funk
      external :: funk

!     Locals

      R_PREC :: a, b
      R_PREC :: del, ddel, tnm, sum, x
      INTG_PREC :: j
      INTG_PREC :: it

!     Statement function

      R_PREC :: func

      func(x) = 2._RKIND*x*funk(aa+x**2)

      save it


      b = sqrt(bb-aa)
      a = 0._RKIND

      if (n.eq.1) then
        s = (b-a)*func(0.5_RKIND*(a+b))
        it = 1

!     2*it points will be added on the next refinement

      else
        tnm = it
        del = (b-a)/(3._RKIND*tnm)
        ddel = del+del
        x = a+0.5_RKIND*del
        sum = 0._RKIND

        do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
        end do

        s = (s+(b-a)*sum/tnm)/3._RKIND

!     the new sum id combined with the old integral to give
!     a refined integral.

        it = 3*it
      endif

      return
      end


      subroutine midsqu(funk,aa,bb,s,n)

c     this routine is an exact replacement for midpnt except that it
c     allows for an inverse square-root singularity int the integrand
c     at the upper limit bb

      implicit none
#include "../enzo/fortran_types.def"

!     Arguments

      R_PREC :: aa, bb, s
      INTG_PREC :: n

      R_PREC :: funk
      external :: funk

!     Locals

      R_PREC :: a, b
      R_PREC :: del, ddel, tnm, sum, x
      INTG_PREC :: j
      INTG_PREC :: it

!     Statement function

      R_PREC :: func

      func(x) = 2._RKIND*x*funk(bb-x**2)

      save it

      b = sqrt(bb-aa)
      a = 0.0_RKIND

      if (n.eq.1) then
        s = (b-a)*func(0.5_RKIND*(a+b))
        it = 1

!     2*it points will be added on the next refinement

      else
        tnm = it
        del = (b-a)/(3.0_RKIND*tnm)
        ddel = del+del
        x = a+0.5_RKIND*del
        sum = 0.0_RKIND

        do j = 1, it
          sum = sum+func(x)
          x = x+ddel
          sum = sum+func(x)
          x = x+del
        end do

        s = (s+(b-a)*sum/tnm)/3.0_RKIND

!     the new sum id combined with the old integral to give
!     a refined integral.

        it = 3*it
      endif

      return
      end
