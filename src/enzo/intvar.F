#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE INTVAR  \\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine intvar(qslice, idim, i1, i2, isteep, steepen, iflatten, 
     &                  flatten, c1, c2, c3, c4, c5, c6, char1, char2,
     &                  c0, dq, ql, qr, q6, qla, qra, ql0, qr0)
c
c  COMPUTES LEFT AND RIGHT EULERIAN INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Greg Bryan
c  date:       March, 1996
c  modified1:
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during a
c    one dimensional sweeps.  This version computes the Eulerian corrections
c    to the left and right states described in section three of Colella &
c    Woodward (1984), JCP.  The routine works on a single variable in
c    one dimension.
c
c  INPUT:
c    qslice   - one dimensional field of quantity q (one of d,e,u,v...)
c    idim     - declared dimension of 1D fields
c    i1, i2   - start and end indexes of active region
c    isteep   - steepening flag (1 = on, 0 = off); only apply to density!
c    steepen    - steepening coefficients
c    iflatten - flattening flag (1 = on, 0 = off)
c    flatten  - flattening coefficients
c    c1-6     - precomputed grid coefficients
c    char1,2  - characteristic distances for +/- waves (for average)
c    c0       - characteristic distance (for lagrangean cell face)
c    dq, ql, qr, q6 - 1D field temporaries
c    
c  OUTPUT:
c    qla, qra - left and right state values (from char1,2)
c    ql0, qr0 - left and right state values (from c0)
c
c  EXTERNALS:
c
c  LOCALS:
c
c  PARAMETERS:
c    ft     - a constant used in eq. 1.124 (=2*2/3)
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  argument declarations
c
      integer idim, i1, i2, iflatten, isteep
      real c1(idim), c2(idim), c3(idim), c4(idim), c5(idim), c6(idim),
     &     char1(idim), char2(idim), c0(idim),
     &     qla(idim), qra(idim), ql0(idim), qr0(idim)
      real qslice(idim), steepen(idim), flatten(idim)
c
c  parameters
c
      real ft
      parameter(ft = 4.0/3.0)
c
c  local declarations (arrays passed as temps)
c
      integer i
      real qplus, qmnus, temp1, temp2, temp3, temp22, temp23
      real dq(idim), ql(idim), qr(idim), q6(idim)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c     Compute average linear slopes (eqn 1.7)
c      Monotonize (eqn 1.8)
c
      do i = i1-2, i2+2
         qplus = qslice(i+1)-qslice(i  )
         qmnus = qslice(i  )-qslice(i-1)
         dq(i) = c1(i)*qplus + c2(i)*qmnus
         temp1 = min(abs(dq(i)), 2.0*abs(qmnus), 2.0*abs(qplus))
         if (qplus*qmnus .gt. 0) then
            dq(i) = temp1*sign(1.0, dq(i))
         else
            dq(i) = 0.0
         endif
      enddo
c     
c     Construct left and right values (eqn 1.6)
c
      do i = i1-1, i2+2
         ql(i) = c3(i)*qslice(i-1) + c4(i)*qslice(i) +
     &           c5(i)*    dq(i-1)   + c6(i)*dq(i)
         qr(i-1) = ql(i)
      enddo
c
c     Steepen if asked for (use precomputed steepening parameter)
c
      if (isteep .ne. 0) then
         do i = i1-1, i2+1
            ql(i) = (1.0-steepen(i))*ql(i) + 
     &              steepen(i)*(qslice(i-1)+0.5*dq(i-1))
            qr(i) = (1.0-steepen(i))*qr(i) + 
     &              steepen(i)*(qslice(i+1)-0.5*dq(i+1))
         enddo
      endif
c
c     Monotonize again (eqn 1.10)
c
      do i=i1-1,i2+1
         temp1 = (qr(i)-qslice(i))*(qslice(i)-ql(i))
         temp2 = qr(i)-ql(i)
         temp3 = 6.0*(qslice(i)-0.5*(qr(i)+ql(i)))
         if (temp1 .le. 0.0) then
            ql(i) = qslice(i)
            qr(i) = qslice(i)
         endif
         temp22 = temp2**2
         temp23 = temp2*temp3
         if (temp22 .lt. temp23)
     &        ql(i) = 3.0*qslice(i) - 2.0*qr(i)
         if (temp22 .lt. -temp23)
     &        qr(i) = 3.0*qslice(i) - 2.0*ql(i)
      enddo
c
c     If requested, flatten slopes with flatteners calculated in calcdiss (4.1)
c
      if (iflatten .ne. 0) then
         do i = i1-1, i2+1
            ql(i) = qslice(i)*flatten(i) + ql(i)*(1.0-flatten(i))
            qr(i) = qslice(i)*flatten(i) + qr(i)*(1.0-flatten(i))
         enddo
      endif
c
c    Now construct left and right interface values (eqn 1.12 and 3.3)
c
      do i = i1-1, i2+1
         q6(i) = 6.0*(qslice(i)-0.5*(ql(i)+qr(i)))
         dq(i) = qr(i) - ql(i)
      enddo
c
      do i = i1, i2+1
        qla(i)= qr(i-1)-char1(i-1)*(dq(i-1)-(1.0-ft*char1(i-1))*q6(i-1))
        qra(i)= ql(i  )+char2(i  )*(dq(i  )+(1.0-ft*char2(i  ))*q6(i  ))
      enddo
c
      do i=i1,i2+1
         ql0(i) = qr(i-1)-c0(i-1)*(dq(i-1)-(1.0-ft*c0(i-1))*q6(i-1))
         qr0(i) = ql(i  )-c0(i  )*(dq(i  )+(1.0+ft*c0(i  ))*q6(i  ))
      enddo
c
c
      return
      end

