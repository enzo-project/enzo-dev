#include "fortran.def"
#include "error.def"
c=======================================================================
c///////////////////////  SUBROUTINE INTPRIM  \\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine intprim(dslice, uslice, vslice, wslice, pslice, 
     &     idim, i1, i2, isteep, steepen, iflatten, flatten, 
     &     iposrec, c1, c2, c3, c4, c5, c6, char1, char2, c0, 
     &     dd, dl, dr, d6, dla, dra, dl0, dr0,
     &     du, ul, ur, u6, ula, ura, ul0, ur0,
     &     dv, vl, vr, v6, vla, vra, vl0, vr0,
     &     dw, wl, wr, w6, wla, wra, wl0, wr0,
     &     dp, pl, pr, p6, pla, pra, pl0, pr0,
     &     lem, rem, dxdt2, gamma)
c
c  COMPUTES LEFT AND RIGHT EULERIAN INTERFACE VALUES FOR RIEMANN SOLVER
c
c  written by: Greg Bryan
c  date:       March, 1996
c  modified1:  June, 2010 by John Wise -- differencing and monotonizing 
c              done in characteristic variables.  Works on all primitive
c              variables.  Modified intvar.src.
c
c  PURPOSE:  Uses piecewise parabolic interpolation to compute left-
c    and right interface values to be fed into Riemann solver during a
c    one dimensional sweeps.  The routine works on all primitive 
c    variables in one dimension.
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
#include "fortran_types.def"
c
c  argument declarations
c
      INTG_PREC idim, i1, i2, iflatten, isteep, iposrec
      R_PREC gamma
      R_PREC c1(idim), c2(idim), c3(idim), c4(idim), c5(idim), c6(idim),
     &     char1(idim), char2(idim), c0(idim)
      R_PREC dslice(idim), uslice(idim), vslice(idim), wslice(idim),
     $     pslice(idim)
      R_PREC dla(idim), dra(idim), dl0(idim), dr0(idim),
     $     ula(idim), ura(idim), ul0(idim), ur0(idim),
     $     vla(idim), vra(idim), vl0(idim), vr0(idim),
     $     wla(idim), wra(idim), wl0(idim), wr0(idim),
     $     pla(idim), pra(idim), pl0(idim), pr0(idim)
      R_PREC dd(idim), dl(idim), dr(idim), d6(idim),
     $     du(idim), ul(idim), ur(idim), u6(idim),
     $     dv(idim), vl(idim), vr(idim), v6(idim),
     $     dw(idim), wl(idim), wr(idim), w6(idim),
     $     dp(idim), pl(idim), pr(idim), p6(idim)
      R_PREC steepen(idim), flatten(idim)
      R_PREC lem(idim,5,5), rem(idim,5,5), dxdt2(idim)
c
c  local declarations (arrays passed as temps)
c
      INTG_PREC i, m, n
      R_PREC temp1, temp2, dxdt, gamma1, pc, rhoc, du2, dplim, 
     $       limit, one
      R_PREC qplus(idim,5), qmnus(idim,5), qcent(idim,5), qvanl(idim,5),
     $     qplusc(idim,5), qmnusc(idim,5), qcentc(idim,5), 
     $     qvanlc(idim,5), dqc(idim,5)
      parameter (one=1._RKIND)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c     Compute average linear slopes (eqn 1.7)
c      JHW (June 2010): Added van Leer slopes (idea taken ATHENA)
c
      call linslope(dslice, c1, c2, i1, i2, idim, qplus(1,1),
     $     qmnus(1,1), qcent(1,1), qvanl(1,1))
      call linslope(uslice, c1, c2, i1, i2, idim, qplus(1,2),
     $     qmnus(1,2), qcent(1,2), qvanl(1,2))
      call linslope(vslice, c1, c2, i1, i2, idim, qplus(1,3),
     $     qmnus(1,3), qcent(1,3), qvanl(1,3))
      call linslope(wslice, c1, c2, i1, i2, idim, qplus(1,4),
     $     qmnus(1,4), qcent(1,4), qvanl(1,4))
      call linslope(pslice, c1, c2, i1, i2, idim, qplus(1,5),
     $     qmnus(1,5), qcent(1,5), qvanl(1,5))
c
c     Project primitive variable differences along characteristics
c     
      do n = 1, 5
         do i = i1-2, i2+2
            qplusc(i,n) = qplus(i,1) * lem(i,1,n)
            qmnusc(i,n) = qmnus(i,1) * lem(i,1,n)
            qcentc(i,n) = qcent(i,1) * lem(i,1,n)
            qvanlc(i,n) = qvanl(i,1) * lem(i,1,n)
         enddo
      enddo
      do n = 1, 5
         do m = 2, 5
            do i = i1-2, i2+2
               qplusc(i,n) = qplusc(i,n) + qplus(i,m) * lem(i,m,n)
               qmnusc(i,n) = qmnusc(i,n) + qmnus(i,m) * lem(i,m,n)
               qcentc(i,n) = qcentc(i,n) + qcent(i,m) * lem(i,m,n)
               qvanlc(i,n) = qvanlc(i,n) + qvanl(i,m) * lem(i,m,n)
            enddo
         enddo
      enddo
c
c     Monotonize characteristic variables (eqn 1.8)
c     
      do n = 1, 5
         do i = i1-2, i2+2
            if (qmnusc(i,n)*qplusc(i,n) .gt. 0.0) then
               temp1 = min(abs(qmnusc(i,n)), abs(qplusc(i,n)))
               temp2 = min(abs(qcentc(i,n)), abs(qvanlc(i,n)))
               dqc(i,n) = min(2._RKIND*temp1, temp2) * 
     &              sign(one, qcentc(i,n))
            else
               dqc(i,n) = 0._RKIND
            endif
         enddo
      enddo
c
c     If requested, limit the slopes to keep a positive pressure
c     (U-reconstruction; Waagan 2009, JCP, 228, 8609, Eqn 3.17).
c     Recall that i = 0->4 = rho, u, v, w, p.
c
      if (iposrec .eq. 1) then
         gamma1 = gamma - 1._RKIND
         ! Density
         do i = i1-2, i2+2
            limit = 0.5_RKIND * dslice(i)
            if (dqc(i,1) .gt. limit) then
               print*, 'Enforcing positive pressure by limiting ' ,
     $              'density slope', i, dslice(i), limit
               WARNING_MESSAGE
               dqc(i,1) = limit
            endif
         enddo
         ! Velocity1
         do i = i1-2, i2+2
            if (dqc(i,2) .gt. dxdt2(i)) then
               print*, 'Enforcing positive pressure by limiting ' ,
     $              'velocity slope', i, uslice(i), dxdt2(i)
               WARNING_MESSAGE
               dqc(i,2) = dxdt2(i)
            endif
         enddo
         ! Pressure
         do i = i1-2, i2+2
            du2 = dqc(i,2)**2
            do n = 3, 4
               du2 = du2 + dqc(i,n)**2
            enddo
            rhoc = dslice(i) - dqc(i,1)
            pc = pslice(i) - dqc(i,5)
            dplim = 2._RKIND * abs(0.125_RKIND*rhoc*du2*gamma1+pc)
            if (dqc(i,5) .gt. dplim) then
               print*, 'Enforcing positive pressure by limiting ' ,
     $              'pressure slope', i, pslice(i), dplim
               WARNING_MESSAGE
               dqc(i,5) = dplim
            endif
         enddo
      endif
c
c     Project monotonized slopes back to the primitive variables
c
      call char2prim(1_IKIND, dqc, rem, i1, i2, idim, dd)
      call char2prim(2_IKIND, dqc, rem, i1, i2, idim, du)
      call char2prim(3_IKIND, dqc, rem, i1, i2, idim, dv)
      call char2prim(4_IKIND, dqc, rem, i1, i2, idim, dw)
      call char2prim(5_IKIND, dqc, rem, i1, i2, idim, dp)
c
c     Construct left and right values (eqn 1.6)
c
      call lr_states(dslice, c3, c4, c5, c6, dd, i1, i2, idim, dl, dr)
      call lr_states(uslice, c3, c4, c5, c6, du, i1, i2, idim, ul, ur)
      call lr_states(vslice, c3, c4, c5, c6, dv, i1, i2, idim, vl, vr)
      call lr_states(wslice, c3, c4, c5, c6, dw, i1, i2, idim, wl, wr)
      call lr_states(pslice, c3, c4, c5, c6, dp, i1, i2, idim, pl, pr)
c
c     Steepen density if asked for (use precomputed steepening parameter)
c
      if (isteep .ne. 0) then
         do i = i1-1, i2+1
            dl(i) = (1._RKIND-steepen(i))*dl(i) + 
     &              steepen(i)*(dslice(i-1)+0.5_RKIND*dd(i-1))
            dr(i) = (1._RKIND-steepen(i))*dr(i) + 
     &              steepen(i)*(dslice(i+1)-0.5_RKIND*dd(i+1))
         enddo
      endif
c
c     Monotonize again, flatten, and check for boundness
c
      call mono_n_flat(dslice, dl, dr, iflatten, flatten, i1, i2, idim)
      call mono_n_flat(uslice, ul, ur, iflatten, flatten, i1, i2, idim)
      call mono_n_flat(vslice, vl, vr, iflatten, flatten, i1, i2, idim)
      call mono_n_flat(wslice, wl, wr, iflatten, flatten, i1, i2, idim)
      call mono_n_flat(pslice, pl, pr, iflatten, flatten, i1, i2, idim)
c
c    Now construct left and right interface values (eqn 1.12 and 3.3)
c
      call lr_interface(dslice, dl, dr, char1, char2, c0, i1, i2, idim,
     $     dd, d6, dla, dra, dl0, dr0)
      call lr_interface(uslice, ul, ur, char1, char2, c0, i1, i2, idim,
     $     du, u6, ula, ura, ul0, ur0)
      call lr_interface(vslice, vl, vr, char1, char2, c0, i1, i2, idim,
     $     dv, v6, vla, vra, vl0, vr0)
      call lr_interface(wslice, wl, wr, char1, char2, c0, i1, i2, idim,
     $     dw, w6, wla, wra, wl0, wr0)
      call lr_interface(pslice, pl, pr, char1, char2, c0, i1, i2, idim,
     $     dp, p6, pla, pra, pl0, pr0)
c
      return
      end

c-----------------------------------------------------------------------

      subroutine linslope(qslice, c1, c2, i1, i2, idim, qplus, qmnus, 
     $     qcent, qvanl)
c      
      implicit none
#include "fortran_types.def"
c
      INTG_PREC i1, i2, idim
      R_PREC qslice(idim), qplus(idim), qmnus(idim), qcent(idim), 
     $     qvanl(idim), c1(idim), c2(idim)
      
      INTG_PREC i

      do i = i1-2, i2+2
         qplus(i) = qslice(i+1)-qslice(i  )
         qmnus(i) = qslice(i  )-qslice(i-1)
         qcent(i) = c1(i)*qplus(i) + c2(i)*qmnus(i)
         if (qplus(i)*qmnus(i) .gt. 0) then
            qvanl(i) = 2._RKIND*qplus(i)*qmnus(i)/(qmnus(i)+qplus(i))
         else
            qvanl(i) = 0._RKIND
         endif
      enddo

      end

c-----------------------------------------------------------------------

      subroutine char2prim(n, dqc, rem, i1, i2, idim, dq)
c      
      implicit none
#include "fortran_types.def"
c
      INTG_PREC n, i1, i2, idim
      R_PREC dqc(idim,5), dq(idim), rem(idim,5,5)
      
      INTG_PREC i, m

      do i = i1-2, i2+2
         dq(i) = dqc(i,1) * rem(i,1,n)
         do m = 2, 5
            dq(i) = dq(i) + dqc(i,m) * rem(i,m,n)
         enddo
      enddo

      end

c-----------------------------------------------------------------------

      subroutine lr_states(qslice, c3, c4, c5, c6, dq, i1, i2, idim, 
     $     ql, qr)
c      
      implicit none
#include "fortran_types.def"
c
      INTG_PREC i1, i2, idim
      R_PREC qslice(idim), c3(idim), c4(idim), c5(idim), c6(idim)
      R_PREC dq(idim), ql(idim), qr(idim)

      INTG_PREC i

      do i = i1-1, i2+2
         ql(i) = c3(i)*qslice(i-1) + c4(i)*qslice(i) +
     &           c5(i)*    dq(i-1)   + c6(i)*dq(i)
         qr(i-1) = ql(i)
      enddo

      end

c-----------------------------------------------------------------------

      subroutine mono_n_flat(qslice, ql, qr, iflatten, flatten, 
     $     i1, i2, idim)
c
      implicit none
#include "fortran_types.def"
c
      INTG_PREC i1, i2, idim, iflatten
      R_PREC qslice(idim), ql(idim), qr(idim), flatten(idim)

      INTG_PREC i
      R_PREC temp1, temp2, temp3, temp22, temp23
c
c     Monotonize again (eqn 1.10)
c
      do i=i1-1,i2+1
         temp1 = (qr(i)-qslice(i))*(qslice(i)-ql(i))
         temp2 = qr(i)-ql(i)
         temp3 = 6._RKIND*(qslice(i)-0.5_RKIND*(qr(i)+ql(i)))
         if (temp1 .le. 0._RKIND) then
            ql(i) = qslice(i)
            qr(i) = qslice(i)
         endif
         temp22 = temp2**2
         temp23 = temp2*temp3
         if (temp22 .lt. temp23)
     &        ql(i) = 3._RKIND*qslice(i) - 2._RKIND*qr(i)
         if (temp22 .lt. -temp23)
     &        qr(i) = 3._RKIND*qslice(i) - 2._RKIND*ql(i)
      enddo
c
c     If requested, flatten slopes with flatteners calculated in calcdiss (4.1)
c
      if (iflatten .ne. 0) then
         do i = i1-1, i2+1
            ql(i) = qslice(i)*flatten(i) + ql(i)*(1._RKIND-flatten(i))
            qr(i) = qslice(i)*flatten(i) + qr(i)*(1._RKIND-flatten(i))
         enddo
      endif
c
c     Ensure that the L/R values lie between neighboring cell-centered 
c     values (Taken from ATHENA, lr_states)
c
#define CHECK_LR
#ifdef CHECK_LR
      do i = i1-1, i2+1
         ql(i) = max(min(qslice(i), qslice(i-1)), ql(i))
         ql(i) = min(max(qslice(i), qslice(i-1)), ql(i))
         qr(i) = max(min(qslice(i), qslice(i+1)), qr(i))
         qr(i) = min(max(qslice(i), qslice(i+1)), qr(i))
      enddo
#endif

      end

c-----------------------------------------------------------------------

      subroutine lr_interface(qslice, ql, qr, char1, char2, c0,
     $     i1, i2, idim, dq, q6, qla, qra, ql0, qr0)
c
      implicit none
#include "fortran_types.def"
c
      INTG_PREC i1, i2, idim
      R_PREC qslice(idim), ql(idim), qr(idim), char1(idim), char2(idim),
     $     c0(idim), dq(idim), q6(idim), qla(idim), qra(idim), 
     $     ql0(idim), qr0(idim)

      INTG_PREC i
      R_PREC ft
      parameter(ft = 4._RKIND/3._RKIND)

      do i = i1-1, i2+1
         q6(i) = 6._RKIND*(qslice(i)-0.5_RKIND*(ql(i)+qr(i)))
         dq(i) = qr(i) - ql(i)
      enddo
c
      do i = i1, i2+1
        qla(i) = qr(i-1)-char1(i-1)*(dq(i-1) - 
     &        (1._RKIND-ft*char1(i-1))*q6(i-1))
        qra(i) = ql(i  )+char2(i  )*(dq(i  ) + 
     &       (1._RKIND-ft*char2(i  ))*q6(i  ))
      enddo
c
      do i=i1,i2+1
         ql0(i)=qr(i-1)-c0(i-1)*(dq(i-1)-(1._RKIND-ft*c0(i-1))*q6(i-1))
         qr0(i)=ql(i  )-c0(i  )*(dq(i  )+(1._RKIND+ft*c0(i  ))*q6(i  ))
      enddo

      end
