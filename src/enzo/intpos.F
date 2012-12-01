#include "fortran.def"
#include "error.def"
c=======================================================================
c///////////////////////  SUBROUTINE INTPOS  \\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine intpos(dslice, uslice, pslice, idim, i1, i2,
     $     dxdt2, gamma,
     $     drho, rhol, rhor, dla, dra, dl0, dr0,
     $     du, ul, ur, ula, ura, ul0, ur0,
     $     dp, pl, pr, pla, pra, pl0, pr0)
c
c  ENFORCES POSTIVITY IN THE PRIMITIVE VARIABLES
c
c  written by: John Wise
c  date:       July, 2010
c  modified1:
c
c  PURPOSE:  Forces the pressure to be positive, using the method
c            described in Waagan (2009), JCP, 228, 8609.  Because we're
c            using the primitive variables, use W-reconstruction 
c            (equation 3.14).  Uses piecewise linear reconstruction.
c
c  INPUT:
c    dslice   - one dimensional field of density
c    uslice   - one dimensional field of velocity
c    pslice   - one dimensional field of pressure
c    idim     - declared dimension of 1D fields
c    i1, i2   - start and end indexes of active region
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
      INTG_PREC idim, i1, i2, iflatten, isteep
      R_PREC gamma
      R_PREC dslice(idim), uslice(idim), pslice(idim), dxdt2(idim),
     $     drho(idim), rhol(idim), rhor(idim), dla(idim), dra(idim),
     $     dl0(idim), dr0(idim),
     $     du(idim), ul(idim), ur(idim), ula(idim), ura(idim), 
     $     ul0(idim), ur0(idim),
     $     dp(idim), pl(idim), pr(idim), pla(idim), pra(idim), 
     $     pl0(idim), pr0(idim)

c
c  local declarations (arrays passed as temps)
c
      INTG_PREC i
      R_PREC limit, gamma1i
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c     Density
c     
      do i = i1, i1+1
         if (drho(i) .gt. 0.5_RKIND*dslice(i)) then
            drho(i) = 0.5_RKIND*dslice(i)
            dla(i) = dslice(i) - 0.5_RKIND*drho(i)
            dra(i) = dslice(i) + 0.5_RKIND*drho(i)
            dl0(i) = dla(i)
            dr0(i) = dra(i)
            print*, 'Enforcing positive pressure by limiting ' ,
     $           'density slope', i, dslice(i), drho(i)
            WARNING_MESSAGE
         endif
      enddo
c
c     Velocity1
c     
      gamma1i = 1._RKIND / (gamma + 1._RKIND)
      do i = i1, i1+1
         limit = 2._RKIND * gamma1i * dxdt2(i)
         if (du(i) .gt. limit) then
            du(i) = limit
            ula(i) = uslice(i) - 0.5_RKIND*du(i)
            ura(i) = uslice(i) + 0.5_RKIND*du(i)
            ul0(i) = ula(i)
            ur0(i) = ura(i)
            print*, 'Enforcing positive pressure by limiting ' ,
     $           'velocity slope', i, uslice(i), du(i)
            WARNING_MESSAGE
         endif
      enddo
c
c     Pressure
c     
      do i = i1, i1+1
         limit = pslice(i) * gamma1i
         if (du(i) .gt. limit) then
            dp(i) = limit
            pla(i) = pslice(i) - 0.5_RKIND*dp(i)
            pra(i) = pslice(i) + 0.5_RKIND*dp(i)
            pl0(i) = pla(i)
            pr0(i) = pra(i)
            print*, 'Enforcing positive pressure by limiting ' ,
     $           'pressure slope', i, pslice(i), dp(i)
            WARNING_MESSAGE
         endif
      enddo

      return
      end

