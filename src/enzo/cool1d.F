#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c////////////////////////  SUBROUTINE COOL1D  \\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cool1d(
     &                d, e, ge, u, v, w,
     &                in, jn, kn, nratec, idual, idim, imethod, 
     &                iter, igammah,
     &                is, ie, j, k, 
     &                temstart, temend, fh, utem, urho, aye,
     &                eta1, eta2, gamma, coola, gammaha,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d, cool, mu
     &                     )
c
c  COMPUTE RADIATIVE COOLING/HEATING RATE (DE/DT IN CODE UNITS)
c
c  written by: Greg Bryan
c  date:       March, 1997
c  modified1:
c
c  PURPOSE:
c    Computes the radiativen cooling rate for a 2d slice, based the
c      provided cooling curve.
c
c  INPUTS:
c    is,ie   - start and end indicies of active region (zero-based!)
c
c  PARAMETERS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c  Arguments
c
      INTG_PREC in, jn, kn, is, ie, j, k,
     &        idual, nratec, idim, imethod, iter, igammah
      R_PREC  temstart, temend, fh, utem, urho, aye, gammaha,
     &        eta1, eta2, gamma, coola(nratec), mu
      R_PREC  d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn)
c
c  Parameters
c
      R_PREC pmin
      parameter (pmin = tiny)
      real*8 mh

      parameter (mh = mass_h)  ! DPC

c
c  Locals
c
      INTG_PREC i
      R_PREC logtem0, logtem9, dlogtem, dom, dom_inv
c
c  Slice locals
c 
      INTG_PREC indixe(in)
      R_PREC t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &     tgas(in), tgasold(in), cool(in)
      R_PREC edot(in)
c
c  Set the mean molecular mass
c     (see also Grid_ComputeTemperatureField.C)
c
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Set log values of start and end of lookup tables
c
      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem= (log(temend) - log(temstart))/REAL(nratec-1,RKIND)

c
c     Set units
c
      dom = urho*(aye**3)/mh
      dom_inv  = 1.0/dom
c
c     Compute Pressure
c
      if (imethod .eq. 2) then
c
c        Zeus - e() is really gas energy
c
         do i = is+1, ie+1
            p2d(i) = (gamma - 1._RKIND)*d(i,j,k)*e(i,j,k)
         enddo
      else
         if (idual .eq. 1) then
c
c           PPM with dual energy -- use gas energy
c
            do i = is+1, ie+1
               p2d(i) = (gamma - 1._RKIND)*d(i,j,k)*ge(i,j,k)
            enddo
         else
c
c           PPM without dual energy -- use total energy
c
            do i = is+1, ie+1
               p2d(i) = e(i,j,k) - 0.5_RKIND*u(i,j,k)**2
               if (idim .gt. 1) p2d(i) = p2d(i) - 0.5_RKIND*v(i,j,k)**2
               if (idim .gt. 2) p2d(i) = p2d(i) - 0.5_RKIND*w(i,j,k)**2
               p2d(i) = max((gamma - 1._RKIND)*d(i,j,k)*p2d(i), tiny)
            enddo
         endif
      endif
c
c     Compute temperature
c
      do i = is+1, ie+1
         tgas(i) = max(p2d(i)*utem*mu/d(i,j,k), temstart)
      enddo
c
c     If this is the first time through, just set tgasold to tgas
c
      if (iter .eq. 1) then
         do i = is+1, ie+1
            tgasold(i) = tgas(i)
         enddo
      endif
c
c     Loop over a slice
c
      do i = is+1, ie+1
c
c        Compute log temperature and truncate if above/below table max/min
c
         logtem(i) = log(0.5_RKIND*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)
c
c        Compute index into the table and precompute parts of linear interp
c
         indixe(i) = min(nratec-1, max(1,
     &        int((logtem(i)-logtem0)/dlogtem,IKIND)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = t2(i) - t1(i)
c
c        Lookup cooling values and do a linear temperature in log(T)
c
         cool(i) = coola(indixe(i)) + (logtem(i) - t1(i))
     .         *(coola(indixe(i)+1) -coola(indixe(i)))/tdef(i)
c
c        Compute the cooling function (assuming completely ionized)
c              (the factors of mh have been incorporated into coolinit).
c              (rate of change of specific energy)
c
         edot(i) = -cool(i)*(0.4_RKIND*(fh+1)*d(i,j,k))
c
c        Photoelectric heating  
c
         edot(i) = edot(i) + REAL(igammah,RKIND)*gammaha*dom_inv
c
      enddo
c
c     Compute (external) radiative heating terms
c
#ifdef RADIATION
#ifdef UNUSED
c
      do i = is+1, ie+1
         edot(i) = edot(i) + (
c
     .                  0._RKIND
c
     .      )/dom
      enddo
c
#endif /* UNUSED */
#endif /* RADIATION */
c
c     Set tgasold
c
      do i=is+1, ie+1
         tgasold(i) = tgas(i)
      enddo
c
      return
      end

