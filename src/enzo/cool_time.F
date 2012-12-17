#include "fortran.def"
c=======================================================================
c//////////////////////  SUBROUTINE COOL_TIME \ \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cool_time(
     &                d, e, ge, u, v, w, cooltime,
     &                in, jn, kn, nratec, iexpand, imethod,
     &                idual, idim, igammah,
     &                is, js, ks, ie, je, ke, 
     &                dt, aye, temstart, temend, fh, 
     &                utem, urho,
     &                eta1, eta2, gamma, coola, gammaha, mu)
c
c  CALCULATE COOLING TIME IN CODE UNITS (NON-MULTI SPECIES VERSION)
c
c  written by: Greg Bryan
c  date:       March, 1997
c  modified1:
c
c  PURPOSE:
c    returns the cooling time in code units on the grid
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
      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, 
     &        idual, iexpand, nratec, idim, imethod, igammah
      R_PREC    dt, aye, temstart, temend, fh, utem, urho, gammaha,
     &        eta1, eta2, gamma, coola(nratec), mu
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        cooltime(in,jn,kn)
c
c  Parameters
c
      INTG_PREC ijk
      parameter (ijk = MAX_ANY_SINGLE_DIRECTION)
c
c  Locals
c
      INTG_PREC i, j, k
      R_PREC    energy
c
c  Slice locals
c 
      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &     dtit(ijk), ttot(ijk), edot(ijk), tgas(ijk), tgasold(ijk),
     &     cool(ijk)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Loop over slices (in the k-direction)
c
      do k = ks+1, ke+1
       do j = js+1, je+1
c
c        Set time elapsed to zero for each cell (and convert dens to abs)
c
         do i = is+1, ie+1
            d(i,j,k) = d(i,j,k)/aye**3
         enddo
c
c        Compute the cooling rate on a slice
c
         call cool1d(d, e, ge, u, v, w,
     &               in, jn, kn, nratec, idual, idim, imethod, 
     &               1_IKIND, igammah,
     &               is, ie, j, k, 
     &               temstart, temend, fh, utem, urho, aye,
     &               eta1, eta2, gamma, coola, gammaha,
     &               indixe, t1, t2, logtem, tdef, edot,
     &               tgas, tgasold, p2d, cool, mu)
c
c        Compute the cooling time on the slice
c
         do i = is+1, ie+1
            energy = max(p2d(i)/d(i,j,k)/(gamma-1._RKIND), tiny)
            cooltime(i,j,k) = abs(energy/edot(i))
         enddo
c
c       Set the density back to comoving
c
         do i = is+1, ie+1
            d(i,j,k) = d(i,j,k)*aye**3
         enddo
c
c     Next slice
c
        enddo
      enddo
c
      return
      end

