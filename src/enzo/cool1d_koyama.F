#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c////////////////////////  SUBROUTINE cool1d_koyama  \\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cool1d_koyama(d,e,ge,u,v,w ,
     &                  in,jn, kn, nratec, idual, idim, imethod, iter,
     &                  is, ie, j,k, temstart, temend,
     &                  utem, uxyz, urho, utim, gamma, coola,
     &                  edot, tgas, tgasold, p2d, cool
     &                  )

    
      implicit NONE
#include "fortran_types.def"
c
c  From Koyama & Inutsuka 2006, arxiv: asrto-ph 0605528v1
c
c  Arguments
c
      INTG_PREC in, jn, kn, is, ie, j, k,
     &        idual, nratec, idim, imethod, iter
      R_PREC    temstart, temend, utem, uxyz, urho, utim,
     &        gamma, coola(nratec)
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
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
      R_PREC thalf
c
c  Slice locals
c 
      R_PREC p2d(in), tgas(in), tgasold(in), cool(in)
      R_PREC edot(in), Lambda, Heating

c
c     Note that the temperature definition here is slightly different
c     than that in cool1d.src.
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
         tgas(i) = max(p2d(i)*utem/d(i,j,k), temstart)
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
        thalf = 0.5_RKIND*(tgas(i) + tgasold(i))
c        thalf = tgas(i)                              !kludge
        edot(i) = 0._RKIND
        if( thalf .gt. temstart ) then
            Heating = 2.e-26_RKIND
            Lambda = -(1.e7_RKIND 
     &           * exp(-1.184e5_RKIND/(thalf+1000._RKIND)) +
     &           14.e-3_RKIND*sqrt(thalf)*exp(-92._RKIND/thalf))
            Lambda = Lambda * Heating
            edot(i) = d(i,j,k)*urho/(mh*mh)*Lambda
            edot(i) = edot(i) + Heating/mh
            edot(i) = edot(i) * utim /(utem * kboltz)*mh

        endif


      enddo

      do i=is+1, ie+1
         tgasold(i) = tgas(i)
      enddo
      end
c
