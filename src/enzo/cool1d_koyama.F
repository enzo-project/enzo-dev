
#include "fortran.def"
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
c
c  From Koyama & Inutsuka 2006, arxiv: asrto-ph 0605528v1
c
c  Arguments
c
      integer in, jn, kn, is, ie, j, k,
     &        idual, nratec, idim, imethod, iter
      real    temstart, temend, utem, uxyz, urho, utim,
     &        gamma, coola(nratec)
      real    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn)
c
c  Parameters
c
      real pmin
      parameter (pmin = tiny)
      double precision mh

      parameter (mh = 1.67262d-24)  ! DPC

c
c  Locals
c
      integer i
      real thalf
c
c  Slice locals
c 
      real p2d(in), tgas(in), tgasold(in), cool(in)
      real edot(in), kboltz, Lambda, Heating
      parameter ( kboltz = 1.30865d-16 ) !erg / Kelvin

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
            p2d(i) = (gamma - 1.0)*d(i,j,k)*e(i,j,k)
         enddo
      else
         if (idual .eq. 1) then
c
c           PPM with dual energy -- use gas energy
c
            do i = is+1, ie+1
               p2d(i) = (gamma - 1.0)*d(i,j,k)*ge(i,j,k)
            enddo
         else
c
c           PPM without dual energy -- use total energy
c
            do i = is+1, ie+1
               p2d(i) = e(i,j,k) - 0.5*u(i,j,k)**2
               if (idim .gt. 1) p2d(i) = p2d(i) - 0.5*v(i,j,k)**2
               if (idim .gt. 2) p2d(i) = p2d(i) - 0.5*w(i,j,k)**2
               p2d(i) = max((gamma - 1.0)*d(i,j,k)*p2d(i), tiny)
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
        thalf = 0.5*(tgas(i) + tgasold(i))
c        thalf = tgas(i)                              !kludge
        edot(i) = 0
        if( thalf .gt. temstart ) then
            Heating = 2.d-26
            Lambda = -(1d7 * exp(-1.184d5/(thalf+1000)) +
     &          14d-3*sqrt(thalf)*exp(-92.0/thalf))
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
