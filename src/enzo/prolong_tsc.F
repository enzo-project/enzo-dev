c=======================================================================
c/////////////////////////  SUBROUTINE PROLONG_TSC  \\\\\\\\\\\\\\\\\\\\
c
      subroutine prolong_tsc(source, dest, ndim, sdim1, sdim2, sdim3,
     &                       ddim1, ddim2, ddim3, 
     &                       start1, start2, start3,
     &                       refine1, refine2, refine3)
c
c  PROLONG FROM SOURCE TO DEST
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source       - source field
c     sdim1-3      - source dimension
c     ddim1-3      - destination dimension
c     ndim         - rank of fields
c
c  OUTPUT ARGUMENTS: 
c     dest         - restricted field
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim,
     &        start1, start2, start3, refine1, refine2, refine3
      R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
c
c  locals
c
      INTG_PREC i, j, k, i1, j1, k1
      R_PREC    fact1, fact2, fact3, x, y, z, dxm, dym, dzm, dx0, dy0,
     &        dz0, dxp, dyp, dzp
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Precompute some things
c
      fact1 = 1._RKIND/REAL(refine1,RKIND)
      fact2 = 1._RKIND/REAL(refine2,RKIND)
      fact3 = 1._RKIND/REAL(refine3,RKIND)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=1, ddim1
            x = (REAL(i+start1,RKIND)-0.5_RKIND)*fact1 + 0.5_RKIND
            i1 = int(x,IKIND) + 1
            dxm = 0.5_RKIND*(        -x+REAL(i1,RKIND))
            dxp = 0.5_RKIND*(1._RKIND+x-REAL(i1,RKIND))
            dx0 = 1._RKIND - dxp - dxm
            dest(i,1,1) = source(i1-1,1,1)*dxm + 
     &                    source(i1  ,1,1)*dx0 +
     &                    source(i1+1,1,1)*dxp
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=1, ddim2
            y = (REAL(j+start2,RKIND)-0.5_RKIND)*fact2 + 0.5_RKIND
            j1 = int(y,IKIND) + 1
            dym = 0.5_RKIND*(        -y+REAL(j1,RKIND))
            dyp = 0.5_RKIND*(1._RKIND+y-REAL(j1,RKIND))
            dy0 = 1._RKIND - dyp - dym
            do i=1, ddim1
               x = (REAL(i+start1,RKIND)-0.5_RKIND)*fact1 + 0.5_RKIND
               i1 = int(x,IKIND) + 1
               dxm = 0.5_RKIND*(        -x+REAL(i1,RKIND))
               dxp = 0.5_RKIND*(1._RKIND+x-REAL(i1,RKIND))
               dx0 = 1._RKIND - dxp - dxm
               dest(i,j,1) = source(i1-1,j1-1,1)*dxm*dym + 
     &                       source(i1  ,j1-1,1)*dx0*dym +
     &                       source(i1+1,j1-1,1)*dxp*dym +
     &                       source(i1-1,j1  ,1)*dxm*dy0 + 
     &                       source(i1  ,j1  ,1)*dx0*dy0 +
     &                       source(i1+1,j1  ,1)*dxp*dy0 +
     &                       source(i1-1,j1+1,1)*dxm*dyp + 
     &                       source(i1  ,j1+1,1)*dx0*dyp +
     &                       source(i1+1,j1+1,1)*dxp*dyp
            enddo
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         do k=1, ddim3
            z = (REAL(k+start3,RKIND)-0.5_RKIND)*fact3 + 0.5_RKIND
            k1 = int(z,IKIND) + 1
            dzm = 0.5_RKIND*(        -z+REAL(k1,RKIND))**2
            dzp = 0.5_RKIND*(1._RKIND+z-REAL(k1,RKIND))**2
            dz0 = 1._RKIND - dzp - dzm
            do j=1, ddim2
               y = (REAL(j+start2,RKIND)-0.5_RKIND)*fact2 + 0.5_RKIND
               j1 = int(y,IKIND) + 1
               dym = 0.5_RKIND*(        -y+REAL(j1,RKIND))**2
               dyp = 0.5_RKIND*(1._RKIND+y-REAL(j1,RKIND))**2
               dy0 = 1._RKIND - dyp - dym
               do i=1, ddim1
                  x = (REAL(i+start1,RKIND)-0.5_RKIND)*fact1 + 0.5_RKIND
                  i1 = int(x,IKIND) + 1
                  dxm = 0.5_RKIND*(        -x+REAL(i1,RKIND))**2
                  dxp = 0.5_RKIND*(1._RKIND+x-REAL(i1,RKIND))**2
                  dx0 = 1._RKIND - dxp - dxm
                  dest(i,j,k) = source(i1-1,j1-1,k1-1)*dxm*dym*dzm + 
     &                          source(i1  ,j1-1,k1-1)*dx0*dym*dzm +
     &                          source(i1+1,j1-1,k1-1)*dxp*dym*dzm +
     &                          source(i1-1,j1  ,k1-1)*dxm*dy0*dzm + 
     &                          source(i1  ,j1  ,k1-1)*dx0*dy0*dzm +
     &                          source(i1+1,j1  ,k1-1)*dxp*dy0*dzm +
     &                          source(i1-1,j1+1,k1-1)*dxm*dyp*dzm + 
     &                          source(i1  ,j1+1,k1-1)*dx0*dyp*dzm +
     &                          source(i1+1,j1+1,k1-1)*dxp*dyp*dzm +

     &                          source(i1-1,j1-1,k1  )*dxm*dym*dz0 + 
     &                          source(i1  ,j1-1,k1  )*dx0*dym*dz0 +
     &                          source(i1+1,j1-1,k1  )*dxp*dym*dz0 +
     &                          source(i1-1,j1  ,k1  )*dxm*dy0*dz0 + 
     &                          source(i1  ,j1  ,k1  )*dx0*dy0*dz0 +
     &                          source(i1+1,j1  ,k1  )*dxp*dy0*dz0 +
     &                          source(i1-1,j1+1,k1  )*dxm*dyp*dz0 + 
     &                          source(i1  ,j1+1,k1  )*dx0*dyp*dz0 +
     &                          source(i1+1,j1+1,k1  )*dxp*dyp*dz0 +

     &                          source(i1-1,j1-1,k1+1)*dxm*dym*dzp + 
     &                          source(i1  ,j1-1,k1+1)*dx0*dym*dzp +
     &                          source(i1+1,j1-1,k1+1)*dxp*dym*dzp +
     &                          source(i1-1,j1  ,k1+1)*dxm*dy0*dzp + 
     &                          source(i1  ,j1  ,k1+1)*dx0*dy0*dzp +
     &                          source(i1+1,j1  ,k1+1)*dxp*dy0*dzp +
     &                          source(i1-1,j1+1,k1+1)*dxm*dyp*dzp + 
     &                          source(i1  ,j1+1,k1+1)*dx0*dyp*dzp +
     &                          source(i1+1,j1+1,k1+1)*dxp*dyp*dzp
               enddo
            enddo
         enddo
      endif
c
      return
      end
