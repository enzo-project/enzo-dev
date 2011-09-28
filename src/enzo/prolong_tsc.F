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
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim,
     &        start1, start2, start3, refine1, refine2, refine3
      real    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
c
c  locals
c
      integer i, j, k, i1, j1, k1
      real    fact1, fact2, fact3, x, y, z, dxm, dym, dzm, dx0, dy0,
     &        dz0, dxp, dyp, dzp
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Precompute some things
c
      fact1 = 1.d0/real(refine1)
      fact2 = 1.d0/real(refine2)
      fact3 = 1.d0/real(refine3)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=1, ddim1
            x = (real(i+start1)-0.5)*fact1 + 0.5d0
            i1 = int(x) + 1
            dxm = 0.5d0*(    -x+real(i1))
            dxp = 0.5d0*(1.d0+x-real(i1))
            dx0 = 1.d0 - dxp - dxm
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
            y = (real(j+start2)-0.5d0)*fact2 + 0.5d0
            j1 = int(y) + 1
            dym = 0.5d0*(    -y+real(j1))
            dyp = 0.5d0*(1.d0+y-real(j1))
            dy0 = 1.d0 - dyp - dym
            do i=1, ddim1
               x = (real(i+start1)-0.5d0)*fact1 + 0.5d0
               i1 = int(x) + 1
               dxm = 0.5d0*(    -x+real(i1))
               dxp = 0.5d0*(1.d0+x-real(i1))
               dx0 = 1.d0 - dxp - dxm
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
            z = (real(k+start3)-0.5)*fact3 + 0.5d0
            k1 = int(z) + 1
            dzm = 0.5d0*(    -z+real(k1))**2
            dzp = 0.5d0*(1.d0+z-real(k1))**2
            dz0 = 1.d0 - dzp - dzm
            do j=1, ddim2
               y = (real(j+start2)-0.5d0)*fact2 + 0.5d0
               j1 = int(y) + 1
               dym = 0.5d0*(    -y+real(j1))**2
               dyp = 0.5d0*(1.d0+y-real(j1))**2
               dy0 = 1.d0 - dyp - dym
               do i=1, ddim1
                  x = (real(i+start1)-0.5d0)*fact1 + 0.5d0
                  i1 = int(x) + 1
                  dxm = 0.5d0*(   -x+real(i1))**2
                  dxp = 0.5d0*(1.d0+x-real(i1))**2
                  dx0 = 1.d0 - dxp - dxm
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
