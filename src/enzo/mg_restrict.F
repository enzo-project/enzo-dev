c=======================================================================
c/////////////////////////  SUBROUTINE MG_RESTRICT  \\\\\\\\\\\\\\\\\\\\
c
      subroutine mg_restrict(source, dest, ndim, sdim1, sdim2, sdim3,
     &                       ddim1, ddim2, ddim3)
c
c  MULTIGRID: RESTRICT FROM SOURCE TO DEST
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
      integer ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim
      real    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
c
c  locals
c
      integer i, j, k, i1, j1, k1
      real    fact1, fact2, fact3, x, y, z, dxm, dym, dzm, dx0, dy0,
     &        dz0, dxp, dyp, dzp, coef1, coef3
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     These coefficients modify the restricted values for 1 and 2D problems
c       to speed up convergence.
c
      coef1 = 2.d0
      coef3 = 0.52d0
c
c     Precompute some things
c
      fact1 = real(sdim1-1)/real(ddim1-1)
      if (ndim .ge. 2) fact2 = real(sdim2-1)/real(ddim2-1)
      if (ndim .ge. 3) fact3 = real(sdim3-1)/real(ddim3-1)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=2, ddim1-1
            x = real(i-1)*fact1 + 0.5d0
            i1 = int(x) + 1
            dxm = 0.5d0*(   -x+real(i1))
            dxp = 0.5d0*(1.d0+x-real(i1))
            dx0 = 1.d0 - dxp - dxm
            dest(i,1,1) = source(i1-1,1,1)*dxm + 
     &                    source(i1  ,1,1)*dx0 +
     &                    source(i1+1,1,1)*dxp
            dest(i,1,1) = coef1*dest(i,1,1)
         enddo
         dest(    1,1,1) = source(    1,1,1)
         dest(ddim1,1,1) = source(sdim1,1,1)
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=2, ddim2-1
            y = real(j-1)*fact2 + 0.5d0
            j1 = int(y) + 1
            dym = 0.5d0*(   -y+real(j1))
            dyp = 0.5d0*(1.d0+y-real(j1))
            dy0 = 1.d0 - dyp - dym
            do i=2, ddim1-1
               x = real(i-1)*fact1 + 0.5d0
               i1 = int(x) + 1
               dxm = 0.5d0*(   -x+real(i1))
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
            dest(1    ,j,1) = source(1    ,j1,1)
            dest(ddim1,j,1) = source(sdim1,j1,1)
         enddo
         do i=1, ddim1
            i1 = min(max(int(real(i-1)*fact1 + 0.5d0) + 1, 1), sdim1)
            dest(i,    1,1) = source(i1,    1,1)
            dest(i,ddim2,1) = source(i1,sdim2,1)
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         do k=2, ddim3-1
            z = real(k-1)*fact3 + 0.5d0
            k1 = int(z) + 1
            dzm = 0.5d0*(   -z+real(k1))**2
            dzp = 0.5d0*(1.d0+z-real(k1))**2
            dz0 = 1.d0 - dzp - dzm
            do j=2, ddim2-1
               y = real(j-1)*fact2 + 0.5d0
               j1 = int(y) + 1
               dym = 0.5d0*(   -y+real(j1))**2
               dyp = 0.5d0*(1.d0+y-real(j1))**2
               dy0 = 1.d0 - dyp - dym
               do i=2, ddim1-1
                  x = real(i-1)*fact1 + 0.5d0
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
                  dest(i,j,k) = coef3*dest(i,j,k)
               enddo
               dest(1    ,j,k) = source(1    ,j1,k1)
               dest(ddim1,j,k) = source(sdim1,j1,k1)
            enddo
            do i=1, ddim1
               i1 = min(max(int(real(i-1)*fact1 + 0.5d0) + 1, 1), sdim1)
               dest(i,    1,k) = source(i1,    1,k1)
               dest(i,ddim2,k) = source(i1,sdim2,k1)
            enddo
         enddo
         do j=1, ddim2
            j1 = min(max(int(real(j-1)*fact2 + 0.5d0) + 1, 1), sdim2)
            do i=1, ddim1
               i1 = min(max(int(real(i-1)*fact1 + 0.5d0) + 1, 1), sdim1)
               dest(i,j,    1) = source(i1,j1,    1)
               dest(i,j,ddim3) = source(i1,j1,sdim3)
            enddo
         enddo
      endif
c
      return
      end
