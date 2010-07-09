c=======================================================================
c/////////////////////////  SUBROUTINE MG_PROLONG  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine mg_prolong(source, dest, ndim, sdim1, sdim2, sdim3,
     &                      ddim1, ddim2, ddim3)
c
c  MULTIGRID: PROLONG FROM SOURCE TO DEST
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
c     dest         - prolonged field
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
      real    fact1, fact2, fact3, x, y, z, dx, dy, dz, 
     &        edge1, edge2, edge3, half
      parameter (half = 0.5001)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Precompute some things
c
      fact1 = real(sdim1-1)/real(ddim1-1)
      if (ndim .ge. 2) fact2 = real(sdim2-1)/real(ddim2-1)
      if (ndim .ge. 3) fact3 = real(sdim3-1)/real(ddim3-1)
      edge1 = real(sdim1) - half
      edge2 = real(sdim2) - half
      edge3 = real(sdim3) - half
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=1, ddim1
            x = min(max(real(i-1)*fact1 + 0.5, half), edge1)
            i1 = int(x + 0.5)
            dx = real(i1) + 0.5 - x
            dest(i,1,1) = source(i1,1,1)*dx + source(i1+1,1,1)*(1.0-dx)
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=1, ddim2
            y = min(max(real(j-1)*fact2 + 0.5, half), edge2)
            j1 = int(y + 0.5)
            dy = real(j1) + 0.5 - y
            do i=1, ddim1
               x = min(max(real(i-1)*fact1 + 0.5, half), edge1)
               i1 = int(x + 0.5)
               dx = real(i1) + 0.5 - x
               dest(i,j,1) = source(i1  ,j1  ,1)*     dx *     dy  + 
     &                       source(i1+1,j1  ,1)*(1.0-dx)*     dy  +
     &                       source(i1  ,j1+1,1)*     dx *(1.0-dy) +
     &                       source(i1+1,j1+1,1)*(1.0-dx)*(1.0-dy)
            enddo
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         do k=1, ddim3
            z = min(max(real(k-1)*fact3 + 0.5, half), edge3)
            k1 = int(z + 0.5)
            dz = real(k1) + 0.5 - z
            do j=1, ddim2
               y = min(max(real(j-1)*fact2 + 0.5, half), edge2)
               j1 = int(y + 0.5)
               dy = real(j1) + 0.5 - y
               do i=1, ddim1
                  x = min(max(real(i-1)*fact1 + 0.5, half), edge1)
                  i1 = int(x + 0.5)
                  dx = real(i1) + 0.5 - x
                  dest(i,j,k) = 
     &              source(i1  ,j1  ,k1  )*     dx *     dy *     dz  +
     &              source(i1+1,j1  ,k1  )*(1.0-dx)*     dy *     dz  +
     &              source(i1  ,j1+1,k1  )*     dx *(1.0-dy)*     dz  +
     &              source(i1+1,j1+1,k1  )*(1.0-dx)*(1.0-dy)*     dz  +
     &              source(i1  ,j1  ,k1+1)*     dx *     dy *(1.0-dz) +
     &              source(i1+1,j1  ,k1+1)*(1.0-dx)*     dy *(1.0-dz) +
     &              source(i1  ,j1+1,k1+1)*     dx *(1.0-dy)*(1.0-dz) +
     &              source(i1+1,j1+1,k1+1)*(1.0-dx)*(1.0-dy)*(1.0-dz)
               enddo
            enddo
         enddo
      endif
c
      return
      end
