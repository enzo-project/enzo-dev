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
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim
      R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3)
c
c  locals
c
      INTG_PREC i, j, k, i1, j1, k1
      R_PREC    fact1, fact2, fact3, x, y, z, dx, dy, dz, 
     &        edge1, edge2, edge3, half
      parameter (half = 0.5001_RKIND)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Precompute some things
c
      fact1 = REAL(sdim1-1,RKIND)/REAL(ddim1-1,RKIND)
      if (ndim .ge. 2) fact2 = REAL(sdim2-1,RKIND)/REAL(ddim2-1,RKIND)
      if (ndim .ge. 3) fact3 = REAL(sdim3-1,RKIND)/REAL(ddim3-1,RKIND)
      edge1 = REAL(sdim1,RKIND) - half
      edge2 = REAL(sdim2,RKIND) - half
      edge3 = REAL(sdim3,RKIND) - half
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=1, ddim1
            x = min(max(REAL(i-1,RKIND)*fact1 + 0.5_RKIND, half), edge1)
            i1 = int(x + 0.5_RKIND,IKIND)
            dx = REAL(i1,RKIND) + 0.5_RKIND - x
            dest(i,1,1) = source(i1,1,1)*dx + 
     &           source(i1+1,1,1)*(1._RKIND-dx)
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=1, ddim2
            y = min(max(REAL(j-1,RKIND)*fact2 + 0.5_RKIND, half), edge2)
            j1 = int(y + 0.5_RKIND,IKIND)
            dy = REAL(j1,RKIND) + 0.5_RKIND - y
            do i=1, ddim1
               x = min(max(REAL(i-1,RKIND)*fact1 + 0.5_RKIND, half), 
     &              edge1)
               i1 = int(x + 0.5_RKIND,IKIND)
               dx = REAL(i1,RKIND) + 0.5_RKIND - x
               dest(i,j,1) = source(i1  ,j1  ,1)*     dx *     dy  + 
     &              source(i1+1,j1  ,1)*(1._RKIND-dx)*     dy  +
     &              source(i1  ,j1+1,1)*     dx *(1._RKIND-dy) +
     &              source(i1+1,j1+1,1)*(1._RKIND-dx)*(1._RKIND-dy)
            enddo
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         do k=1, ddim3
            z = min(max(REAL(k-1,RKIND)*fact3 + 0.5_RKIND, half), edge3)
            k1 = int(z + 0.5_RKIND,IKIND)
            dz = REAL(k1,RKIND) + 0.5_RKIND - z
            do j=1, ddim2
               y = min(max(REAL(j-1,RKIND)*fact2 + 0.5_RKIND, half), 
     &              edge2)
               j1 = int(y + 0.5_RKIND,IKIND)
               dy = REAL(j1,RKIND) + 0.5_RKIND - y
               do i=1, ddim1
                  x = min(max(REAL(i-1,RKIND)*fact1 + 0.5_RKIND, half), 
     &                 edge1)
                  i1 = int(x + 0.5_RKIND,IKIND)
                  dx = REAL(i1,RKIND) + 0.5_RKIND - x
                  dest(i,j,k) = 
     &                 source(i1  ,j1  ,k1  )
     &                    *     dx *     dy *     dz  +
     &                 source(i1+1,j1  ,k1  )
     &                    *(1._RKIND-dx)*     dy *     dz  +
     &                 source(i1  ,j1+1,k1  )
     &                    *     dx *(1._RKIND-dy)*     dz  +
     &                 source(i1+1,j1+1,k1  )
     &                    *(1._RKIND-dx)*(1._RKIND-dy)*     dz +
     &                 source(i1  ,j1  ,k1+1)
     &                    *     dx *     dy *(1._RKIND-dz) +
     &                 source(i1+1,j1  ,k1+1)
     &                    *(1._RKIND-dx)*     dy *(1._RKIND-dz)+
     &                 source(i1  ,j1+1,k1+1)
     &                    *     dx *(1._RKIND-dy)*(1._RKIND-dz)+
     &                 source(i1+1,j1+1,k1+1)
     &                    *(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)
               enddo
            enddo
         enddo
      endif
c
      return
      end
