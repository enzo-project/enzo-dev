#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE DEP_GRID_CIC  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine dep_grid_cic(source, dest, temp, velx, vely, velz, 
     &                    dt, rfield, ndim, ihydro, delx, dely, delz,
     &                    sdim1, sdim2, sdim3, 
     &                    sstart1, sstart2, sstart3, 
     &                    send1, send2, send3,
     &                    offset1, offset2, offset3,
     &                    ddim1, ddim2, ddim3,
     &                    refine1, refine2, refine3)
c
c  DEPOSIT SOURCE GRID INTO DEST GRID USING CIC INTERPOLATION
c
c  written by: Greg Bryan
c  date:       March, 1999
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source       - source field
c     rfield       - source-like field indicating if cell is refined 
c                       (1=no, 0=yes)
c     sdim1-3      - source dimension
c     ddim1-3      - destination dimension
c     ndim         - rank of fields
c     refine1-3    - refinement factors
c     sstart1-3    - source start index
c     send1-3      - source end index
c     offset1-3     - offset from this grid edge to dest grid edge
c                    (>= 0, in dest cell units)
c     velx,y,z     - velocities
c     dt           - time step
c     delx         - cell size of source grid
c     temp         - temporary field, 4*size of dest
c     ihydro       - hydro method (2 - zeus, velocity is cell centered)
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
      integer ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim, ihydro,
     &        refine1, refine2, refine3, sstart1, sstart2, sstart3,
     &        send1, send2, send3
      real    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3),
     &        rfield(sdim1, sdim2, sdim3),
     &        velx(sdim1, sdim2, sdim3), vely(sdim1, sdim2, sdim3),
     &        velz(sdim1, sdim2, sdim3), dt, delx, dely, delz,
     &        offset1, offset2, offset3,
     &        temp(ddim1, ddim2, ddim3, 4)
c
c  locals
c
      integer i, j, k, i1, j1, k1, n
      real    fact1, fact2, fact3, x, y, z, dx, dy, dz, weight, mass,
     &        coef1, coef2, coef3, shift1, shift2, shift3,
     &        start1, start2, start3, half, edge1, edge2, edge3, temp1
      parameter (half = 0.5001)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Clear dest and temp_vel fields.
c

!     write(0,'("grid_cic: ",9i4)') sdim1,sdim2,sdim3,ddim1,ddim2,ddim3,
!    &                              offset1,offset2,offset3

      do k=1,ddim3
         do j=1,ddim2
            do i=1,ddim1
               dest(i,j,k) = 0
            enddo
            do n=1,ndim+1
               do i=1,ddim1
                  temp(i,j,k,n) = 0
               enddo
            enddo
         enddo
      enddo
c
c     Precompute some things
c
      fact1 = 1.0/real(refine1)
      fact2 = 1.0/real(refine2)
      fact3 = 1.0/real(refine3)
c
                       coef1 = dt/delx*fact1
      if (ndim .gt. 1) coef2 = dt/dely*fact2
      if (ndim .gt. 2) coef3 = dt/delz*fact3
c     
      start1 = -sstart1 - 0.5 + offset1*real(refine1)
      start2 = -sstart2 - 0.5 + offset2*real(refine2)
      start3 = -sstart3 - 0.5 + offset3*real(refine3)
c
      edge1 = real(ddim1) - half
      edge2 = real(ddim2) - half
      edge3 = real(ddim3) - half

!     write(0,'("grid_cic: ",6f12.4)') start1,start2,start3,
!    &                                 edge1,edge2,edge3

c
c     a) 1D
c
      if (ndim .eq. 1) then
         weight = fact1
c
c        compute density and mass-weighted velocity field
c
         do i=sstart1+1, send1+1
            x = min(max((start1 + i)*fact1, half), edge1)
            i1 = int(x + 0.5)
            dx = real(i1) + 0.5 - x
            mass = (1.0 - rfield(i,1,1))*weight*source(i,1,1)
c
            temp(i1  ,1  ,1, 1) = temp(i1  ,1  ,1, 1) + 
     &                                 mass*     dx
            temp(i1+1,1  ,1, 1) = temp(i1+1,1  ,1, 1) + 
     &                                 mass*(1.0-dx)
c
            temp1 = velx(i,1,1)*mass
            if (ihydro .eq. 2) temp1 = 
     &           0.5*(velx(i,1,1)+velx(i+1,1,1))*mass
            temp(i1  ,1  ,1, 2) = temp(i1  ,1  ,1, 2) + 
     &                     temp1*     dx
            temp(i1+1,1  ,1, 2) = temp(i1+1,1  ,1, 2) + 
     &                     temp1*(1.0-dx)
         enddo
c
c        Use velocity and mass field to generate mass field advanced by dt
c
         do i=1,ddim1
            shift1 = temp(i,1,1,2)/max(temp(i,1,1,1),tiny)*coef1
            x = min(max((i - 0.5 + shift1), half), edge1)
            i1 = int(x + 0.5)
            dx = real(i1) + 0.5 - x
            mass = temp(i,1,1,1)
            dest(i1  ,1  ,1) = dest(i1  ,1  ,1) + mass*     dx
            dest(i1+1,1  ,1) = dest(i1+1,1  ,1) + mass*(1.0-dx)
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         weight = fact1*fact2
c
c        compute density and mass-weighted velocity field
c
         do j=sstart2+1, send2+1
            y = min(max((start2 + j)*fact2, half), edge2)
            j1 = int(y + 0.5)
            dy = real(j1) + 0.5 - y
            do i=sstart1+1, send1+1
               x = min(max((start1 + i)*fact1, half), edge1)
               i1 = int(x + 0.5)
               dx = real(i1) + 0.5 - x
               mass = (1.0 - rfield(i,j,1))*weight*source(i,j,1)
c
               temp(i1  ,j1  ,1,1) = temp(i1  ,j1  ,1,1) +
     &                             mass*     dx *     dy
               temp(i1+1,j1  ,1,1) = temp(i1+1,j1  ,1,1) +
     &                             mass*(1.0-dx)*     dy
               temp(i1  ,j1+1,1,1) = temp(i1  ,j1+1,1,1) +
     &                             mass*     dx *(1.0-dy)
               temp(i1+1,j1+1,1,1) = temp(i1+1,j1+1,1,1) +
     &                             mass*(1.0-dx)*(1.0-dy)
c
               temp1 = velx(i,j,1)*mass
               if (ihydro .eq. 2) temp1 = 
     &              0.5*(velx(i,j,1)+velx(i+1,j,1))*mass
               temp(i1  ,j1  ,1,2) = temp(i1  ,j1  ,1,2) +
     &                 temp1*     dx *     dy
               temp(i1+1,j1  ,1,2) = temp(i1+1,j1  ,1,2) +
     &                 temp1*(1.0-dx)*     dy
               temp(i1  ,j1+1,1,2) = temp(i1  ,j1+1,1,2) +
     &                 temp1*     dx *(1.0-dy)
               temp(i1+1,j1+1,1,2) = temp(i1+1,j1+1,1,2) +
     &                 temp1*(1.0-dx)*(1.0-dy)
c
               temp1 = vely(i,j,1)*mass
               if (ihydro .eq. 2) temp1 = 
     &              0.5*(vely(i,j,1)+vely(i,j+1,1))*mass
               temp(i1  ,j1  ,1,3) = temp(i1  ,j1  ,1,3) +
     &                 temp1*     dx *     dy
               temp(i1+1,j1  ,1,3) = temp(i1+1,j1  ,1,3) +
     &                 temp1*(1.0-dx)*     dy
               temp(i1  ,j1+1,1,3) = temp(i1  ,j1+1,1,3) +
     &                 temp1*     dx *(1.0-dy)
               temp(i1+1,j1+1,1,3) = temp(i1+1,j1+1,1,3) +
     &                 temp1*(1.0-dx)*(1.0-dy)
c
            enddo
         enddo
c
c        Use velocity and mass field to generate mass field advanced by dt
c
         do j=1, ddim2
            do i=1, ddim1
               shift2 = temp(i,j,1,3)/
     &                  max(temp(i,j,1,1),tiny)*coef2
               y = min(max((j - 0.5 + shift2), half), edge2)
               j1 = int(y + 0.5)
               dy = real(j1) + 0.5 - y
               shift1 = temp(i,j,1,2)/
     &                  max(temp(i,j,1,1),tiny)*coef1
               x = min(max((i - 0.5 + shift1), half), edge1)
               i1 = int(x + 0.5)
               dx = real(i1) + 0.5 - x
               mass = temp(i,j,1,1)
               dest(i1  ,j1  ,1) = dest(i1  ,j1  ,1) +
     &                             mass*     dx *     dy
               dest(i1+1,j1  ,1) = dest(i1+1,j1  ,1) +
     &                             mass*(1.0-dx)*     dy
               dest(i1  ,j1+1,1) = dest(i1  ,j1+1,1) +
     &                             mass*     dx *(1.0-dy)
               dest(i1+1,j1+1,1) = dest(i1+1,j1+1,1) +
     &                             mass*(1.0-dx)*(1.0-dy)
            enddo
         enddo
      endif
c
c     c) 3D
c
      if (ndim .eq. 3) then
         weight = fact1*fact2*fact3
c
c        compute density and mass-weighted velocity field
c
         do k=sstart3+1, send3+1
            z = min(max((start3 + k)*fact3, half), edge3)
            k1 = int(z + 0.5)
            dz = real(k1) + 0.5 - z
            do j=sstart2+1, send2+1
               y = min(max((start2 + j)*fact2, half), edge2)
               j1 = int(y + 0.5)
               dy = real(j1) + 0.5 - y
               do i=sstart1+1, send1+1
                  x = min(max((start1 + i)*fact1, half), edge1)
                  i1 = int(x + 0.5)
                  dx = real(i1) + 0.5 - x
                  mass = (1.0 - rfield(i,j,k))*weight*source(i,j,k)
c
                  temp(i1  ,j1  ,k1  ,1) = temp(i1  ,j1  ,k1  ,1) +
     &                       mass*     dx *     dy *     dz
                  temp(i1+1,j1  ,k1  ,1) = temp(i1+1,j1  ,k1  ,1) +
     &                       mass*(1.0-dx)*     dy *     dz
                  temp(i1  ,j1+1,k1  ,1) = temp(i1  ,j1+1,k1  ,1) +
     &                       mass*     dx *(1.0-dy)*     dz
                  temp(i1+1,j1+1,k1  ,1) = temp(i1+1,j1+1,k1  ,1) +
     &                       mass*(1.0-dx)*(1.0-dy)*     dz
                  temp(i1  ,j1  ,k1+1,1) = temp(i1  ,j1  ,k1+1,1) +
     &                       mass*     dx *     dy *(1.0-dz)
                  temp(i1+1,j1  ,k1+1,1) = temp(i1+1,j1  ,k1+1,1) +
     &                       mass*(1.0-dx)*     dy *(1.0-dz)
                  temp(i1  ,j1+1,k1+1,1) = temp(i1  ,j1+1,k1+1,1) +
     &                       mass*     dx *(1.0-dy)*(1.0-dz)
                  temp(i1+1,j1+1,k1+1,1) = temp(i1+1,j1+1,k1+1,1) +
     &                       mass*(1.0-dx)*(1.0-dy)*(1.0-dz)
c
                  temp1 = velx(i,j,k)*mass
                  if (ihydro .eq. 2) temp1 = 
     &                 0.5*(velx(i,j,k)+velx(i+1,j,k))*mass
                  temp(i1  ,j1  ,k1  ,2) = temp(i1  ,j1  ,k1  ,2) +
     &                       temp1*     dx *     dy *     dz
                  temp(i1+1,j1  ,k1  ,2) = temp(i1+1,j1  ,k1  ,2) +
     &                       temp1*(1.0-dx)*     dy *     dz
                  temp(i1  ,j1+1,k1  ,2) = temp(i1  ,j1+1,k1  ,2) +
     &                       temp1*     dx *(1.0-dy)*     dz
                  temp(i1+1,j1+1,k1  ,2) = temp(i1+1,j1+1,k1  ,2) +
     &                       temp1*(1.0-dx)*(1.0-dy)*     dz
                  temp(i1  ,j1  ,k1+1,2) = temp(i1  ,j1  ,k1+1,2) +
     &                       temp1*     dx *     dy *(1.0-dz)
                  temp(i1+1,j1  ,k1+1,2) = temp(i1+1,j1  ,k1+1,2) +
     &                       temp1*(1.0-dx)*     dy *(1.0-dz)
                  temp(i1  ,j1+1,k1+1,2) = temp(i1  ,j1+1,k1+1,2) +
     &                       temp1*     dx *(1.0-dy)*(1.0-dz)
                  temp(i1+1,j1+1,k1+1,2) = temp(i1+1,j1+1,k1+1,2) +
     &                       temp1*(1.0-dx)*(1.0-dy)*(1.0-dz)
c
                  temp1 = vely(i,j,k)*mass
                  if (ihydro .eq. 2) temp1 = 
     &                 0.5*(vely(i,j,k)+vely(i,j+1,k))*mass
                  temp(i1  ,j1  ,k1  ,3) = temp(i1  ,j1  ,k1  ,3) +
     &                       temp1*     dx *     dy *     dz
                  temp(i1+1,j1  ,k1  ,3) = temp(i1+1,j1  ,k1  ,3) +
     &                       temp1*(1.0-dx)*     dy *     dz
                  temp(i1  ,j1+1,k1  ,3) = temp(i1  ,j1+1,k1  ,3) +
     &                       temp1*     dx *(1.0-dy)*     dz
                  temp(i1+1,j1+1,k1  ,3) = temp(i1+1,j1+1,k1  ,3) +
     &                       temp1*(1.0-dx)*(1.0-dy)*     dz
                  temp(i1  ,j1  ,k1+1,3) = temp(i1  ,j1  ,k1+1,3) +
     &                       temp1*     dx *     dy *(1.0-dz)
                  temp(i1+1,j1  ,k1+1,3) = temp(i1+1,j1  ,k1+1,3) +
     &                       temp1*(1.0-dx)*     dy *(1.0-dz)
                  temp(i1  ,j1+1,k1+1,3) = temp(i1  ,j1+1,k1+1,3) +
     &                       temp1*     dx *(1.0-dy)*(1.0-dz)
                  temp(i1+1,j1+1,k1+1,3) = temp(i1+1,j1+1,k1+1,3) +
     &                       temp1*(1.0-dx)*(1.0-dy)*(1.0-dz)
c
                  temp1 = velz(i,j,k)*mass
                  if (ihydro .eq. 2) temp1 = 
     &                 0.5*(velz(i,j,k)+velz(i,j,k+1))*mass
                  temp(i1  ,j1  ,k1  ,4) = temp(i1  ,j1  ,k1  ,4) +
     &                       temp1*     dx *     dy *     dz
                  temp(i1+1,j1  ,k1  ,4) = temp(i1+1,j1  ,k1  ,4) +
     &                       temp1*(1.0-dx)*     dy *     dz
                  temp(i1  ,j1+1,k1  ,4) = temp(i1  ,j1+1,k1  ,4) +
     &                       temp1*     dx *(1.0-dy)*     dz
                  temp(i1+1,j1+1,k1  ,4) = temp(i1+1,j1+1,k1  ,4) +
     &                       temp1*(1.0-dx)*(1.0-dy)*     dz
                  temp(i1  ,j1  ,k1+1,4) = temp(i1  ,j1  ,k1+1,4) +
     &                       temp1*     dx *     dy *(1.0-dz)
                  temp(i1+1,j1  ,k1+1,4) = temp(i1+1,j1  ,k1+1,4) +
     &                       temp1*(1.0-dx)*     dy *(1.0-dz)
                  temp(i1  ,j1+1,k1+1,4) = temp(i1  ,j1+1,k1+1,4) +
     &                       temp1*     dx *(1.0-dy)*(1.0-dz)
                  temp(i1+1,j1+1,k1+1,4) = temp(i1+1,j1+1,k1+1,4) +
     &                       temp1*(1.0-dx)*(1.0-dy)*(1.0-dz)
c
               enddo
            enddo
         enddo
c
c        Use velocity and mass field to generate mass field advanced by dt
c
         do k=1, ddim3
            do j=1, ddim2
               do i=1, ddim1
                  shift3 = temp(i,j,k,4)/
     &                     max(temp(i,j,k,1),tiny)*coef3
                  z = min(max((k - 0.5 + shift3), half), edge3)
                  k1 = int(z + 0.5)
                  dz = real(k1) + 0.5 - z
c
                  shift2 = temp(i,j,k,3)/
     &                     max(temp(i,j,k,1),tiny)*coef2
                  y = min(max((j - 0.5 + shift2), half), edge2)
                  j1 = int(y + 0.5)
                  dy = real(j1) + 0.5 - y
c
                  shift1 = temp(i,j,k,2)/
     &                     max(temp(i,j,k,1),tiny)*coef1
                  x = min(max((i - 0.5 + shift1), half), edge1)
                  i1 = int(x + 0.5)
                  dx = real(i1) + 0.5 - x
c
c                  if (i1 .lt. 1 .or. i1 .ge. ddim1 .or.
c     &                j1 .lt. 1 .or. j1 .ge. ddim2 .or.
c     &                k1 .lt. 1 .or. k1 .ge. ddim3    )
c     &                write(6,*) i1,j1,k1,ddim1,ddim2,ddim3
c
                  mass = temp(i,j,k,1)
                  dest(i1  ,j1  ,k1  ) = dest(i1  ,j1  ,k1  ) +
     &                       mass*     dx *     dy *     dz
                  dest(i1+1,j1  ,k1  ) = dest(i1+1,j1  ,k1  ) +
     &                       mass*(1.0-dx)*     dy *     dz
                  dest(i1  ,j1+1,k1  ) = dest(i1  ,j1+1,k1  ) +
     &                       mass*     dx *(1.0-dy)*     dz
                  dest(i1+1,j1+1,k1  ) = dest(i1+1,j1+1,k1  ) +
     &                       mass*(1.0-dx)*(1.0-dy)*     dz
                  dest(i1  ,j1  ,k1+1) = dest(i1  ,j1  ,k1+1) +
     &                       mass*     dx *     dy *(1.0-dz)
                  dest(i1+1,j1  ,k1+1) = dest(i1+1,j1  ,k1+1) +
     &                       mass*(1.0-dx)*     dy *(1.0-dz)
                  dest(i1  ,j1+1,k1+1) = dest(i1  ,j1+1,k1+1) +
     &                       mass*     dx *(1.0-dy)*(1.0-dz)
                  dest(i1+1,j1+1,k1+1) = dest(i1+1,j1+1,k1+1) +
     &                       mass*(1.0-dx)*(1.0-dy)*(1.0-dz)
               enddo
            enddo
         enddo
c
      endif
c
      return
      end



c=======================================================================
c///////////////////////  SUBROUTINE INT_GRID_CIC  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine int_grid_cic(source, dest, ndim,
     &                        sdim1, sdim2, sdim3,
     &                        ddim1, ddim2, ddim3,
     &                        start1, start2, start3,
     &                        refine1, refine2, refine3)
c
c  INTERPOLATE FROM SOURCE GRID INTO DEST
c
c  written by: Greg Bryan
c  date:       March, 1999
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     source       - source field
c     sdim1-3      - source dimension
c     ddim1-3      - destination dimension
c     ndim         - rank of fields
c     start1-3     - start of dest field in dest grid cells (real)
c     refine1-3    - refinement factors
c     sstart1-3    - source start index
c     send1-3      - source end index
c
c  OUTPUT ARGUMENTS:
c     dest         - interpolated field
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
     &        refine1, refine2, refine3
      real    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3),
     &        start1, start2, start3
c
c  locals
c
      integer i, j, k, i1, j1, k1
      real    fact1, fact2, fact3, x, y, z, dx, dy, dz
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Precompute some things
c
      fact1 = 1.0/real(refine1)
      fact2 = 1.0/real(refine2)
      fact3 = 1.0/real(refine3)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         do i=1, ddim1
            x = (start1 + i - 0.5)*fact1
            i1 = int(x + 0.5)
            dx = real(i1) + 0.5 - x
            dest(i,1,1) = source(i1  ,1  ,1)*     dx +
     &                    source(i1+1,1  ,1)*(1.0-dx)
         enddo
      endif
c
c     b) 2D
c
      if (ndim .eq. 2) then
         do j=1, ddim2
            y = (start2 + j - 0.5)*fact2
            j1 = int(y + 0.5)
            dy = real(j1) + 0.5 - y
            do i=1, ddim1
               x = (start1 + i - 0.5)*fact1
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
            z = (start3 + k - 0.5)*fact3
            k1 = int(z + 0.5)
            dz = real(k1) + 0.5 - z
            do j=1, ddim2
               y = (start2 + j - 0.5)*fact2
               j1 = int(y + 0.5)
               dy = real(j1) + 0.5 - y
               do i=1, ddim1
                  x = (start1 + i - 0.5)*fact1
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
