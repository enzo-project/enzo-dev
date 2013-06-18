#include "fortran.def"
c=======================================================================
c///////////////////////  SUBROUTINE DEP_GRID_NGP  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine dep_grid_ngp(source, dest, temp, velx, vely, velz, 
     &                    dt, rfield, ndim, ihydro, delx, dely, delz,
     &                    sdim1, sdim2, sdim3, 
     &                    sstart1, sstart2, sstart3, 
     &                    send1, send2, send3,
     &                    offset1, offset2, offset3,
     &                    ddim1, ddim2, ddim3,
     &                    refine1, refine2, refine3)
c
c  DEPOSIT SOURCE GRID INTO DEST GRID USING NGP INTERPOLATION
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
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC ddim1, ddim2, ddim3, sdim1, sdim2, sdim3, ndim, ihydro,
     &        refine1, refine2, refine3, sstart1, sstart2, sstart3,
     &        send1, send2, send3
      R_PREC    source(sdim1, sdim2, sdim3), dest(ddim1, ddim2, ddim3),
     &        rfield(sdim1, sdim2, sdim3),
     &        velx(sdim1, sdim2, sdim3), vely(sdim1, sdim2, sdim3),
     &        velz(sdim1, sdim2, sdim3), dt, delx, dely, delz,
     &        offset1, offset2, offset3,
     &        temp(ddim1, ddim2, ddim3, 4)
c
c  locals
c
      INTG_PREC i, j, k, i1, j1, k1, n
      R_PREC    fact1, fact2, fact3, weight, mass,
     &          start1, start2, start3
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
!      write(6,*) 'sdim:',sdim1,sdim2,sdim3
!      write(6,*) 'ddim:',ddim1,ddim2,ddim3
!      write(6,*) 'offset:',offset1,offset2,offset3
!      write(6,*) 'sstart:',sstart1,sstart2,sstart3
c
c     Clear dest and temp_vel fields.
c
      do k=1,ddim3
         do j=1,ddim2
            do i=1,ddim1
               dest(i,j,k) = 0._RKIND
            enddo
         enddo
      enddo
c
c     Precompute some things
c
      fact1 = 1._RKIND/REAL(refine1,RKIND)
      fact2 = 1._RKIND/REAL(refine2,RKIND)
      fact3 = 1._RKIND/REAL(refine3,RKIND)
c     
      start1 = -sstart1 + offset1*REAL(refine1,RKIND)
      start2 = -sstart2 + offset2*REAL(refine2,RKIND)
      start3 = -sstart3 + offset3*REAL(refine3,RKIND)
c
c     a) 1D
c
      if (ndim .eq. 1) then
         weight = fact1
c
c        compute density and mass-weighted velocity field
c
         do i=sstart1+1, send1+1
            i1 = min(max(int((start1 + i - 0.5_RKIND)*fact1) + 1, 1), 
     &           ddim1)
c
            mass = (1._RKIND - rfield(i,1,1))*weight*source(i,1,1)
c
            dest(i1, 1, 1) = dest(i1, 1, 1) + mass
c
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
            j1 = min(max(int((start2 + j - 0.5_RKIND)*fact2) + 1, 1), 
     &              ddim2)
            do i=sstart1+1, send1+1
               i1 = min(max(int((start1 + i - 0.5_RKIND)*fact1) + 1, 1),
     &                 ddim1)
c
               mass = (1._RKIND - rfield(i,j,1))*weight*source(i,j,1)
c
               dest(i1  ,j1  ,1) = dest(i1  ,j1  ,1) + mass
c
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
            k1 = min(max(int((start3 + k - 0.5_RKIND)*fact3) + 1, 
     &           1), ddim3)
            do j=sstart2+1, send2+1
               j1 = min(max(int((start2 + j - 0.5_RKIND)*fact2) + 1, 
     &              1), ddim2)
               do i=sstart1+1, send1+1
                  i1 = min(max(int((start1 + i - 0.5_RKIND)*fact1) + 1, 
     &                 1), ddim1)
c
                  mass = (1._RKIND - rfield(i,j,k))*weight*source(i,j,k)
c
                  dest(i1  ,j1  ,k1  ) = dest(i1  ,j1  ,k1  ) + mass
c
               enddo
            enddo
         enddo
c
c         write(6,*) ddim1, ddim2, ddim3, refine1
c         do i=1, ddim1
c            write(6,*) i, dest(i,ddim2/2, ddim3/3)
c         enddo
c
      endif
c
      return
      end

