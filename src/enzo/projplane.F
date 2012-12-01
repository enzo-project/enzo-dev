#include "fortran.def"
c=======================================================================
c////////////////////////  SUBROUTINE PROJECT  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine projplane(grid1, grid2, flaggrid, iflag, ismooth,
     &                     gdim1, gdim2, gdim3, gcellsize,
     &                     plane, pdim1, pdim2, pcellsize,
     &                     projdim, ifield, weight,
     &                     gleft, gfarleft, gright, pleft, pright,
     &                     npstart, npend, fracleft, fracright)
c
c  PROJECTS A 3D FLOAT FIELD TO A 2D PLANE
c
c  written by: Greg Bryan
c  date:       April, 1996
c  modified1:
c
c  PURPOSE:
c
c  INPUTS:
c     grid1      - 3D grid to be projected
c     grid2      - 3D grid to be used if ifield .ne. 1
c     flaggrid   - 3D field of real flags (0/1 - do/do not assign to this pt)
c     iflag      - flag indicating if flaggrid should be used (1 - yes)
c     ismooth    - flag indicating if smoothing should be used (1 - yes)
c     gdim1,2,3  - declared dimensions of grid1,2
c     gcellsize  - cell size of grid
c     plane      - 2D projection plane
c     pdim1,2    - declared dimensions if plane
c     pcellsize  - cell size of plane
c     projdim    - projection direction (0 - x, 1 - y, 2 - z)
c     ifield     - flag to specify projection calc method
c                     1 - multiply field1 by weight
c                     2 - multiply sqrt(field1)*field2**2 by weight
c     weight     - weighting value (i.e. cell volume) for plane
c     gleft      - start of active region in grid
c     gfarleft   - start of grid
c     gright     - end of active region in grid
c     pleft      - start of projection plane (in 3D)
c     pright     - end of projection plane (in 3D)
c     npstart    - start index in projection direction for grid (0 based)
c     npend      - end index in projection direction for grid (0 based)
c     fracleft   - the fraction of the leftmost plane to be used (0-1)
c     fracright  - the fraction of the rightmost plane to be used (0-1)
c
c  OUTPUT ARGUMENTS: 
c     line     - projected line
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
      INTG_PREC gdim1, gdim2, gdim3, ifield, iflag, ismooth, 
     &        pdim1, pdim2, projdim, npstart, npend
      R_PREC    weight, fracleft, fracright
      R_PREC    grid1(gdim1,gdim2,gdim3), grid2(gdim1,gdim2,gdim3),
     &        flaggrid(gdim1,gdim2,gdim3), plane(pdim1,pdim2)
      P_PREC gleft(3), gfarleft(3), gright(3), pleft(3), pright(3),
     &        pcellsize, gcellsize
c
c  locals
c
      INTG_PREC i, i0, j, j0, k0, n, istart(3), iend(3)
      R_PREC    sample1, sample2, dx, dy, dz, weight1
      R_PREC  xpos, ypos, zpos
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
c     Compute start and end indexes (in plane space) for the projected region
c
c      write(6,*) 'projplane: projdim,npstrt,end =',projdim,npstart,npend
c      write(6,*) 'projplane: p,gcellsize,weight =',pcellsize,gcellsize,
c     &                                             weight
      do n = 1, 3
         istart(n) = nint((max(gleft(n), pleft(n))   - pleft(n))/
     &               pcellsize,IKIND)
         iend(n)   = nint((min(gright(n), pright(n)) - pleft(n))/
     &               pcellsize,IKIND)
c         write(6,*) 'projplane: i,istart,iend =',n,istart(n),iend(n)
c         write(6,*) 'projplane: gleft, gright =',gleft(n),gright(n)
c         write(6,*) 'projplane: pleft, pright =',pleft(n),pright(n)
      enddo
c
c     Loop over the grid indexes along the projection axis
c
      do n = npstart+1, npend+1
c
      weight1 = weight
      if (n .eq. npstart+1) weight1 = weight * fracleft
      if (n .eq. npend  +1) weight1 = weight * fracright
c
       if (ismooth .eq. 1) then
c
c=======================================================================
c      linear interpolation
c
c ----------------------------------------------------------------------
c        project in i-direction
c
         if (projdim .eq. 0) then
c
c           Loop over the projected region of the plane
c
            do j = istart(3)+1, iend(3)
               do i = istart(2)+1, iend(2)
c
c                 Compute position in grid space
c
                  ypos = pleft(2) + (REAL(i,RKIND)-0.5_RKIND)*pcellsize
                  zpos = pleft(3) + (REAL(j,RKIND)-0.5_RKIND)*pcellsize
c
c                 Compute indexes in grid (for base cell)
c
                  j0 = int((ypos-gfarleft(2))/gcellsize+0.5_RKIND,IKIND)
                  k0 = int((zpos-gfarleft(3))/gcellsize+0.5_RKIND,IKIND)
c
c                 Compute weights for interpolation
c
                  dy = REAL(j0,RKIND) - (ypos - gfarleft(2))/gcellsize
     &                 + 0.5_RKIND
                  dz = REAL(k0,RKIND) - (zpos - gfarleft(3))/gcellsize
     &                 + 0.5_RKIND
c
c                 Do 2D linear interpolation into grid field
c
                  sample1 = grid1(n,j0  ,k0  )*     dy *     dz  +
     &                 grid1(n,j0  ,k0+1)*     dy *(1._RKIND-dz) +
     &                 grid1(n,j0+1,k0  )*(1._RKIND-dy)*     dz  +
     &                 grid1(n,j0+1,k0+1)*(1._RKIND-dy)*(1._RKIND-dz)
c
                  if (iflag .eq. 1) then
                     if (flaggrid(n,
     &                    int((ypos-gfarleft(2))/gcellsize,IKIND)+1,
     &                    int((zpos-gfarleft(3))/gcellsize,IKIND)+1) 
     &                   .ne. 0.0) sample1 = 0._RKIND
                  endif
c
                  if (ifield .eq. 1) then
                     plane(i,j) = plane(i,j) + sample1 * weight1
                  else
c
                     sample2 = grid2(n,j0  ,k0  )*     dy *     dz  +
     &                    grid2(n,j0  ,k0+1)*     dy *(1._RKIND-dz) +
     &                    grid2(n,j0+1,k0  )*(1._RKIND-dy)*     dz  +
     &                    grid2(n,j0+1,k0+1)*(1._RKIND-dy)*(1._RKIND-dz)
c
                     if (ifield .eq. 2) plane(i,j) = plane(i,j) + 
     &                    sqrt(sample1)*sample2**2 * weight1
                     if (ifield .eq. 3) plane(i,j) = plane(i,j) + 
     &                    sample1**1.5_RKIND*sample2**2 *weight1
c
                  endif
c     
               enddo
            enddo
c     
         endif
c
c ----------------------------------------------------------------------
c   project in j-direction
c
         if (projdim .eq. 1) then
c
            do j = istart(3)+1, iend(3)
               do i = istart(1)+1, iend(1)
c
                  xpos = pleft(1) + (REAL(i,RKIND)-0.5_RKIND)*pcellsize
                  zpos = pleft(3) + (REAL(j,RKIND)-0.5_RKIND)*pcellsize
c
                  i0 = int((xpos-gfarleft(1))/gcellsize+0.5_RKIND,IKIND)
                  k0 = int((zpos-gfarleft(3))/gcellsize+0.5_RKIND,IKIND)
c
                  dx = REAL(i0,RKIND) - (xpos - gfarleft(1))/gcellsize 
     &                 + 0.5_RKIND
                  dz = REAL(k0,RKIND) - (zpos - gfarleft(3))/gcellsize 
     &                 + 0.5_RKIND
c
                  sample1 = grid1(i0  ,n,k0  )*     dx *     dz  +
     &                 grid1(i0  ,n,k0+1)*     dx *(1._RKIND-dz) +
     &                 grid1(i0+1,n,k0  )*(1._RKIND-dx)*     dz  +
     &                 grid1(i0+1,n,k0+1)*(1._RKIND-dx)*(1._RKIND-dz)
c
                  if (iflag .eq. 1) then
                     if (flaggrid(
     &                    int((xpos-gfarleft(1))/gcellsize,IKIND)+1,n,
     &                    int((zpos-gfarleft(3))/gcellsize,IKIND)+1) 
     &                   .ne. 0._RKIND) sample1 = 0._RKIND
                  endif
c
                  if (ifield .eq. 1) then
                     plane(i,j) = plane(i,j) + sample1 * weight1
                  else
c
                     sample2 = grid2(i0  ,n,k0  )*     dx *     dz  +
     &                    grid2(i0  ,n,k0+1)*     dx *(1._RKIND-dz) +
     &                    grid2(i0+1,n,k0  )*(1._RKIND-dx)*     dz  +
     &                    grid2(i0+1,n,k0+1)*(1._RKIND-dx)*(1._RKIND-dz)
c
                     if (ifield .eq. 2) plane(i,j) = plane(i,j) + 
     &                    sqrt(sample1)*sample2**2 * weight1
                     if (ifield .eq. 3) plane(i,j) = plane(i,j) + 
     &                    sample1**1.5_RKIND*sample2**2 *weight1
c
                  endif
c     
               enddo
            enddo
c     
         endif
c
c ----------------------------------------------------------------------
c   project in k-direction
c
         if (projdim .eq. 2) then
c
            do j = istart(2)+1, iend(2)
               do i = istart(1)+1, iend(1)
c
                  xpos = pleft(1) + (REAL(i,RKIND)-0.5)*pcellsize
                  ypos = pleft(2) + (REAL(j,RKIND)-0.5)*pcellsize
c
                  i0 = int((xpos-gfarleft(1))/gcellsize+0.5_RKIND,IKIND)
                  j0 = int((ypos-gfarleft(2))/gcellsize+0.5_RKIND,IKIND)
c
                  dx = REAL(i0,RKIND) - (xpos - gfarleft(1))/gcellsize 
     &                 + 0.5_RKIND
                  dy = REAL(j0,RKIND) - (ypos - gfarleft(2))/gcellsize 
     &                 + 0.5_RKIND
c
                  sample1 = grid1(i0  ,j0  ,n)*     dx *     dy  +
     &                 grid1(i0  ,j0+1,n)*     dx *(1._RKIND-dy) +
     &                 grid1(i0+1,j0  ,n)*(1._RKIND-dx)*     dy  +
     &                 grid1(i0+1,j0+1,n)*(1._RKIND-dx)*(1._RKIND-dy)
c     
                  if (iflag .eq. 1) then
                     if (flaggrid(
     &                    int((xpos-gfarleft(1))/gcellsize,IKIND)+1,
     &                    int((ypos-gfarleft(2))/gcellsize,IKIND)+1,n) 
     &                   .ne. 0._RKIND) sample1 = 0._RKIND
                  endif
c
                  if (ifield .eq. 1) then
                     plane(i,j) = plane(i,j) + sample1 * weight1
                  else
c
                     sample2 = grid2(i0  ,j0  ,n)*     dx *     dy  +
     &                    grid2(i0  ,j0+1,n)*     dx *(1._RKIND-dy) +
     &                    grid2(i0+1,j0  ,n)*(1._RKIND-dx)*     dy  +
     &                    grid2(i0+1,j0+1,n)*(1._RKIND-dx)*(1._RKIND-dy)
c     
                     if (ifield .eq. 2) plane(i,j) = plane(i,j) + 
     &                    sqrt(sample1)*sample2**2 * weight1
                     if (ifield .eq. 3) plane(i,j) = plane(i,j) + 
     &                    sample1**1.5_RKIND*sample2**2 *weight1
c
                  endif
c     
               enddo
            enddo
c     
         endif
c
        else
c
c=======================================================================
c       Zero-level interpolation
c
c ----------------------------------------------------------------------
c        project in i-direction
c
         if (projdim .eq. 0) then
c
c           Loop over the projected region of the plane
c
            do j = istart(3)+1, iend(3)
               do i = istart(2)+1, iend(2)
c
                  j0 = int((pleft(2) + (REAL(i,RKIND)-0.5_RKIND)
     &                 * pcellsize - gfarleft(2))/gcellsize,IKIND) + 1
                  k0 = int((pleft(3) + (REAL(j,RKIND)-0.5_RKIND)
     &                 * pcellsize - gfarleft(3))/gcellsize,IKIND) + 1
c
                  sample1 = grid1(n,j0,k0)
c
                  if (iflag .eq. 1) then
                     if (flaggrid(n,j0,k0) .ne. 0.0) sample1 = 0._RKIND
                  endif
c
                  if (ifield .eq. 1) plane(i,j) = plane(i,j) + 
     &                 sample1 * weight1
                  if (ifield .eq. 2) plane(i,j) = plane(i,j) + 
     &                 sqrt(sample1) * grid2(n,j0,k0)**2 * weight1
                  if (ifield .eq. 3) plane(i,j) = plane(i,j) + 
     &                 sample1**1.5_RKIND *grid2(n,j0,k0)**2 * weight1
c     
               enddo
            enddo
c     
         endif
c
c ----------------------------------------------------------------------
c        project in j-direction
c
         if (projdim .eq. 1) then
c
c           Loop over the projected region of the plane
c
            do j = istart(3)+1, iend(3)
               do i = istart(1)+1, iend(1)
c
                  i0 = int((pleft(1) + (REAL(i,RKIND)-0.5_RKIND)
     &                 * pcellsize - gfarleft(1))/gcellsize,IKIND) + 1
                  k0 = int((pleft(3) + (REAL(j,RKIND)-0.5_RKIND)
     &                 * pcellsize - gfarleft(3))/gcellsize,IKIND) + 1
c
                  sample1 = grid1(i0,n,k0)
c
                  if (iflag .eq. 1) then
                     if (flaggrid(i0,n,k0) .ne. 0._RKIND) 
     &                    sample1 = 0._RKIND
                  endif
c
                  if (ifield .eq. 1) plane(i,j) = plane(i,j) + 
     &                 sample1 * weight1
                  if (ifield .eq. 2) plane(i,j) = plane(i,j) + 
     &                 sqrt(sample1)  *grid2(i0,n,k0)**2 * weight1
                  if (ifield .eq. 3) plane(i,j) = plane(i,j) + 
     &                 sample1**1.5_RKIND *grid2(i0,n,k0)**2 * weight1
c     
               enddo
            enddo
c     
         endif
c
c ----------------------------------------------------------------------
c        project in k-direction
c
         if (projdim .eq. 2) then
c
c           Loop over the projected region of the plane
c
            do j = istart(2)+1, iend(2)
               do i = istart(1)+1, iend(1)
c
                  i0 = int((pleft(1) + (REAL(i,RKIND)-0.5_RKIND)
     &                 * pcellsize - gfarleft(1))/gcellsize,IKIND) + 1
                  j0 = int((pleft(2) + (REAL(j,RKIND)-0.5_RKIND)
     &                 * pcellsize - gfarleft(2))/gcellsize,IKIND) + 1
c
                  sample1 = grid1(i0,j0,n)
c
                  if (iflag .eq. 1) then
                     if (flaggrid(i0,j0,n) .ne. 0._RKIND) 
     &                    sample1 = 0._RKIND
                  endif
c
                  if (ifield .eq. 1) plane(i,j) = plane(i,j) + 
     &                 sample1 * weight1
                  if (ifield .eq. 2) plane(i,j) = plane(i,j) + 
     &                 sqrt(sample1)  *grid2(i0,j0,n)**2 * weight1
                  if (ifield .eq. 3) plane(i,j) = plane(i,j) + 
     &                 sample1**1.5_RKIND *grid2(i0,j0,n)**2 * weight1
c     
               enddo
            enddo
c     
         endif
c          
c
        endif
c
c     Next plane in projection
c
      enddo
c
      return
      end
