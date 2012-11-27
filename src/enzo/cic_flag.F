#include "fortran.def"
c=======================================================================
c//////////////////////  SUBROUTINE CIC_FLAG  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cic_flag(posx, posy, posz, partmass, ndim, npositions, 
     &                    itype, ffield, leftedge, 
     &                    dim1, dim2, dim3, cellsize, imatch1, 
     &                    imatch2, minmassmust, buffersize)
#ifndef CONFIG_PFLOAT_16
c
c  PERFORMS 1/2/3D CLOUD-IN-CELL MARKING OF FLAGGING FIELD
c
c  written by: Greg Bryan
c  date:       September, 2003
c  modified1:
c
c  PURPOSE: This routine performs a 1,2 or 3 dimension setting of a
c           flagging field for must-refine particles.

c
c  INPUTS:
c     ndim       - dimensionality
c     cellsize   - the cell size of field
c     dim1,2,3   - real dimensions of field
c     leftedge   - the left edge(s) of field
c     npositions - number of particles
c     posx,y,z   - particle positions
c     partmass   - particle masses
c     itype      - 1D field (length npositions) of types
c     imatch1,2  - INTG_PREC(s) indicating type of particle to match
c     minmassmust - minimum mass for MustRefine particles to refine around
c     buffersize - how many cells to flag around the particle
c
c  OUTPUT ARGUMENTS: 
c     ffield      - field to be deposited to
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
      INTG_PREC dim1, dim2, dim3, npositions, ndim, imatch1, imatch2
      R_PREC minmassmust
      INTG_PREC buffersize
      P_PREC posx(npositions), posy(npositions), posz(npositions), 
     &        partmass(npositions), leftedge(3)
      P_PREC cellsize
      INTG_PREC itype(npositions)
      INTG_PREC ffield(dim1, dim2, dim3)
c
c  locals
c
      INTG_PREC i1, j1, k1, n, ii, jj, kk
      INTG_PREC istart, iend, jstart, jend, kstart, kend
      R_PREC    xpos, ypos, zpos
      P_PREC edge1, edge2, edge3, half, fact
      parameter (half = 0.5001_PKIND)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      fact = 1._PKIND/cellsize
      edge1 = REAL(dim1,PKIND) - half
      edge2 = REAL(dim2,PKIND) - half
      edge3 = REAL(dim3,PKIND) - half
c
c     1D
c
      if (ndim .eq. 1) then
c
         do n=1, npositions
c
c           only do this for must-refine particles
c
            if (itype(n) .eq. imatch1 .or. itype(n) .eq. imatch2) then
c
c	    only for particles above minimum mass
c
	      if (partmass(n) .gt. minmassmust) then
c
c             Compute the position of the central cell
c
              xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
c
c             Convert this into an INTG_PREC index
c
              i1  = int(xpos + 0.5_RKIND,IKIND)
c
c             set flagging field
c
c             keep default option unrolled!
              if (buffersize .eq. 1) then
                 ffield(i1  ,1,1) = 1
                 ffield(i1+1,1,1) = 1
              else
                 istart = max(i1-buffersize+1, 1)
                 iend = min(i1+buffersize, dim1)
                 do ii = istart, iend
                    ffield(ii,1,1) = 1
                 enddo
              endif
c
            endif
	    endif
c
         enddo
c
      endif
c
c     2D
c
      if (ndim .eq. 2) then
c
         do n=1, npositions
c
c           only do this for must-refine particles
c
            if (itype(n) .eq. imatch1 .or. itype(n) .eq. imatch2) then
c
c	    only for particles above minimum mass
c
	      if (partmass(n) .gt. minmassmust) then
c
c             Compute the position of the central cell
c
              xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
              ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
c
c             Convert this into an INTG_PREC index
c
              i1  = int(xpos + 0.5_RKIND,IKIND)
              j1  = int(ypos + 0.5_RKIND,IKIND)
c
c             Interpolate from field into sumfield
c
c             keep default option unrolled!
              if (buffersize .eq. 1) then
                 ffield(i1  ,j1  ,1) = 1
                 ffield(i1+1,j1  ,1) = 1
                 ffield(i1  ,j1+1,1) = 1
                 ffield(i1+1,j1+1,1) = 1
              else
                 istart = max(i1-buffersize+1, 1)
                 iend = min(i1+buffersize, dim1)
                 jstart = max(j1-buffersize+1, 1)
                 jend = min(j1+buffersize, dim2)
                 do jj = jstart, jend
                    do ii = istart, iend
                       ffield(ii,jj,1) = 1
                    enddo
                 enddo
              endif
c
              endif
            endif
c
         enddo
c
      endif
c
c     3D
c
      if (ndim .eq. 3) then
c
         do n=1, npositions
c
c           only do this for must-refine particles
c
            if (itype(n) .eq. imatch1 .or. itype(n) .eq. imatch2) then
c
c	      only for particles above minimum mass
c
	      if (partmass(n) .gt. minmassmust) then 
c
c             Compute the position of the central cell
c
              xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
              ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
              zpos = min(max((posz(n) - leftedge(3))*fact, half), edge3)
c  
c             Convert this into an INTG_PREC index
c  
              i1  = int(xpos + 0.5_RKIND,IKIND)
              j1  = int(ypos + 0.5_RKIND,IKIND)
              k1  = int(zpos + 0.5_RKIND,IKIND)
c
c             Set flagging field
c     
c             keep default option unrolled!
              if (buffersize .eq. 1) then
                 ffield(i1  ,j1  ,k1  ) = 1
                 ffield(i1+1,j1  ,k1  ) = 1
                 ffield(i1  ,j1+1,k1  ) = 1
                 ffield(i1+1,j1+1,k1  ) = 1
                 ffield(i1  ,j1  ,k1+1) = 1
                 ffield(i1+1,j1  ,k1+1) = 1
                 ffield(i1  ,j1+1,k1+1) = 1
                 ffield(i1+1,j1+1,k1+1) = 1
              else
                 istart = max(i1-buffersize+1, 1)
                 iend = min(i1+buffersize, dim1)
                 jstart = max(j1-buffersize+1, 1)
                 jend = min(j1+buffersize, dim2)
                 kstart = max(k1-buffersize+1, 1)
                 kend = min(k1+buffersize, dim3)
c	         print*,'cic_flag: buffersize is:',buffersize
                 do kk = kstart, kend
                    do jj = jstart, jend
                       do ii = istart, iend
                          ffield(ii,jj,kk) = 1
                       enddo
                    enddo
                 enddo
              endif
c
              endif
	  endif
c
         enddo
c
      endif
c
      return
#endif
      end
