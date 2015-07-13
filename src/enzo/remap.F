#include "fortran.def"
#include "error.def"

c=======================================================================
c///////////////////////////  SUBROUTINE REMAP  \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine remap(
     &          dslice, eslice, uslice, vslice, wslice, 
     &          geslice, dxnslice, dx,
     &          idim, jdim, i1, i2, j1, j2, idual,
     &          df, ef, uf, vf, wf, gef
     &                )
c
c  PERFORMS REMAP IN SWEEP-DIRECTION
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  November, 1994 by Greg Bryan; switched to slicewise
c  modified2:
c
c  PURPOSE:  Remaps the fundamental variables d,e,u,v,w back to the
c    Eulerian mesh after Lagrangean update.  Fluxes for each variable
c    computed in INTRMP are used.
c
c  INPUTS:
c    dslice - extracted 2d slice of the density, d
c    dxnslice - distance between Lagrangean zone edges in sweep direction
c    dx     - distance between Eulerian zone edges in sweep direction
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    uslice - extracted 2d slice of the 1-velocity, u
c    vslice - extracted 2d slice of the 2-velocity, v
c    wslice - extracted 2d slice of the 3-velocity, w
c    df     - density flux
c    ef     - total specific energy flux (e*dx*d)
c    uf     - 1-direction momentum flux (u*dx*d)
c    vf     - 2-direction momentum flux (u*dx*d)
c    wf     - 3-direction momentum flux (u*dx*d)
c
c  OUTPUTS:
c    dslice - updated density slice
c    eslice - updated total specific energy slice
c    uslice - updated velocity-1 slice
c    vslice - updated velocity-2 slice
c    wslice - updated velocity-3 slice
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
      INTG_PREC ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC i1, i2, idim, idual, j1, j2, jdim
      R_PREC  dslice(idim,jdim),      dx(idim     ),  eslice(idim,jdim),
     &        uslice(idim,jdim),  vslice(idim,jdim),  wslice(idim,jdim),
     &        dxnslice(idim,jdim), geslice(idim,jdim)
      R_PREC  df(idim,jdim),      ef(idim,jdim),      uf(idim,jdim),
     &        vf(idim,jdim),      wf(idim,jdim),     gef(idim,jdim)
c
c  local declarations
c
      INTG_PREC i, j
      R_PREC    dm(ijkn), dmnu(ijkn)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Update zones, based on fluxes in/out of cell
c
      do j=j1, j2
c
         do i=i1, i2
            dm(i)       = dslice(i,j)*dxnslice(i,j)
            dmnu(i)     = dm(i) - (df(i+1,j)-df(i,j))
            if (dmnu(i) .le. 0._RKIND) then
               write(6,*) 'remap',i,j,i2,dslice(i,j),dxnslice(i,j)
               write(6,*) df(i+1,j),df(i,j)
               write(6,*) dslice(i-1,j),dslice(i,j),dslice(i+1,j)
               write(6,*) uslice(i-1,j),uslice(i,j),uslice(i+1,j)
               ERROR_MESSAGE
            endif
            dslice(i,j) = dmnu(i)/dx(i)
            eslice(i,j) = (dm(i)*eslice(i,j) - (ef(i+1,j)-ef(i,j)))/
     &                    dmnu(i)
            uslice(i,j) = (dm(i)*uslice(i,j) - (uf(i+1,j)-uf(i,j)))/
     &                    dmnu(i)
            vslice(i,j) = (dm(i)*vslice(i,j) - (vf(i+1,j)-vf(i,j)))/
     &                    dmnu(i)
            wslice(i,j) = (dm(i)*wslice(i,j) - (wf(i+1,j)-wf(i,j)))/
     &                    dmnu(i)
         enddo
c
         if (idual .eq. 1) then
            do i=i1, i2
               if (geslice(i,j) .le. 0._RKIND) write(6,*) 'remap_a',i
               geslice(i,j) = max(0.5_RKIND*geslice(i,j), 
     &                         (dm(i)*geslice(i,j) - 
     &                         (gef(i+1,j)-gef(i,j)))/ dmnu(i) )
               if (geslice(i,j) .le. 0._RKIND) write(6,*) 'remap_b',i
            enddo
         endif
c
      enddo
c
      return
      end
