#include "fortran.def"
#include "error.def"

c=======================================================================
c///////////////////////////  SUBROUTINE LGRG  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine lgrg(
     &            dslice, eslice, geslice, uslice, dx, diffcoef,
     &            idim, jdim, i1, i2, j1, j2, dt, gamma, idiff, 
     &            idual, etol2, vslice, wslice,
     &            dls, drs, pls, prs, uls, urs, pbar, ubar,
     &            dxlslice, dxnslice, gravity, grslice
     &               )
c
c  SOLVES THE LAGRANGEAN CONSERVATION LAWS USING FLUXES FROM THE RIEMASS SOLVER
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  November, 1994 by Greg Bryan; switch to slicewise and put all 
c                     the information into the argument line
c  modified2:  November, 1994 by GB; moved diffusion coefficient to own routine
c  modified3:
c
c  PURPOSE:  Updates the conservation laws in Lagrangean form using
c    fluxes in the sweep-direction computed by the Riemann solver.
c    This versions works on a single two dimensional slice.  It also adds
c    diffusive fluxes, if requested.
c
c  INPUTS:
c    diffcoef - diffusion coefficient in slice k
c    dslice - extracted 2d slice of the density, d
c    dt     - timestep in problem time
c    dl,rs  - density at left and right edges of each cell
c    dx     - distance between Eulerian zone edges in sweep direction
c    eslice - extracted 2d slice of the energy, e
c    gamma  - parameter in ideal gas law
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    idiff  - INTG_PREC flag for standard artificial diffusion (0 = off)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    pl,rs  - pressure at left and right edges of each cell
c    ul,rs  - 1-velocity at left and right edges of each cell
c    uslice - extracted 2d slice of the 1-velocity, u
c
c  OUTPUT:
c    dslice - extracted 2d slice of the density, d
c    eslice - extracted 2d slice of the energy, e
c    uslice - extracted 2d slice of the 1-velocity, u
c
c  EXTERNALS:
c
#define NO_GRAVITY_SECOND_ORDER_CORRECTION
c
c  LOCALS:
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
      INTG_PREC gravity, i1, i2, idiff, idim, idual, j1, j2, jdim
      R_PREC  dt, etol2, gamma
      R_PREC  diffcoef(idim,jdim),  dslice(idim,jdim),      dx(idim   ),
     &        eslice(idim,jdim),  uslice(idim,jdim), grslice(idim,jdim),
     &        dxlslice(idim,jdim), dxnslice(idim,jdim), 
     &        geslice(idim,jdim)
      R_PREC  dls(idim,jdim),      drs(idim,jdim),  vslice(idim,jdim),
     &        pbar(idim,jdim),      pls(idim,jdim),  wslice(idim,jdim),
     &        prs(idim,jdim),     ubar(idim,jdim),
     &        uls(idim,jdim),      urs(idim,jdim)
c
c  local declarations
c
      INTG_PREC i, j
      R_PREC    qa, qb, qc, nu, di(ijkn), upb(ijkn), erat(ijkn),
     &        dadx(ijkn), dddx(ijkn), pcent(ijkn)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
      do j=j1, j2
c
c compute ratio of thermal to total energy
c
      if (idual .eq. 1) then
         do i=i1-1, i2+1
            erat(i) = (eslice(i,j) - 0.5_RKIND*(uslice(i,j)**2 +
     &                                    vslice(i,j)**2 +
     &                                    wslice(i,j)**2) )/eslice(i,j)
         enddo
      endif
c
c  Compute diffusive fluxes (using coefficients from diffcoef)
c
      if (idiff .ne. 0) then
         di(i1-1) = 1._RKIND/dslice(i1-1,j)
         do i=i1, i2+1
            nu         = min(dslice(i,j),dslice(i-1,j))*diffcoef(i,j)
            di(i)      = 1._RKIND/dslice(i,j)
            ubar (i,j) = ubar(i,j)           + nu*(di(i  )-di(i-1  ))
            pbar (i,j) = pbar(i,j)           -
     &                                 nu*(uslice(i,j)-uslice(i-1,j))
            upb(i    ) = ubar(i,j)*pbar(i,j) -
     &                                 nu*(eslice(i,j)-eslice(i-1,j))
         enddo
      else
         do i=i1, i2+1
            upb(i)     = ubar(i,j)*pbar(i,j)
         enddo
      endif
c
c Compute pcent, the cell-centered "time-averaged" pressure.
c   Perform Riemann correction if necessary.
c
      if (idual .eq. 1) then
         do i=i1, i2
            pcent(i) = max((gamma - 1._RKIND)*geslice(i,j)*dslice(i,j),
     &                 0.5_RKIND*(pbar(i,j)+pbar(i+1,j)))
            if (max(erat(i-1),erat(i),erat(i+1)) .gt. etol2)
     &          pcent(i) = (gamma - 1._RKIND)*geslice(i,j)*dslice(i,j)
            pcent(i) = min(pcent(i), 
     &                     0.5_RKIND*geslice(i,j)*dslice(i,j)*dx(i)/
     &                     (max(ubar(i+1,j)-ubar(i,j), tiny)*dt))
         enddo
      endif
c
c  Compute new Lagrangean grid positions
c
      do i=i1, i2
         dxnslice(i,j) = dx(i) + dt*(ubar(i+1,j) - ubar(i,j))
         if (dxnslice(i,j) .le. 0._RKIND) then
            write(6,*) 'lgrg1',i,j,i2,j2,dx(i),dt
            write(6,*) dxnslice(i-1,j),dxnslice(i,j),dxnslice(i+1,j)
            write(6,*) ubar(i-1,j),ubar(i+1,j),ubar(i,j)
            write(6,*) dslice(i-1,j),dslice(i,j),dslice(i+1,j)
            write(6,*) uslice(i-1,j),uslice(i,j),uslice(i+1,j)
            write(6,*) eslice(i-1,j),eslice(i,j),eslice(i+1,j)
            write(6,*) geslice(i-1,j),geslice(i,j),geslice(i+1,j)
            ERROR_MESSAGE
         endif
      enddo
c
c  Compute distance each cell edge has travelled
c
      do i=i1, i2+1
         dxlslice(i,j) = dt*ubar(i,j)
      enddo
c
c  Set boundary values of dxn/dxl.  This should be improved.
c
      do i=i1-3, i1-1
         dxnslice(i,j) = dx(i) + dt*(ubar(i1+1,j) - ubar(i1,j))
      enddo
c
      do i=i2+1, i2+3
         dxnslice(i,j) = dx(i) + dt*(ubar(i2+1,j) - ubar(i2,j))
      enddo
c
c  If there is gravity, the compute the second order correction to the
c   acceleration due a slope in both density and acceleration.
c
#ifdef GRAVITY_SECOND_ORDER_CORRECTION
c
c      Compute slopes and enforce limited monotonoctiy on dddx
c
      if (gravity .eq. 1) then
         do i=i1,i2
            dadx(i) = grslice(i+1,j) - grslice(i-1,j)
            dddx(i) =  dslice(i+1,j) -  dslice(i-1,j)
c     
            dddx(i) = 2._RKIND*( dslice(i,j) - max(dslice(i,j) -
     &           0.5_RKIND*dddx(i), min(dslice(i,j), dslice(i-1,j))))
            dddx(i) = 2._RKIND*(-dslice(i,j) + max(dslice(i,j) +
     &           0.5_RKIND*dddx(i), min(dslice(i,j), dslice(i+1,j))))
         enddo
c
         do i=i1, i2
            grslice(i,j) = grslice(i,j) + 
     &                     0.5_RKIND*dadx(i)*dddx(i)/(12.0*dslice(i,j))
         enddo
      endif
c
#endif /* GRAVITY_SECOND_ORDER_CORRECTION */
c
c  Update conservation laws
c
      do i=i1, i2
        qa = dt/(dslice(i,j)*dx(i))
        eslice(i,j) = eslice(i,j) - qa*(upb (i+1  ) - upb (i  ))
        uslice(i,j) = uslice(i,j) - qa*(pbar(i+1,j) - pbar(i,j))
      enddo
c
c  Update internal energy source term
c
      if (idual .eq. 1) then
         do i=i1, i2
            geslice(i,j) = geslice(i,j) - dt*pcent(i)/dslice(i,j)*
     &                       (ubar(i+1,j) - ubar(i,j))/dx(i)
         enddo
      endif
c
c  Add gravity
c
      if (gravity .eq. 1) then
         do i=i1, i2
            eslice(i,j) = eslice(i,j) + dt*grslice(i,j)*(uslice(i,j) +
     &                                  0.5_RKIND*dt*grslice(i,j))
            uslice(i,j) = uslice(i,j) + dt*grslice(i,j)
         enddo
      endif
c
c  Set new density
c
      do i=i1, i2
         dslice(i,j) = dslice(i,j)*dx(i)/dxnslice(i,j)
         if (dslice(i,j) .le. 0._RKIND) then
            write(6,*) 'lgrg',i,j,i2,dslice(i,j),dx(i),dxnslice(i,j)
            write(6,*) ubar(i-1,j),ubar(i+1,j),ubar(i,j)
            write(6,*) dslice(i-1,j),dslice(i,j),dslice(i+1,j)
            write(6,*) uslice(i-1,j),uslice(i,j),uslice(i+1,j)
            ERROR_MESSAGE
         endif
      enddo
c
c  Next slice
c
      enddo
c
      return
      end
