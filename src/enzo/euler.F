#include "fortran.def"
#include "error.def"

c=======================================================================
c//////////////////////////  SUBROUTINE EULER  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine euler(
     &            dslice, eslice, grslice, geslice,
     &            uslice, vslice, wslice, dx, diffcoef,
     &            idim, jdim, i1, i2, j1, j2, dt, 
     &            gamma, idiff, gravity, idual, eta1, eta2,
     &            df, ef, uf, vf, wf, gef, ges,
     &            ncolor, colslice, colf, dfloor
     &                 )
c
c  SOLVES THE EULERIAN CONSERVATION LAWS USING FLUXES FROM THE RIEMANN SOLVER
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  June, 1993 by Greg Bryan (Lagrange+remap --> Eulerian)
c  modified2:  July, 1994 by GB; switched to slicewise and put all the
c                                information into the argument line
c  modified3:  July, 1994 by GB; moved diffusion coefficient to own routine
c  modified4:  Febuary, 2005 by Alexei Kritsuk; fixed a bug in RAREFACTION1
c                            and a bug in update eq. 3.1 with idiff != 0.
c  modified5:  Sept. 2005 by DC, fixed the flux array to work with cosmology 
c                         and AMR
c  modified6:  June, 2010 by JHW; put two-shock flux calculations into 
c                         flux_twoshock()      
c
c  PURPOSE:  Updates the conservation laws in Eulerian form using
c    fluxes in the sweep-direction computed by the Riemann solver.  This
c    versions works on a single two dimensional slice.  It also adds
c    diffusive fluxes, if requested.
c
c  INPUT:
c    diffcoef - diffusion coefficient in slice k
c    dslice - extracted 2d slice of the density, d
c    dt     - timestep in problem time
c    dl,rs  - density at left and right edges of each cell
c    dx     - distance between Eulerian zone edges in sweep direction
c    eslice - extracted 2d slice of the energy, e
c    eta1   - (dual) selection parameter for gas energy (typically ~0.001)
c    eta2   - (dual) selection parameter for total energy (typically ~0.1)
c    gamma  - parameter in ideal gas law
c    geslice - extracted 2d slice of the gas energy, ge
c    gravity - gravity flag (0 = off)
c    grslice - acceleration in this dimension in this slice
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    idiff  - INTG_PREC flag for standard artificial diffusion (0 = off)
c    idual  - dual energy formalism flag (0 = off)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    pl,rs  - pressure at left and right edges of each cell
c    ul,rs  - 1-velocity at left and right edges of each cell
c    uslice - extracted 2d slice of the 1-velocity, u
c    vl,rs  - 2-velocity at left and right edges of each cell
c    vslice - extracted 2d slice of the 2-velocity, v
c    wl,rs  - 3-velocity at left and right edges of each cell
c    wslice - extracted 2d slice of the 3-velocity, w
c    dfloor - density floor applied if non-zero
c
c  OUTPUT:
c    dslice - extracted 2d slice of the density, d
c    geslice - extracted 2d slice of the gas energy, ge
c    eslice - extracted 2d slice of the energy, e
c    uslice - extracted 2d slice of the 1-velocity, u
c    vslice - extracted 2d slice of the 2-velocity, v
c    wslice - extracted 2d slice of the 3-velocity, w
c
c  LOCALS:
c
#define NO_GRAVITY_SECOND_ORDER_CORRECTION
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
      INTG_PREC ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC gravity, i1, i2, idiff, idim, idual, j1, j2, 
     &       jdim, ncolor
      R_PREC    dt, eta1, eta2, gamma, dfloor
      R_PREC diffcoef(idim,jdim),  dslice(idim,jdim),      dx(idim   ),
     &       eslice(idim,jdim), grslice(idim,jdim), geslice(idim,jdim),
     &       uslice(idim,jdim),  vslice(idim,jdim),  wslice(idim,jdim)
      R_PREC df(idim,jdim),      ef(idim,jdim),      uf(idim,jdim),
     &     vf(idim,jdim),      wf(idim,jdim),     gef(idim,jdim),
     &     ges(idim,jdim)
      R_PREC colslice(idim,jdim,ncolor), colf(idim,jdim,ncolor)
c
c  local declarations
c
      INTG_PREC i, j, n
      R_PREC qa, qb, qc, frac, uadvect, eold
      R_PREC dnu(ijkn), dnuinv(ijkn), uold(ijkn), pcent(ijkn), 
     &     eratio(ijkn), dadx(ijkn),  dddx(ijkn)
      R_PREC min_color
      parameter(min_color=1e-5_RKIND*tiny)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     write(6,*) 'EULER: dt =',dt,' idiff =',idiff
c     write(6,*) 'EULER: idim =',idim,' jdim =',jdim
c     write(6,*) 'EULER: i1   =',i1,  ' i2   =',i2
c     write(6,*) 'EULER: j1   =',j1,  ' j2   =',j2
c
c  Set constants
c
      qa = (gamma + 1._RKIND)/(2._RKIND*gamma)
      qb = (gamma - 1._RKIND)/(gamma + 1._RKIND)
c
c  Loop over sweep lines (in this slice)
c
      do j=j1, j2
c     
c  Update conservation laws  (eq. 3.1)
c
       do i=i1, i2
          qc = dt/dx(i)
          dnu(i) = dslice(i,j) + (df(i,j) - df(i+1,j))
          dnuinv(i) = 1._RKIND/dnu(i)
c          write(21,*) i, dnu(i), dub(i), db(i), ub(i)
c          uadvect = (dnu(i)-dslice(i,j))/(qc*(db(i+1)-db(i)))
c          if (uslice(i,j)/(uadvect) .gt. 1.8 .and. uadvect .ne. 0.0)
c     &        write(6,*) i,j,uadvect,uslice(i,j),dt
c          if (dnu(i) .le. 0.0 .or. eslice(i,j)*dslice(i,j) +
c     &      qc*(dueb(i)-dueb(i+1)+upb(i)-upb(i+1)) .lt. 0.0) then


c
c  Apply density floor (if used) and update inv density
c
          if (dfloor .gt. 0._RKIND) then
              dnu(i)    = max( dnu(i), dfloor)
              dnuinv(i) = 1._RKIND / dnu(i)
          endif

          uold(i)     = uslice(i,j)
c          write(20,*) i, uslice(i,j), (uslice(i,j)*dslice(i,j) +
c     &         qc*(duub(i)-duub(i+1) +  pb(i) -  pb(i+1)))/dnu(i),
c     &               dnu(i), dslice(i,j)
          uslice(i,j) = (uslice(i,j)*dslice(i,j) +
     $         (uf(i,j) - uf(i+1,j))) * dnuinv(i)
          vslice(i,j) = (vslice(i,j)*dslice(i,j) +
     $         (vf(i,j) - vf(i+1,j))) * dnuinv(i)
          wslice(i,j) = (wslice(i,j)*dslice(i,j) +
     $         (wf(i,j) - wf(i+1,j))) * dnuinv(i)

c          if ((eslice(i,j)*dslice(i,j) +
c     &         qc*(dueb(i)-dueb(i+1) + upb(i) - upb(i+1)))/dnu(i) .lt.
c     &         0.5*eslice(i,j)) then
c            write(6,*) i,j,qc*(dueb(i)-dueb(i+1)),
c     &                     qc*(upb(i)-upb(i+1)), eslice(i,j)*dslice(i,j)
c            write(6,*) dueb(i),dueb(i+1),db(i),db(i+1)
c            write(6,*) eb(i),eb(i+1),ub(i),ub(i+1)
c          endif
          eold = eslice(i,j)
          eslice(i,j) = max(0.1_RKIND*eslice(i,j),
     &                  (eslice(i,j)*dslice(i,j) +
     &         (ef(i,j) - ef(i+1,j))) * dnuinv(i))
       enddo
c
c      Colour variables (note: colf already multiplied by dt)
c
       do n=1,ncolor
          do i=i1, i2
c             colslice(i,j,n) = (colslice(i,j,n)*dslice(i,j) +
c     &         (colf(i,j,n)-colf(i+1,j,n))/dx(i)      )/dnu(i) ! c*d conserved
             colslice(i,j,n) =  colslice(i,j,n)              +
     &            (colf(i,j,n)-colf(i+1,j,n))/dx(i) ! c conserved
             colslice(i,j,n) = max(colslice(i,j,n), min_color)
          enddo
       enddo
c
c      Conservation law for gas energy, if using the dual energy formalism
c           (this includes both the flux term and a source term - yuck).
c         Here, we compute the ratio of thermal energies derived the
c            two different ways and then use that ratio to choose how to
c            compute the pressure at the center of zone i.  This is needed
c            for the source term in the gas energy equation.
c
       if (idual .eq. 1) then
          do i=i1, i2
             if (geslice(i,j) .lt. 0._RKIND) 
     $            write(6,*) 'euler: geslice<0:', i,j,geslice(i,j)
             geslice(i,j) = max((geslice(i,j)*dslice(i,j) +
     $            (gef(i,j) - gef(i+1,j)) + ges(i,j)) * dnuinv(i)
     $            ,0.5_RKIND*geslice(i,j))
c             if (geslice(i,j) .lt. tiny) geslice(i,j) = tiny
             if (geslice(i,j) .lt. 0._RKIND) then
                write(6,*) i,j,dslice(i,j),dnu(i),pcent(i),eslice(i,j),
     &               gef(i,j), gef(i+1,j), ges(i,j),
     &               geslice(i,j)*dslice(i,j),
     &               uslice(i,j),uslice(i+1,j)
c                write(6,*) pb(i),pb(i+1),p0(i),p0(i+1),
c     &               pbar(i,j),pbar(i+1,j),lbar(i),l0(i),
c     &               pls(i,j),prs(i,j),pls(i+1,j),prs(i+1,j)
c                write(6,*) ub(i),ub(i+1),geb(i),geb(i+1),
c     &               db(i),db(i+1),dub(i),dub(i+1)
                write(0,*) 'stop in euler with geslice < 0'
                ERROR_MESSAGE
             endif
          enddo
       endif
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
     &            0.5_RKIND*dddx(i), min(dslice(i,j), dslice(i-1,j))))
             dddx(i) = 2._RKIND*(-dslice(i,j) + max(dslice(i,j) +
     &            0.5_RKIND*dddx(i), min(dslice(i,j), dslice(i+1,j))))
          enddo
c
          do i=i1, i2
             grslice(i,j) = grslice(i,j) + 
     &            0.5_RKIND*dadx(i)*dddx(i)/(12._RKIND*dslice(i,j))
          enddo
       endif
c
#endif /* GRAVITY_SECOND_ORDER_CORRECTION */
c
c  If there is gravity, add the gravity terms here (eq. 3.1 or 3.8).
c    (Note: the acceleration is already time-centered).
c
       if (gravity .eq. 1) then
          do i=i1, i2

#define GRAVITY_METHOD1

#ifdef GRAVITY_METHOD1 /* (right) */
             uslice(i,j) = uslice(i,j) + dt*grslice(i,j)*0.5_RKIND*
     &            (dslice(i,j)*dnuinv(i)+1._RKIND)
             eslice(i,j) = eslice(i,j) + dt*grslice(i,j)*0.5_RKIND*
     &            (uslice(i,j) + uold(i)*dslice(i,j)*dnuinv(i))
             if (eslice(i,j) .le. 0) write(6,*) 'eu1',i,j,eslice(i,j),
     &            dslice(i,j),dnu(i),uslice(i,j),grslice(i,j),
     &            dt*grslice(i,j)*0.5_RKIND*
     &            (dslice(i,j)/dnu(i)+1._RKIND),
     &            dt*grslice(i,j)*0.5_RKIND*
     &            (uslice(i,j) + uold(i)*dslice(i,j)/dnu(i))
             eslice(i,j) = max(eslice(i,j), tiny)
#endif /* GRAVITY_METHOD1 */

#ifdef GRAVITY_METHOD2 /* (wrong) */
             uslice(i,j) = uslice(i,j) + dt*grslice(i,j)
             eslice(i,j) = eslice(i,j) + dt*grslice(i,j)*0.5_RKIND*
     &                                      (uslice(i,j) + uold(i))
#endif /* GRAVITY_METHOD2 */
          enddo
       endif
c
c  Update the new density
c
       do i=i1, i2
          dslice(i,j) = dnu(i)
       enddo
c
      enddo
c
      return
      end

