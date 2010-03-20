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
     &            dls, drs, pls, prs, gels, gers,
     &            uls, urs, vls, vrs, wls, wrs,
     &            pbar, ubar,
     &            df, ef, uf, vf, wf, gef,
     &            ncolor, colslice, colls, colrs, colf
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
c  modified5:  Sept. 2005 by DC, fixed the flux array to work with cosmology and AMR
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
c    idiff  - integer flag for standard artificial diffusion (0 = off)
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
c
      integer ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, i1, i2, idiff, idim, idual, j1, j2, jdim, ncolor
      real    dt, eta1, eta2, gamma
      real diffcoef(idim,jdim),  dslice(idim,jdim),      dx(idim     ),
     &       eslice(idim,jdim), grslice(idim,jdim), geslice(idim,jdim),
     &       uslice(idim,jdim),  vslice(idim,jdim),  wslice(idim,jdim)
      real      dls(idim,jdim),     drs(idim,jdim),
     &         pbar(idim,jdim),     pls(idim,jdim),    gels(idim,jdim),
     &          prs(idim,jdim),    ubar(idim,jdim),    gers(idim,jdim),
     &          uls(idim,jdim),     urs(idim,jdim),     vls(idim,jdim),
     &          vrs(idim,jdim),     wls(idim,jdim),     wrs(idim,jdim),
     &           df(idim,jdim),      ef(idim,jdim),      uf(idim,jdim),
     &           vf(idim,jdim),      wf(idim,jdim),     gef(idim,jdim)
      real colslice(idim,jdim,ncolor),   colls(idim,jdim,ncolor),
     &        colrs(idim,jdim,ncolor),    colf(idim,jdim,ncolor)
c
c  local declarations
c
      integer i, j, n
      real      alpha      ,    c0(ijkn),    cb(ijkn),  cbar(ijkn),
     &             d0(ijkn),    db(ijkn),  dbar(ijkn),   dnu(ijkn),
     &            dub(ijkn),  duub(ijkn),  duvb(ijkn),  duwb(ijkn),
     &           dueb(ijkn),    eb(ijkn),    l0(ijkn),  lbar(ijkn),
     &             p0(ijkn),   geb(ijkn), dugeb(ijkn),    pb(ijkn),
     &             qa      ,    qb      ,    qc      ,    sn(ijkn),
     &             u0(ijkn),    ub(ijkn),   upb(ijkn),    vb(ijkn),
     &             wb(ijkn),    z0(ijkn),
     &           uold(ijkn), pcent(ijkn),eratio(ijkn),  frac      ,
     &           dadx(ijkn),  dddx(ijkn), colb(ijkn,MAX_COLOR)
      real dele, delu, uleft, uright, uadvect, eold
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
      qa = (gamma + 1.0)/(2.0*gamma)
      qb = (gamma - 1.0)/(gamma + 1.0)
c
c  Loop over sweep lines (in this slice)
c
      do j=j1, j2
c
c  Evaluate time-averaged quantities 
c   (see Colella, Siam J Sci Stat Comput 1982, 3, 77.  Appendix)
c
       do i=i1, i2+1
c
          sn(i)    = sign(1.0, -ubar(i,j))
c
c  Collect values of interest depending on which way fluid is flowing
c
          if (sn(i) .lt. 0.0) then
             u0(i) = uls(i,j)
             p0(i) = pls(i,j)
             d0(i) = dls(i,j)
c             vb(i) = vls(i,j)
c             wb(i) = wls(i,j)
c             geb(i) = gels(i,j)
          else
             u0(i) = urs(i,j)
             p0(i) = prs(i,j)
             d0(i) = drs(i,j)
c             vb(i) = vrs(i,j)
c             wb(i) = wrs(i,j)
c             geb(i) = gers(i,j)
          endif
          c0(i)    = sqrt(max(gamma*p0(i)/d0(i), tiny))
          z0(i)    = c0(i)*d0(i)*sqrt(max(
     &               1.0 + qa*(pbar(i,j)/p0(i)-1.0), tiny))
c          write(24,*) i, u0(i), uls(i,j), urs(i,j), c0(i)
       enddo
c
c      Repeat for colour variables (now moved below)
c
#ifdef UNUSED
       do n=1,ncolor
          do i=i1, i2+1
             if (sn(i) .lt. 0.0) then
                colb(i,n) = colls(i,j,n)
             else
                colb(i,n) = colrs(i,j,n)
             endif
          enddo
       enddo
#endif
c
c  Compute equivalent bar (inside shock & rarefaction) values for density
c    and sound speed
c
       do i=i1, i2+1
          dbar(i)  = 1.0/(1.0/d0(i) - (pbar(i,j)-p0(i))/
     &                                max(z0(i)**2, tiny))
          cbar(i)  = sqrt(max(gamma*pbar(i,j)/dbar(i), tiny))
       enddo
c
c
c  Find lambda values for the shock and rarefaction
c
       do i=i1, i2+1
         if (pbar(i,j) .lt. p0(i)) then
            l0(i)   = u0(i)*sn(i) + c0(i)
            lbar(i) = sn(i)*ubar(i,j)+cbar(i)
         else
            l0(i)   = u0(i)*sn(i) + z0(i)/d0(i)
c            l0(i)   = ubar(i,j)*sn(i) + z0(i)/d0(i)
            lbar(i) = l0(i)
         endif
c
c  Compute values for inside a rarefaction fan
c     (method described in Colella, 1982)
c
#define RAREFACTION2 //AK
c
#ifdef RAREFACTION0
c
c  (quasi-)isothermal hydrodynamics => use RAREFACTION0
c  sound speed is assumed to be 1.0 in expression for db(i)
c  (see eqs (3.11) - (3.16) in Balsara (1994) ApJ 420, p. 200)
c
         ub(i)    = -cbar(i)*sn(i)
         db(i)    = d0(i)*exp(-sn(i)*u0(i) - 1.0)
         pb(i)    = db(i)*cbar(i)**2
c
#endif /* RAREFACTION0 */
c
#ifdef RAREFACTION1
c
c  Compute values for inside a rarefaction fan
c  (method described in Colella (1982) SIAM J Sci Stat Comp 3, p. 107)
c
         cb(i)    = (-sn(i)*ubar(i,j) + 2.0*cbar(i)/(gamma - 1.0))*qb
         ub(i)    = -cb(i)*sn(i)
c         alpha    = (c0(i)**2 - cb(i)**2)/(c0(i)**2 + cb(i)**2)/
c     &              (gamma - 1.0)
c         db(i)    = d0(i)*(1.0 - alpha)/(1.0 + alpha)
         db(i)    = d0(i)*(cb(i)/c0(i))**(2.0/(gamma - 1.0)) !AK
         pb(i)    = db(i)*cb(i)**2/gamma
c
#endif /* RAREFACTION1 */
c
#ifdef RAREFACTION2
c
c  Compute values for inside a rarefaction fan
c     (linear interpolation between end states, as suggested in PPM ref)
c
         frac     = l0(i) - lbar(i)
         if (frac .lt. tiny) frac = tiny
         if (frac .gt. 1.0)  frac = 1.0
         frac     = (0.0 - lbar(i))/frac
         frac     = min(max(frac, 0.0), 1.0)
         pb(i)    = p0(i)*frac + pbar(i,j)*(1.0 - frac)
         db(i)    = d0(i)*frac + dbar(i  )*(1.0 - frac)
         ub(i)    = u0(i)*frac + ubar(i,j)*(1.0 - frac)
c
#endif /* RAREFACTION2 */
c
       enddo
c
c  Cull appropriate states depending on where eulerian position is in solution
c    (lbar >= 0 --> inside post-shock region,
c     l0   <  0 --> outside shock/rarefaction wave,
c     otherwise --> inside rarefaction wave).
c
       do i=i1, i2+1
         if (lbar(i) .ge. 0.0) then !AK
            pb(i) = pbar(i,j)
            db(i) = dbar(i  )
            ub(i) = ubar(i,j)
         endif
         if (l0  (i) .lt. 0.0) then !AK
            pb(i) = p0(i)
            db(i) = d0(i)
            ub(i) = u0(i)
         endif
c         write(22,*) i, ub(i), u0(i), ubar(i,j), sn(i)
c
c         if (db(i) .gt. dbar(i) .and. db(i) .gt. d0(i))
c     &      write(6,1000) i, j, db(i), dbar(i), d0(i)
c 1000  format('euler: ', 2(i4), 3(e12.4))
c
       enddo
c
#define NO_CAVITATE_AVOID
#ifdef CAVITATE_AVOID
       do i=i1, i2+1
         if (ubar(i,j)*ub(i) .lt. 0.0 .and. 
     &       -db(i)*ub(i)*dt/dx(i) .gt. dslice(i,j) .or.
     &        db(i)*ub(i)*dt/dx(i) .gt. dslice(i-1,j)     ) then
            write(6,1010) i,j,ubar(i,j),ub(i),db(i),
     &              dslice(i,j),dslice(i-1,j)
 1010       format('cav:',2i6,8g12.4)
            ub(i) = 0.0
         endif
      enddo
#endif
c
c  Collect values of interest depending on which way fluid is flowing
c     (note: new placement for this - uses ub instead of ubar for consistency)
c
       do i=i1, i2+1
          if (ub(i) .gt. 0.0) then !AK
             vb(i) = vls(i,j)
             wb(i) = wls(i,j)
             geb(i) = gels(i,j)
          else
             vb(i) = vrs(i,j)
             wb(i) = wrs(i,j)
             geb(i) = gers(i,j)
          endif
       enddo
c
c      Repeat for colour variables
c
       do n=1,ncolor
          do i=i1, i2+1
             if (ub(i) .gt. 0.0) then
                colb(i,n) = colls(i,j,n)
             else
                colb(i,n) = colrs(i,j,n)
             endif
          enddo
       enddo
c
c  Dual energy formalism: if sound speed squared is less than eta1*v^2 
c    then discard pbar,dbar,ubar in favour of p0,d0,u0.  This amounts
c     to assuming that we are outside the shocked region but the flow is
c     hypersonic so this should be true.  This is inserted because the
c     mechanism above occasionally fails in hypersonic flows.
c
#ifdef UNUSED
       if (idual .eq. 1) then
          do i=i1, i2+1
             if (gamma*pbar(i,j)/min(dls(i,j),drs(i,j)) .lt. 
     &                                         eta1*ubar(i,j)**2) then
                 pb(i) = p0(i)
                 db(i) = d0(i)
                 ub(i) = u0(i)
              endif
           enddo
        endif
#endif /* UNUSED */
c
c  Calculate total specific energy corresponding to this state
c     (and the specific gas energy).
c
       do i=i1, i2+1
         eb(i) = pb(i)/((gamma-1.0)*db(i)) + 
     &           0.5*(ub(i)**2 + vb(i)**2 + wb(i)**2)
c         geb(i) = pb(i)/((gamma-1.0)*db(i))
       enddo
c
c  Compute terms in differenced hydro equations (eq. 3.1)
c
       if (idiff .ne. 0) then
c
c     ...with diffusion
c
         do i=i1,i2+1
            upb(i)  =  pb(i)*ub(i)
c
            dub(i)  =  ub(i)*db(i)  !AK
c
            duub(i) = dub(i)*ub(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*uslice(i-1,j) - dslice(i,j)*uslice(i,j))
c
            duvb(i) = dub(i)*vb(i)
            duwb(i) = dub(i)*wb(i)
c
c     (should we add diffusion to the cross velocities?  I doubt it)
c     I don't. This doubt kills Noh test problem in 2D at high resolution.
c     Diffusion has to be added to cross velocities. !AK        May 2005.
c
            duvb(i) = dub(i)*vb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*vslice(i-1,j) - dslice(i,j)*vslice(i,j))
            duwb(i) = dub(i)*wb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*wslice(i-1,j) - dslice(i,j)*wslice(i,j))
c
            dueb(i) = dub(i)*eb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*eslice(i-1,j) - dslice(i,j)*eslice(i,j))
c
c	This update must be the last !AK
c
            dub(i)  =  dub(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)               - dslice(i,j)            )
         enddo
c
c        If using dual energy formalism, compute dugeb
c
         if (idual .eq. 1) then
            do i=i1,i2+1
              dugeb(i) = dub(i)*geb(i) + diffcoef(i,j)*
     &         (dslice(i-1,j)*geslice(i-1,j) - dslice(i,j)*geslice(i,j))
            enddo
         endif
c
       else
c
c     ...and without
c
         do i=i1, i2+1
            upb(i)  =  pb(i)*ub(i)
            dub(i)  =  ub(i)*db(i)
            duub(i) = dub(i)*ub(i)
            duvb(i) = dub(i)*vb(i)
            duwb(i) = dub(i)*wb(i)
            dueb(i) = dub(i)*eb(i)
         enddo
c
         if (idual .eq. 1) then
            do i=i1, i2+1
               dugeb(i) = dub(i)*geb(i)
            enddo
         endif
c
       endif
c
c  Copy into flux slices (to return to caller)
c

!      do i=i1, i2+1
!         df(i,j) = dt*dub(i)
!         ef(i,j) = dt*(dueb(i) + upb(i))
!         uf(i,j) = dt*(duub(i) + pb(i))
!         vf(i,j) = dt*duvb(i)
!         wf(i,j) = dt*duwb(i)
!      enddo

       do i=i1, i2+1
          qc = dt/dx(i)
          df(i,j) = qc*dub(i)
          ef(i,j) = qc*(dueb(i) + upb(i))
          uf(i,j) = qc*(duub(i) + pb(i))
          vf(i,j) = qc*duvb(i)
          wf(i,j) = qc*duwb(i)
       enddo

c
       do n=1,ncolor
          do i=i1, i2+1
c             colf(i,j,n) = dt*dub(i)*colb(i,n)   ! color*dens conserved
             colf(i,j,n) = dt*ub(i)*colb(i,n)     ! color      conserved
          enddo
       enddo
c
c      Do the same for the gas energy if using the dual energy formalism
c         (note that we do not include the source term)
c
       if (idual .eq. 1) then
          do i=i1, i2+1
             gef(i,j) = dt*dugeb(i)
          enddo
       endif
c     
c  Update conservation laws  (eq. 3.1)
c
       do i=i1, i2
          qc = dt/dx(i)
          dnu(i)      =              dslice(i,j) +
     &                     qc*( dub(i)- dub(i+1))
c          write(21,*) i, dnu(i), dub(i), db(i), ub(i)
c          uadvect = (dnu(i)-dslice(i,j))/(qc*(db(i+1)-db(i)))
c          if (uslice(i,j)/(uadvect) .gt. 1.8 .and. uadvect .ne. 0.0)
c     &        write(6,*) i,j,uadvect,uslice(i,j),dt
c          if (dnu(i) .le. 0.0 .or. eslice(i,j)*dslice(i,j) +
c     &      qc*(dueb(i)-dueb(i+1)+upb(i)-upb(i+1)) .lt. 0.0) then
          if (dnu(i) .le. 0.0 .or. eslice(i,j) .lt. 0.0) then
             if (dnu(i) .le. 0.0) then
                write (6,*) 'euler: dnu <= 0:', i, j
             else
                write (6,*) 'euler: eslice < 0: ', i, j
                write (6,*) qc*(dueb(i)-dueb(i+1)),
     &                     qc*(upb(i)-upb(i+1)), eslice(i,j)*dslice(i,j)
                write (6,*) eb(i-1), eb(i), eb(i+1)
                write (6,*) dueb(i-1), dueb(i), dueb(i+1)
                write (6,*) upb(i-1), upb(i), upb(i+1)
             endif
             write (6,*) 'a', dnu(i), dslice(i,j), qc, dt, lbar(i),l0(i)
             write (6,*) 'b', db(i-1), db(i), db(i+1)
             write (6,*) 'c', ub(i-1), ub(i), ub(i+1)
             write (6,*) 'd', pb(i-1), pb(i), pb(i+1)
             write (6,*) 'e', d0(i-1), d0(i), d0(i+1)
             write (6,*) 'f', u0(i-1), u0(i), u0(i+1)
             write (6,*) 'g', p0(i-1), p0(i), p0(i+1)
             write (6,*) 'h', ubar(i-1,j), ubar(i,j), ubar(i+1,j)
             write (6,*) 'i', pbar(i-1,j), pbar(i,j), pbar(i+1,j)
             write (6,*) 'j', dls(i-1,j), dls(i,j), dls(i+1,j)
             write (6,*) 'k', drs(i-1,j), drs(i,j), drs(i+1,j)
             write (6,*) 'l', uls(i-1,j), uls(i,j), uls(i+1,j)
             write (6,*) 'm', urs(i-1,j), urs(i,j), urs(i+1,j)
             write (6,*) 'n', pls(i-1,j), pls(i,j), pls(i+1,j)
             write (6,*) 'o', prs(i-1,j), prs(i,j), prs(i+1,j)
             write (6,*) 'p', dslice(i-1,j), dslice(i,j), dslice(i+1,j)
             write (6,*) 'q', uslice(i-1,j), uslice(i,j), uslice(i+1,j)
             write (6,*) 'r', geslice(i-1,j), geslice(i,j), 
     $            geslice(i+1,j)
             write (6,*) 's', eslice(i-1,j), eslice(i,j), eslice(i+1,j)
             write (6,*) 't', geslice(i-1,j)*dslice(i-1,j)*(gamma-1.0), 
     &                   geslice(i  ,j)*dslice(i  ,j)*(gamma-1.0), 
     &                   geslice(i+1,j)*dslice(i+1,j)*(gamma-1.0)
             if (gravity .eq. 1) 
     &        write (6,*) grslice(i-1,j), grslice(i,j), grslice(i+1,j)
             write(0,*) 'stop in euler with e < 0'
             ERROR_MESSAGE
c
          endif
c           A small cheat: if idual is on, assume a cosmo sim, and institute
c              a minimum density.  This should be a parameter and pased in.
          if (idual .eq. 1 .and. gravity .eq. 1) 
     &         dnu(i) = max(dnu(i),1.0e-3)
          uold(i)     = uslice(i,j)
c          write(20,*) i, uslice(i,j), (uslice(i,j)*dslice(i,j) +
c     &         qc*(duub(i)-duub(i+1) +  pb(i) -  pb(i+1)))/dnu(i),
c     &               dnu(i), dslice(i,j)
          uslice(i,j) = (uslice(i,j)*dslice(i,j) +
     &         qc*(duub(i)-duub(i+1) +  pb(i) -  pb(i+1)))/dnu(i)
          vslice(i,j) = (vslice(i,j)*dslice(i,j) +
     &         qc*(duvb(i)-duvb(i+1)                    ))/dnu(i)
          wslice(i,j) = (wslice(i,j)*dslice(i,j) +
     &         qc*(duwb(i)-duwb(i+1)                    ))/dnu(i)
c          if ((eslice(i,j)*dslice(i,j) +
c     &         qc*(dueb(i)-dueb(i+1) + upb(i) - upb(i+1)))/dnu(i) .lt.
c     &         0.5*eslice(i,j)) then
c            write(6,*) i,j,qc*(dueb(i)-dueb(i+1)),
c     &                     qc*(upb(i)-upb(i+1)), eslice(i,j)*dslice(i,j)
c            write(6,*) dueb(i),dueb(i+1),db(i),db(i+1)
c            write(6,*) eb(i),eb(i+1),ub(i),ub(i+1)
c          endif
          eold = eslice(i,j)
          eslice(i,j) = max(0.1*eslice(i,j),
     &                  (eslice(i,j)*dslice(i,j) +
     &         qc*(dueb(i)-dueb(i+1) + upb(i) - upb(i+1)))/dnu(i) )
       enddo
c
c      Colour variables (note: colf already multiplied by dt)
c
       do n=1,ncolor
          do i=i1, i2
c             colslice(i,j,n) = (colslice(i,j,n)*dslice(i,j) +
c     &         (colf(i,j,n)-colf(i+1,j,n))/dx(i)      )/dnu(i) ! c*d conserved
             if (colslice(i,j,n)              +
     &         (colf(i,j,n)-colf(i+1,j,n))/dx(i) .lt. 0.0) write(6,*)
     &         'euler_c:',i,j,n,dx(i),
     &          colf(i-1,j,n),colf(i,j,n),colf(i+1,j,n),
     &          colslice(i-1,j,n),colslice(i,j,n),colslice(i+1,j,n),
     &          colb(i-1,n),colb(i,n),colb(i+1,n),
     &          colls(i-1,j,n),colls(i,j,n),colls(i+1,j,n),
     &          colrs(i-1,j,n),colrs(i,j,n),colrs(i+1,j,n),
     &          ub(i-1),ub(i),ub(i+1),
     &          u0(i-1),u0(i),u0(i+1),
     &          ubar(i-1,j),ubar(i,j),ubar(i+1,j)
             colslice(i,j,n) =  colslice(i,j,n)              +
     &         (colf(i,j,n)-colf(i+1,j,n))/dx(i)                ! c conserved
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
          do i=i1-1, i2+1
             eratio(i) = (eslice(i,j) - 
     &         0.5*(uslice(i,j)**2 + vslice(i,j)**2 + wslice(i,j)**2)
     &                    )/eslice(i,j)
          enddo
          do i=i1, i2
c             if (max(eratio(i-1),eratio(i),eratio(i+1)) .gt. eta2) then
                pcent(i) = max((gamma-1.0)*geslice(i,j)*dslice(i,j), 
     &                         tiny)
c             else
c                pcent(i) = 0.5*(pb(i)+pb(i+1))
c             endif
          enddo
          do i=i1, i2
             if (geslice(i,j) .lt. 0.0) write(6,*) i,j,geslice(i,j)
             qc = dt/dx(i)
             geslice(i,j) = max((geslice(i,j)*dslice(i,j) +
     &            qc*(dugeb(i)-dugeb(i+1) + pcent(i)*(ub(i)-ub(i+1))  
     &                                                       ))/dnu(i)
     &                          ,0.5*geslice(i,j))
c             if (geslice(i,j) .lt. tiny) geslice(i,j) = tiny
             if (geslice(i,j) .lt. 0.0) then
                write(6,*) i,j,dslice(i,j),dnu(i),pcent(i),eslice(i,j),
     &               qc*(dugeb(i)-dugeb(i+1)),
     &               qc*(pcent(i)*(ub(i)-ub(i+1))),
     &               geslice(i,j)*dslice(i,j),qc,
     &               uslice(i,j),uslice(i+1,j)
                write(6,*) pb(i),pb(i+1),p0(i),p0(i+1),
     &                     pbar(i,j),pbar(i+1,j),lbar(i),l0(i),
     &                     pls(i,j),prs(i,j),pls(i+1,j),prs(i+1,j)
                write(6,*) ub(i),ub(i+1),geb(i),geb(i+1),
     &                     db(i),db(i+1),dub(i),dub(i+1)
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
             dddx(i) = 2.0*( dslice(i,j) - max(dslice(i,j) -
     &                   0.5*dddx(i), min(dslice(i,j), dslice(i-1,j))))
             dddx(i) = 2.0*(-dslice(i,j) + max(dslice(i,j) +
     &                   0.5*dddx(i), min(dslice(i,j), dslice(i+1,j))))
          enddo
c
          do i=i1, i2
             grslice(i,j) = grslice(i,j) + 
     &                      0.5*dadx(i)*dddx(i)/(12.0*dslice(i,j))
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
             uslice(i,j) = uslice(i,j) + 
     &                     dt*grslice(i,j)*0.5*(dslice(i,j)/dnu(i)+1.0)
             eslice(i,j) = eslice(i,j) + dt*grslice(i,j)*0.5*
     &                     (uslice(i,j) + uold(i)*dslice(i,j)/dnu(i))
             if (eslice(i,j) .le. 0) write(6,*) 'eu1',i,j,eslice(i,j),
     &                     dslice(i,j),dnu(i),uslice(i,j),grslice(i,j),
     &                     dt*grslice(i,j)*0.5*(dslice(i,j)/dnu(i)+1.0),
     &                     dt*grslice(i,j)*0.5*
     &                     (uslice(i,j) + uold(i)*dslice(i,j)/dnu(i))
             eslice(i,j) = max(eslice(i,j), tiny)
#endif /* GRAVITY_METHOD1 */

#ifdef GRAVITY_METHOD2 /* (wrong) */
             uslice(i,j) = uslice(i,j) + dt*grslice(i,j)
             eslice(i,j) = eslice(i,j) + dt*grslice(i,j)*0.5*
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

