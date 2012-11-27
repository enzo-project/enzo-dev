#include "fortran.def"
#include "error.def"
c=======================================================================
c////////////////////////  SUBROUTINE FLUX_HLL \\\\\\\\\\\\\\\\\\\\\\\\

      subroutine flux_hll(
     &     dslice, eslice, geslice, uslice, vslice, wslice, dx,diffcoef,
     &     idim, jdim, i1, i2, j1, j2, dt, gamma,
     &     idiff, idual, eta1, ifallback,
     &     dls, drs, pls, prs, uls, urs, vls, vrs, wls, wrs, gels, gers,
     &     df, uf, vf, wf, ef, gef, ges,
     &     ncolor, colslice, colls, colrs, colf
     &     )

c  RIEMANN SOLVER USING HLL APPROXIMATION
c
c  written by: John Wise
c  date:       July, 2010
c  modified1:  
c
c  PURPOSE:
c     This routine performs an approximate Riemann solution for each
c     zone edge described by dls, drs (etc) where dls(i) is the density
c     on the left side of zone interface (i), and drs(i) is the density
c     on the right side of zone interface (i).  Note that this differs
c     from the zone-centered way of indexing things.  See Toro's book 
c     (Ch. 10) for more details.
c
c  INPUT:
c    dt     - timestep in problem time
c    eta1   - (dual) selection parameter for total energy (typically ~0.001)
c    gamma  - ideal gas law constant
c    gravity - gravity flag (0 = off, 1 = on)
c    grslice - acceleration in this direction in this slice
c    i1,i2  - starting and ending addresses for dimension 1
c    idim   - declared leading dimension of slices
c    idual  - dual energy formalism flag (0 = off)
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    j1,j2  - starting and ending addresses for dimension 2
c    jdim   - declared second dimension of slices
c    dl,rs  - density at left and right edges of each cell
c    pl,rs  - pressure at left and right edges of each cell
c    pmin   - minimum allowed pressure
c    ul,rs  - 1-velocity at left and right edges of each cell
c    
c  OUTPUT:
c    
c  EXTERNALS:
c
c  LOCALS:
c    tolerance - close iteration to tolerance
c
c  PARAMETERS:
c    numiter - max number of Newton iterations to perform [8]
c
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c
      INTG_PREC ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC i1, i2, j1, j2, idim, jdim, ncolor, idiff, idual, 
     $     ifallback
      R_PREC gamma, dt, eta1
      R_PREC dx(idim)
      R_PREC uls(idim,jdim), vls(idim,jdim), wls(idim,jdim),
     $     urs(idim,jdim), vrs(idim,jdim), wrs(idim,jdim),
     $     pls(idim,jdim), dls(idim,jdim), gels(idim,jdim),
     $     prs(idim,jdim), drs(idim,jdim), gers(idim,jdim),
     $     colls(idim,jdim,ncolor), colrs(idim,jdim,ncolor),
     $     diffcoef(idim,jdim), dslice(idim,jdim), 
     $     uslice(idim,jdim), vslice(idim,jdim), wslice(idim,jdim),
     $     eslice(idim,jdim), geslice(idim,jdim),
     $     colslice(idim,jdim,ncolor)
      R_PREC df(idim,jdim), ef(idim,jdim), gef(idim,jdim), 
     $     ges(idim,jdim), uf(idim,jdim), vf(idim,jdim), 
     $     wf(idim,jdim), colf(idim,jdim,ncolor)
c
c  local declarations
c
      INTG_PREC i, j, n
      R_PREC sqrtdl, sqrtdr, gamma1, isdlpdr, vroe1, vroe2, vroe3, v2,
     $     hroe, tl, tr, dl, dr, q1, dub, delu, energy, qc,
     $     gamma1i, csl0, csr0
      R_PREC cs(ijkn), char1(ijkn), char2(ijkn), csl(ijkn), csr(ijkn),
     $     cw(ijkn), cp(ijkn), diffd(ijkn), diffuu(ijkn), diffuv(ijkn),
     $     diffuw(ijkn), diffue(ijkn), diffuge(ijkn), 
     $     bp(ijkn), bm(ijkn), el(ijkn), er(ijkn), pcent(ijkn),
     $     bm0(ijkn), bp0(ijkn)
      R_PREC dfl(ijkn), dfr(ijkn), ufl(ijkn), ufr(ijkn), 
     $     vfl(ijkn), vfr(ijkn), wfl(ijkn), wfr(ijkn),
     $     efl(ijkn), efr(ijkn), gefl(ijkn), gefr(ijkn),
     $     gesl(ijkn), gesr(ijkn),
     $     colfl(ijkn,MAX_COLOR), colfr(ijkn,MAX_COLOR)
      R_PREC sl(ijkn), sr(ijkn), sm(ijkn)
      R_PREC dubl, dubr

      gamma1 = gamma - 1._RKIND
      gamma1i = 1._RKIND / gamma1
c
c  Loop through this slice, doing 1 sweep at a time
c
      do j = j1,j2
c
c     Compute the Roe-averaged data from the left and right states
c
      do i = i1,i2+1
         sqrtdl = sqrt(dls(i,j))
         sqrtdr = sqrt(drs(i,j))
         isdlpdr = 1._RKIND / (sqrtdl + sqrtdr)
         vroe1 = (sqrtdl * uls(i,j) + sqrtdr * urs(i,j)) * isdlpdr
         vroe2 = (sqrtdl * vls(i,j) + sqrtdr * vrs(i,j)) * isdlpdr
         vroe3 = (sqrtdl * wls(i,j) + sqrtdr * wrs(i,j)) * isdlpdr
         v2 = vroe1**2 + vroe2**2 + vroe3**2
c
c     The enthalpy H=(E+P)/d is averaged for adiabatic flows, rather
c     than E or P directly.  
c     sqrtql*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
c         
         el(i) = gamma1i * pls(i,j) + 0.5_RKIND*dls(i,j)*
     $        (uls(i,j)**2 + vls(i,j)**2 + wls(i,j)**2)
         er(i) = gamma1i * prs(i,j) + 0.5_RKIND*drs(i,j)*
     $        (urs(i,j)**2 + vrs(i,j)**2 + wrs(i,j)**2)
         hroe = ((el(i) + pls(i,j))/sqrtdl + 
     $        (er(i) + prs(i,j))/sqrtdr) * isdlpdr
c
c     Compute characteristics using Roe-averaged values
c
         cs(i) = sqrt(gamma1*max((hroe - 0.5_RKIND*v2), tiny))
         char1(i) = vroe1 - cs(i)
         char2(i) = vroe1 + cs(i)
      enddo
c
c     Compute the min/max wave speeds
c
      do i = i1,i2+1

         csl0 = sqrt(gamma*pls(i,j)/dls(i,j))
         csr0 = sqrt(gamma*prs(i,j)/drs(i,j))

         csl(i) = min(uls(i,j)-csl0, char1(i))
         csr(i) = max(urs(i,j)+csr0, char2(i))

         bm(i) = min(csl(i), 0._RKIND)
         bp(i) = max(csr(i), 0._RKIND)
         
         bm0(i) = uls(i,j) - bm(i)
         bp0(i) = urs(i,j) - bp(i)

      enddo
c
c     Weights for different (left, right, star) regions in the cell
c
      do i = i1,i2+1
         q1 = (bp(i) + bm(i)) / (bp(i) - bm(i))
         sl(i) = 0.5_RKIND * (1._RKIND + q1)
         sr(i) = 0.5_RKIND * (1._RKIND - q1)
      enddo
c
c     Compute diffusion terms
c
      if (idiff .ne. 0) then
         do i = i1,i2+1
            diffd(i) = diffcoef(i,j) * (dslice(i-1,j) - dslice(i,j))
            diffuu(i) = diffcoef(i,j) * 
     $           (dslice(i-1,j)*uslice(i-1,j) - dslice(i,j)*uslice(i,j))
            diffuv(i) = diffcoef(i,j) * 
     $           (dslice(i-1,j)*vslice(i-1,j) - dslice(i,j)*vslice(i,j))
            diffuw(i) = diffcoef(i,j) * 
     $           (dslice(i-1,j)*wslice(i-1,j) - dslice(i,j)*wslice(i,j))
            diffue(i) = diffcoef(i,j) * 
     $           (dslice(i-1,j)*eslice(i-1,j) - dslice(i,j)*eslice(i,j))
         enddo
         if (idual .eq. 1) then
            do i = i1,i2+1
               diffuge(i) = 
     $              diffcoef(i,j) * (dslice(i-1,j)*geslice(i-1,j) -
     $              dslice(i,j)*geslice(i,j))
            enddo
         endif
      else
         do i = i1,i2+1
            diffd(i) = 0._RKIND
            diffuu(i) = 0._RKIND
            diffuv(i) = 0._RKIND
            diffuw(i) = 0._RKIND
            diffue(i) = 0._RKIND
            diffuge(i) = 0._RKIND
         enddo
      endif
c
c     Compute the left and right fluxes along the characteristics.
c     Compute the color fluxes afterwards (loop n, then i) for cache
c     efficiency.  For dual energy formalism, do not include the source
c     term.  Treat like an advected quantity.
c
      do i = i1,i2+1

         dubl = dls(i,j) * uls(i,j)
         dubr = drs(i,j) * urs(i,j)

         dfl(i) = dls(i,j) * bm0(i)
         dfr(i) = drs(i,j) * bp0(i)

         ufl(i) = dubl * bm0(i) + pls(i,j)
         ufr(i) = dubr * bp0(i) + prs(i,j)

         vfl(i) = dls(i,j)*vls(i,j) * bm0(i)
         vfr(i) = drs(i,j)*vrs(i,j) * bp0(i)

         wfl(i) = dls(i,j)*wls(i,j) * bm0(i)
         wfr(i) = drs(i,j)*wrs(i,j) * bp0(i)

         efl(i) = el(i) * bm0(i) + pls(i,j)*uls(i,j)
         efr(i) = er(i) * bp0(i) + prs(i,j)*urs(i,j)

      enddo                     ! i
c
c     Gas energy is treated like an advected quantity
c
      if (idual .eq. 1) then
         do i = i1,i2+1
            gefl(i) = bm0(i) * gels(i,j) * dls(i,j)
            gefr(i) = bp0(i) * gers(i,j) * drs(i,j)
            gesl(i) = bm0(i)
            gesr(i) = bp0(i)
         enddo
      endif
c
c     Now for the colors
c
      do n=1,ncolor
         do i=i1,i2+1
            colfl(i,n) = bm0(i) * colls(i,j,n)
            colfr(i,n) = bp0(i) * colrs(i,j,n)
         enddo                  ! i
      enddo                     ! colors
c
c     Compute HLL flux at interface plus diffusion
c
      do i = i1,i2+1
         df(i,j) = sl(i)*dfl(i) + sr(i)*dfr(i) + diffd(i)
         uf(i,j) = sl(i)*ufl(i) + sr(i)*ufr(i) + diffuu(i)
         vf(i,j) = sl(i)*vfl(i) + sr(i)*vfr(i) + diffuv(i)
         wf(i,j) = sl(i)*wfl(i) + sr(i)*wfr(i) + diffuw(i)
         ef(i,j) = sl(i)*efl(i) + sr(i)*efr(i) + diffue(i)
      enddo

      if (idual .eq. 1) then
         do i = i1,i2+1
            gef(i,j) = sl(i)*gefl(i) + sr(i)*gefr(i) + diffuge(i)
         enddo
      endif

      do n = 1, ncolor
         do i = i1,i2+1
            colf(i,j,n) = sl(i)*colfl(i,n) + sr(i)*colfr(i,n)
         enddo
      enddo
c
c     If in the star region (sm<0), add flux from the intermediate waves
c
#ifdef UNUSED
      do i = i1,i2+1
         if (sm(i) .lt. 0._RKIND) then
            df(i,j) = df(i,j) + sm(i) * (drs(i,j) - dls(i,j))
            uf(i,j) = uf(i,j) + sm(i) * 
     $           (drs(i,j)*urs(i,j) - dls(i,j)*uls(i,j))
            vf(i,j) = vf(i,j) + sm(i) * 
     $           (drs(i,j)*vrs(i,j) - dls(i,j)*vls(i,j))
            wf(i,j) = wf(i,j) + sm(i) * 
     $           (drs(i,j)*wrs(i,j) - dls(i,j)*wls(i,j))
            ef(i,j) = ef(i,j) + sm(i) * (er(i) - el(i))
            if (idual .eq. 1) then
               gef(i,j) = gef(i,j) + sm(i) *
     $              (drs(i,j)*gers(i,j) - dls(i,j)*gels(i,j))
            endif
         endif
      enddo

      do n = 1, ncolor
         do i = i1,i2+1
            if (sm(i) .lt. 0._RKIND) then
               colf(i,j,n) = colf(i,j,n) + sm(i) *
     $              (colrs(i,j,n) - colls(i,j,n))
            endif
         enddo
      enddo
#endif
c
c     Calculate source term for gas energy
c
      if (idual .eq. 1) then
         do i=i1, i2+1
            pcent(i) = max((gamma-1._RKIND)*geslice(i,j)*dslice(i,j), 
     $           tiny)
         enddo
         do i=i1, i2+1
            qc = dt/dx(i)
            ges(i,j) = qc * pcent(i) * 
     $           (sl(i  )*gesl(i  ) + sr(i  )*gesr(i  ) -
     $            sl(i+1)*gesl(i+1) - sr(i+1)*gesr(i+1))
         enddo
      endif
c
c     Finally absorb dt/dx into fluxes.  Advected quantities (including
c     the gas energy) are only multiplied by dt.
c
      do i = i1,i2+1
         qc = dt/dx(i)
         df(i,j) = qc*df(i,j)
         uf(i,j) = qc*uf(i,j)
         vf(i,j) = qc*vf(i,j)
         wf(i,j) = qc*wf(i,j)
         ef(i,j) = qc*ef(i,j)
      enddo
      
      if (idual .eq. 1) then
         do i = i1,i2+1
            gef(i,j) = dt*gef(i,j)
         enddo
      endif

      do n = 1,ncolor
         do i = i1,i2+1
            colf(i,j,n) = dt*colf(i,j,n)
         enddo
      enddo
c
c     Check here for negative densities and energies
c
      do i=i1, i2
         if (dslice(i,j)+df(i,j)-df(i+1,j) .le. 0._RKIND) then
            if (eslice(i,j) .lt. 0._RKIND) then
               write (6,*) 'flux_hll: eslice < 0: ', i, j
               write (6,*) qc*(ef(i,j)-ef(i+1,j)),
     &              eslice(i,j)*dslice(i,j)
            else
               write (6,*) 'flux_hll: dnu <= 0:', i, j
            endif
            write (6,*) 'a', dslice(i,j) + df(i,j)-df(i+1,j),
     $           df(i,j)-df(i+1,j), dslice(i,j), dt/dx(i), dt
            write (6,*) 'b', csl(i-1), csl(i), csl(i+1)
            write (6,*) 'c', csr(i-1), csr(i), csr(i+1)
            write (6,*) 'd', char1(i-1), char1(i), char1(i+1)
            write (6,*) 'e', char2(i-1), char2(i), char2(i+1)
            write (6,*) 'f', sl(i-1), sl(i), sl(i+1)
            write (6,*) 'g', sr(i-1), sr(i), sr(i+1)
            write (6,*) 'h1', dfl(i-1), dfl(i), dfl(i+1)
            write (6,*) 'h2', dfr(i-1), dfr(i), dfr(i+1)
            write (6,*) 'i', df(i-1,j), df(i,j), df(i+1,j)
            write (6,*) 'j', uf(i-1,j), uf(i,j), uf(i+1,j)
            write (6,*) 'k', vf(i-1,j), vf(i,j), vf(i+1,j)
            write (6,*) 'l', wf(i-1,j), wf(i,j), wf(i+1,j)
            write (6,*) 'm', ef(i-1,j), ef(i,j), ef(i+1,j)
            write (6,*) 'n', gef(i-1,j), gef(i,j), gef(i+1,j)
            write (6,*) 'o', dls(i-1,j), dls(i,j), dls(i+1,j)
            write (6,*) 'p', drs(i-1,j), drs(i,j), drs(i+1,j)
            write (6,*) 'q', uls(i-1,j), uls(i,j), uls(i+1,j)
            write (6,*) 'r', urs(i-1,j), urs(i,j), urs(i+1,j)
            write (6,*) 's', pls(i-1,j), pls(i,j), pls(i+1,j)
            write (6,*) 't', prs(i-1,j), prs(i,j), prs(i+1,j)
            write (6,*) 'u', dslice(i-1,j), dslice(i,j), dslice(i+1,j)
            write (6,*) 'v', uslice(i-1,j), uslice(i,j), uslice(i+1,j)
            write (6,*) 'w', eslice(i-1,j), eslice(i,j), eslice(i+1,j)
            if (idual .gt. 0) then
               write (6,*) 'x', geslice(i-1,j), geslice(i,j), 
     $              geslice(i+1,j)
               write (6,*) 'y',geslice(i-1,j)*dslice(i-1,j)*(gamma-1.0), 
     &              geslice(i  ,j)*dslice(i  ,j)*(gamma-1.0), 
     &              geslice(i+1,j)*dslice(i+1,j)*(gamma-1.0)
            endif
c     if (gravity .eq. 1) 
c     &        write (6,*) grslice(i-1,j), grslice(i,j), grslice(i+1,j)
            if (ifallback .eq. 1) 
     $           write(0,*) 'Cannot fallback anymore.  Already in HLL.'
            write(0,*) 'stop in euler with e < 0'
            ERROR_MESSAGE
         endif
      enddo
c
c     Now check for negative colors (Warning only)
c
      do n=1, ncolor
         do i=i1, i2
            if (colslice(i,j,n) +
     $           (colf(i,j,n)-colf(i+1,j,n))/dx(i)
     $           .lt. 0._RKIND) then
               write(6,*) 'flux_hll: negative color',i,j,n
               write(6,*) 'a', dx(i), colf(i-1,j,n),colf(i,j,n),
     $              colf(i+1,j,n)
               write(6,*) 'b', colslice(i,j,n)+(colf(i,j,n)-
     $              colf(i+1,j,n))/dx(i),
     $              colslice(i,j,n),
     $              (colf(i,j,n)-colf(i+1,j,n))/dx(i)
               write(6,*) 'c', colslice(i-1,j,n),colslice(i,j,n),
     $              colslice(i+1,j,n)
               write(6,*) 'd', colls(i-1,j,n),colls(i,j,n),
     $              colls(i+1,j,n)
               write(6,*) 'e', colrs(i-1,j,n),colrs(i,j,n),
     $              colrs(i+1,j,n)
               write(6,*) 'f', uls(i-1,j),uls(i,j),uls(i+1,j)
               write(6,*) 'g', urs(i-1,j),urs(i,j),urs(i+1,j)
               WARNING_MESSAGE
            endif
         enddo
      enddo
      
      enddo                     ! j

      end
