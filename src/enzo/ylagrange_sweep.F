#include "fortran.def"

c=======================================================================
c////////////////////  SUBROUTINE YLAGRANGE_SWEEP  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine ylagrange_sweep(i, d, e, u, v, w, ge, in, jn, kn,
     &                        gravity, gr_acc, idual, eta1, eta2,
     &                        is, ie, js, je, ks, ke,
     &                        gamma, pmin, dt, dx, dy, dz,
     &                        idiff, iflatten, isteepen, ipresfree,
     &                        nsubgrids, lface, rface,
     &                        fistart, fiend, fjstart, fjend,
     &                        dindex, eindex, geindex, 
     &                        uindex, vindex, windex, array,
     &                        ncolor, colorpt, coloff, colindex,
     &                        dls, drs, flatten, pbar, 
     &                        pls, prs, pslice, ubar,
     &                        uls, urs, vls, vrs,
     &                        wls, wrs, diffcoef, dslice,
     &                        eslice, uslice, vslice, wslice,
     &                        df, ef, uf, vf, 
     &                        wf, grslice, geslice, gef,
     &                        dynslice, dylslice,
     &                        colslice, colf, colls, colrs
     &                       )
c
c  CONTROL ROUTINE FOR Y-SWEEP LAGRANGE-REMAP VERSION OF PPM
c
c  written by: Greg Bryan
c  date:       June, 1994
c  modified1:
c
c  PURPOSE:  This routine servers as a wrapper for the eulerian version
c            of PPM that works on a two dimensional slice.  We extract
c            a slice, call INTEULR, R_SOLVER and then EULER.  Note
c            that this provides a natural way to do coarse-grain
c            parallelization on three dimension problems.
c
c  INPUTS:
c    d      - density field
c    dt     - timestep
c    dx,dy,dz - grid spacing
c    e      - total specific energy field
c    eta1   - (dual) selection parameter for gas energy (typically ~0.001)
c    eta2   - (dual) selection parameter for total energy (typically ~0.1)
c    gamma  - ideal gas constant
c    ge     - gas specific energy field (used when idual = 1)
c    gravity - gravity flag (0 = off)
c    gr_acc - acceleration due to gravity in this direction
c    i      - current slice position in x-direction
c    idiff  - diffusion flag (0 = off)
c    idual  - dual energy formalism flag (0 = off)
c    ie,je,ke - field active zone end index
c    iflatten - flattening flag (0 = off)
c    in,jn,kn - field dimensions
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    is,js,ks - field active zone start index
c    isteepen - steepening flag (0 = off)
c    pmin   - minimum pressure
c    u      - x-velocity field
c    v      - y-velocity field
c    w      - z-velocity field
c
c  LOCALS: (passed as temporaries)
c    diffcoef - diffusion coefficient in slice i
c    flatten - ammount of flattening (calculated in calcdiss)
c    dl,rs  - density at left and right edges of each cell
c    dslice - extracted 2d slice of the density   , d
c    el,rs  - total specific energy at left and right edges of each cell
c    eslice - extracted 2d slice of the energy    , e
c    geslice - extracted 2d slice of the gas energy, ge
c    pbar   - the pressure at the (left) cell interface 
c             after applying the Riemann solver
c    pl,rs  - pressure at left and right edges of each cell
c    pslice - extracted 2d slice of the pressure  , p
c    ubar   - the (1,2,3) velocity at the (left) cell interface
c             after applying the Riemann solver
c    ul,rs  - 1-velocity at left and right edges of each cell
c    uslice - extracted 2d slice of the 1-velocity, u
c    vl,rs  - 2-velocity at left and right edges of each cell
c    vslice - extracted 2d slice of the 2-velocity, v
c    wl,rs  - 3-velocity at left and right edges of each cell
c    wslice - extracted 2d slice of the 3-velocity, w
c
c  EXTERNALS:
c    pgas2d - computes pressure from equation of state (on a slice)
c    inteuler - computes the Eulerian left and right states for a slice
c    R_SOLVER - Riemann solver (Lagrangean)
c    euler  - converts the lagrangean Riemann results to eulerian
c             coordinates and then computes the Eulerian fluxes for a slice
c    calcdiss - Calculate dissiptation and flattening coefficients
c
c-----------------------------------------------------------------------
      implicit NONE
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer gravity, idiff, idual, iflatten, ipresfree, isteepen, 
     &        i, in, jn, kn, is, ie, js, je, ks, ke, nsubgrids,
     &        ncolor, coloff(ncolor)
      real    dt, eta1, eta2, gamma, pmin, vgy
      real    d(in,jn,kn), e(in,jn,kn), u(in,jn,kn), v(in,jn,kn), 
     &        w(in,jn,kn),ge(in,jn,kn), gr_acc(in,jn,kn),
     &        dx(in), dy(jn), dz(kn)
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      real    array(1), colorpt(1)
c
c  define local slices (passed as temps)
c
      real         dls(jn,kn),        drs(jn,kn),    flatten(jn,kn),
     &            pbar(jn,kn),        pls(jn,kn),
     &             prs(jn,kn),     pslice(jn,kn),       ubar(jn,kn),
     &             uls(jn,kn),        urs(jn,kn),        vls(jn,kn),
     &             vrs(jn,kn),        wls(jn,kn),        wrs(jn,kn),
     &        diffcoef(jn,kn),         df(jn,kn),         ef(jn,kn),
     &              uf(jn,kn),         vf(jn,kn),         wf(jn,kn),
     &             gef(jn,kn)
      real      dslice(jn,kn),     eslice(jn,kn),    grslice(jn,kn),
     &          uslice(jn,kn),     vslice(jn,kn),     wslice(jn,kn),
     &         geslice(jn,kn),   dynslice(jn,kn),   dylslice(jn,kn)
      real colslice(jn,kn,ncolor),  colf(jn,kn,ncolor),
     &        colls(jn,kn,ncolor), colrs(jn,kn,ncolor)
c
c  parameters
c
      integer ysweep
      parameter (ysweep = 2)
c
c  locals
c
      integer ic, idim, j, k, n, nxz, nyz, nzz, offset
      integer k1, k2
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Compute number of active zones
c
      nxz = ie - is + 1
      nyz = je - js + 1
      nzz = ke - ks + 1

      k1 = 1
      k2 = kn

!     k1 = ks
!     k2 = ke

c
      vgy = 0.0
c
c Copy from field to slice
c    (note that we have permuted uvw here and not in the routine calls)
c
      do k=k1, k2
         do j=1, jn
            dslice(j,k) = d(i,j,k)
            eslice(j,k) = e(i,j,k)
            uslice(j,k) = v(i,j,k)
            vslice(j,k) = w(i,j,k)
            wslice(j,k) = u(i,j,k)
         enddo
         if (gravity .eq. 1) then
            do j=1, jn
               grslice(j,k) = gr_acc(i,j,k)
            enddo
         endif
         if (idual .eq. 1) then
            do j=1, jn
               geslice(j,k) = ge(i,j,k)
            enddo
         endif
         do ic=1, ncolor
            do j=1, jn
               colslice(j,k,ic) = 
     &              colorpt(coloff(ic)+((k-1)*jn+j-1)*in+i)
            enddo
         enddo
      enddo
c
c  Compute the pressure on a slice
c
      if (idual .eq. 1) then
         call pgas2d_dual(dslice, eslice, geslice, pslice, 
     &                    uslice, vslice, wslice, eta1, eta2,
     &                    jn, kn, js-3, je+3, k1, k2, gamma, pmin)
      else
         call pgas2d(dslice, eslice, pslice, uslice, vslice, wslice,
     &               jn, kn, js-3, je+3, k1, k2, gamma, pmin)
      endif
c
c  If requested, compute diffusion and slope flattening coefficients
c
      if (idiff .ne. 0 .or. iflatten .ne. 0)
     &   call calcdiss(
     &            dslice, eslice, uslice, w, u, pslice, dy, dz, dx,
     &            jn, kn, in, js, je, ks, ke, i, nxz, ysweep,
     &            in, jn, kn, dt, gamma, idiff, iflatten,
     &            diffcoef, flatten
     &                )
c
c  Compute Lagrangean left and right states at zone edges via interpolation
c
      call intlgrg(
     &            dslice, pslice, uslice, dy, flatten,
     &            jn, kn, js, je, k1, k2,
     &            isteepen, iflatten, dt, gamma,
     &            dls, drs, pls, prs, uls, urs
     &             )
c
c  Compute (Lagrangian part of the) Riemann problem at each zone boundary
c
      call twoshock(
     &            dls, drs, pls, prs, uls, urs, 
     &            jn, kn, js, je+1, k1, k2, dt, gamma, pmin, ipresfree,
     &            pbar, ubar, gravity, grslice, idual, eta1
     &             )
c
c  Compute Lagrangean equations and update zones-centered quantities
c
      call lgrg    (
     &            dslice, eslice, geslice, uslice, dy, diffcoef,
     &            jn, kn, js, je, k1, k2, dt, gamma, idiff,
     &            idual, eta2, vslice, wslice,
     &            dls, drs, pls, prs, uls, urs, pbar, ubar,
     &            dylslice, dynslice, gravity, grslice
     &             )
c
c  Compute the pressure on a slice
c
      if (idual .eq. 1) then
         call pgas2d_dual(dslice, eslice, geslice, pslice, 
     &                    uslice, vslice, wslice, eta1, eta2,
     &                    jn, kn, js-3, je+3, k1, k2, gamma, pmin)
      else
         call pgas2d(dslice, eslice, pslice, uslice, vslice, wslice,
     &               jn, kn, js-3, je+3, k1, k2, gamma, pmin)
      endif
c
c  Compute remap fluxes
c
      call intrmp(
     &          dslice, eslice, pslice, uslice, vslice, wslice,
     &          dynslice, dylslice, geslice,
     &          jn, kn, js, je, k1, k2, i,
     &          isteepen, dt, gamma, vgy, idual,
     &          ysweep, in, kn, 
c nijb, nojb, eijb, eojb,
     &          df, ef, uf, vf, wf, gef
     &           )
c
c  Remap with fluxes
c
      call remap (
     &          dslice, eslice, uslice, vslice, wslice, 
     &          geslice, dynslice, dy,
     &          jn, kn, js, je, k1, k2, idual,
     &          df, ef, uf, vf, wf, gef
     &           )
c
c  Check this slice against the list of subgrids 
c     (all subgrid quantities are zero based)
c    Note that uf/vf/wf are switched to match u/v/wslice.
c
      do n=0, nsubgrids-1
        if (i .ge. fistart(n*3+2)+1 .and. i .le. fiend(n*3+2)+1) then
          idim = fiend(n*3+2) - fistart(n*3+2) + 1
          do k=fjstart(n*3+2)+1, fjend(n*3+2)+1
             offset = i-fistart(n*3+2) + (k-fjstart(n*3+2)-1)*idim
             array(dindex(n*6+3)+offset) = df(lface(n*3+2)+1, k)
             array(dindex(n*6+4)+offset) = df(rface(n*3+2)+2, k)
             array(eindex(n*6+3)+offset) = ef(lface(n*3+2)+1, k)
             array(eindex(n*6+4)+offset) = ef(rface(n*3+2)+2, k)
             if (nxz .gt. 1) then
                array(uindex(n*6+3)+offset) = wf(lface(n*3+2)+1, k)
                array(uindex(n*6+4)+offset) = wf(rface(n*3+2)+2, k)
             endif
             array(vindex(n*6+3)+offset) = uf(lface(n*3+2)+1, k)
             array(vindex(n*6+4)+offset) = uf(rface(n*3+2)+2, k)
             if (nzz .gt. 1) then
                array(windex(n*6+3)+offset) = vf(lface(n*3+2)+1, k)
                array(windex(n*6+4)+offset) = vf(rface(n*3+2)+2, k)
             endif
             if (idual .eq. 1) then
                array(geindex(n*6+3)+offset) = gef(lface(n*3+2)+1, k)
                array(geindex(n*6+4)+offset) = gef(rface(n*3+2)+2, k)
             endif
             do ic=1, ncolor
                array(colindex(n*6+3,ic)+offset) = 
     &                                       colf(lface(n*3+2)+1, k, ic)
                array(colindex(n*6+4,ic)+offset) = 
     &                                       colf(rface(n*3+2)+2, k, ic)
             enddo
          enddo
        endif
      enddo
c
c Copy from slice to field
c
      do k=k1, k2
         do j=1, jn
            d(i,j,k) = dslice(j,k)
            e(i,j,k) = eslice(j,k)
            v(i,j,k) = uslice(j,k)
            w(i,j,k) = vslice(j,k)
            u(i,j,k) = wslice(j,k)
         enddo
         if (idual .eq. 1) then
            do j=1, jn
               ge(i,j,k) = geslice(j,k)
            enddo
         endif
         do ic=1, ncolor
            do j=1, jn
               colorpt(coloff(ic)+((k-1)*jn+j-1)*in+i) = 
     &               colslice(j,k,ic)
            enddo
         enddo
      enddo
c
      return
      end
