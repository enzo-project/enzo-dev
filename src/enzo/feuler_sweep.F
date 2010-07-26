#include "fortran.def"
c=======================================================================
c/////////////////////  SUBROUTINE XEULER_SWEEP  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine feuler_sweep(iter, v, w,
     &                        in, jn, kn, d0n, d1n, d2n, d2nactive,
     &                        gravity, idual, eta1, eta2,
     &                        d0start, d0end, d1start, d1end,
     &                        gamma, pmin, dt, dx, dy, dz,
     &                        idiff, iflatten, isteepen, ipresfree,
     &                        ncolor, sweep_dir,
     &                        dslice, eslice, uslice, vslice, wslice,
     &                        grslice, geslice, colslice,
     &                        df, ef, uf, vf, wf, gef, colf
     &                       )
c
c  CONTROL ROUTINE FOR X-SWEEP DIRECT EULER VERSION OF PPM
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
c    dt     - timestep
c    dx,dy,dz - grid spacing
c    eta1   - (dual) selection parameter for gas energy (typically ~0.001)
c    eta2   - (dual) selection parameter for total energy (typically ~0.1)
c    gamma  - ideal gas constant
c    gravity - gravity flag (0 = off)
c    idiff  - diffusion flag (0 = off)
c    idual  - dual energy formalism flag (0 = off)
c    ie,je,ke - field active zone end index
c    iflatten - flattening flag (0 = off)
c    in,jn,kn - field dimensions
c    ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c    is,js,ks - field active zone start index
c    isteepen - steepening flag (0 = off)
c    k        - current slice position in direction 3
c    pmin   - minimum pressure
c    v      - y-velocity field
c    w      - z-velocity field
c
c    l,rface  - position of left and right flux surfaces (for each dim)
c    fi,jstart,end - start and end indexes for the flux surfaces (for each dim)
c    d...index - indexes for each field into the flux 'array'
c    colorpt - the color 'array' (contains all color fields)
c    ncolor  - the number of color fields
c    coloff  - offsets into the color array (for each color field)
c
c  OUPUTS:
c    array  - array of subgrid fluxes
c    colorpt - color fields
c
c  LOCALS:  (passed as temporaries in argument list)
c    diffcoef - diffusion coefficient in slice k
c    df     - density flux
c    flatten - ammount of flattening (calculated in calcdiss)
c    dl,rs  - density at left and right edges of each cell
c    dslice - extracted 2d slice of the density   , d
c    dx,dy,dz - grid dimension
c    ef     - total energy flux
c    el,rs  - total specific energy at left and right edges of each cell
c    eslice - extracted 2d slice of the energy    , e
c    geslice - extracted 2d slice of the gas energy, ge
c    pbar   - the pressure at the (left) cell interface 
c             after applying the Riemann solver
c    pl,rs  - pressure at left and right edges of each cell
c    pslice - extracted 2d slice of the pressure  , p
c    ubar   - the (1,2,3) velocity at the (left) cell interface
c             after applying the Riemann solver
c    uf     - 1-momuntum flux
c    ul,rs  - 1-velocity at left and right edges of each cell
c    uslice - extracted 2d slice of the 1-velocity, u
c    vf     - 2-momentum flux
c    vl,rs  - 2-velocity at left and right edges of each cell
c    vslice - extracted 2d slice of the 2-velocity, v
c    wf     - 3-momentum flux
c    wl,rs  - 3-velocity at left and right edges of each cell
c    wslice - extracted 2d slice of the 3-velocity, w
c
c    colf   - color flux
c    coll,rs - color at left and right edges of each cell
c    colslice - extracted 2d slice of the color
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
      integer gravity, idiff, idual, iflatten, isteepen, ipresfree,
     &        in, jn, kn, d0n, d1n, d2n, d0start, d0end, d1start, d1end,
     &        iter, ncolor, d2nactive
      real    dt, eta1, eta2, gamma, pmin
      real    v(in,jn,kn), w(in,jn,kn), 
     &        dx(d0n), dy(d1n), dz(d2n)
      real dslice(d0n,d1n), eslice(d0n,d1n), grslice(d0n,d1n),
     &     uslice(d0n,d1n), vslice(d0n,d1n),  wslice(d0n,d1n), 
     &    geslice(d0n,d1n), colslice(d0n,d1n,ncolor)
      real     df(d0n,d1n),     ef(d0n,d1n),
     &         uf(d0n,d1n),     vf(d0n,d1n),     wf(d0n,d1n),
     &        gef(d0n,d1n),   colf(d0n,d1n,ncolor)
      integer sweep_dir
c
c  locals
c
      integer d2, d0s, d1s, d0e, d1e
      integer i
c
c  define local slices (declared here as temps)
c
      real    dls(d0n,d1n),    drs(d0n,d1n),flatten(d0n,d1n),
     &       pbar(d0n,d1n),    pls(d0n,d1n),
     &        prs(d0n,d1n), pslice(d0n,d1n),   ubar(d0n,d1n),
     &        uls(d0n,d1n),    urs(d0n,d1n),    vls(d0n,d1n),
     &        vrs(d0n,d1n),    wls(d0n,d1n),    wrs(d0n,d1n),
     &   diffcoef(d0n,d1n),   gels(d0n,d1n),   gers(d0n,d1n),
     &   colls(d0n,d1n,ncolor), colrs(d0n,d1n,ncolor)

c
c  convert to base-1
c
      d2 = iter+1
      d0s = d0start+1
      d1s = d1start+1
      d0e = d0end+1
      d1e = d1end+1

c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Compute the pressure on a slice
c
      if (idual .eq. 1) then
         call pgas2d_dual(dslice, eslice, geslice, pslice,
     &                    uslice, vslice, wslice, eta1, eta2,
     &                    d0n, d1n, d0s-3, d0e+3, d1s, d1e, gamma, pmin)
      else
         call pgas2d(dslice, eslice, pslice, uslice, vslice, wslice,
     &               d0n, d1n, d0s-3, d0e+3, d1s, d1e, gamma, pmin)
      endif
c
c  If requested, compute diffusion and slope flattening coefficients
c
      if (idiff .ne. 0 .or. iflatten .ne. 0)
     &   call calcdiss(
     &            dslice, eslice, uslice, v, w, pslice, dx, dy, dz,
     &            d0n, d1n, d2n, d0s, d0e, d1s, d1e, d2, d2nactive, 
     &            sweep_dir, in, jn, kn, dt, gamma, idiff, iflatten,
     &            diffcoef, flatten
     &                )
c
c  Compute Eulerian left and right states at zone edges via interpolation
c
      call inteuler(
     &            dslice, pslice, gravity, grslice, geslice,
     &            uslice, vslice, wslice, dx, flatten,
     &            d0n, d1n, d0s, d0e, d1s, d1e, idual, eta1, eta2,
     &            isteepen, iflatten, dt, gamma, ipresfree,
     &            dls, drs, pls, prs, gels, gers,
     &            uls, urs, vls, vrs, wls, wrs,
     &            ncolor, colslice, colls, colrs
     &             )
c
c  Compute (Lagrangian part of the) Riemann problem at each zone boundary
c
      call twoshock(
     &            dls, drs, pls, prs, uls, urs, 
     &            d0n, d1n, d0s, d0e+1, d1s, d1e, 
     &            dt, gamma, pmin, ipresfree,
     &            pbar, ubar, gravity, grslice, idual, eta1
     &             )
c
c  Compute Eulerian fluxes and update zones-centered quantities
c
      call euler   (
     &            dslice, eslice, grslice, geslice,
     &            uslice, vslice, wslice, dx, diffcoef,
     &            d0n, d1n, d0s, d0e, d1s, d1e, dt, 
     &            gamma, idiff, gravity, idual, eta1, eta2,
     &            dls, drs, pls, prs, gels, gers,
     &            uls, urs, vls, vrs, wls, wrs,
     &            pbar, ubar,
     &            df, ef, uf, vf, wf, gef,
     &            ncolor, colslice, colls, colrs, colf
     &             )
c
c  If necessary, recompute the pressure to correctly set ge and e
c
      if (idual .eq. 1)
     &    call pgas2d_dual(dslice, eslice, geslice, pslice,
     &                    uslice, vslice, wslice, eta1, eta2,
     &                    d0n, d1n, d0s-3, d0e+3, d1s, d1e, gamma, pmin)
c
c      write(*,*) 'd(4,4,4) = ',d(4,4,4)
c      do j=js,je
c        do i=1,in
c          write(*,*) 'ds(',i,j,') = ',dslice(i,j)
c          write(*,*) 'd(',i,j,k,') = ',d(i,j,k)
c          write(*,*) 'd(',i,') = ',dslice(i,1),uslice(i,1),dt
c        enddo
c      enddo
c      write(*,*) 'colslice(4,4,1) = ',colslice(4,4,1)
c
      return
      end
