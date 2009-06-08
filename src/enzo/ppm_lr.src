#include "fortran.def"
#include "error.def"

c=======================================================================
c/////////////////////////  SUBROUTINE PPM_LR  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine ppm_lr(d, e, u, v, w, ge,
     &                  gravity, gr_xacc, gr_yacc, gr_zacc,
     &                  gamma, dt, nhy, dx, dy, dz,
     &                  rank, in, jn, kn, start, end,
     &                  gridvel, iflatten, ipresfree,
     &                  idiff, isteepen, idual, eta1, eta2,
     &                  nsubgrids, lface, rface,
     &                  fistart, fiend, fjstart, fjend,
     &                  array, dindex, eindex,
     &                  uindex, vindex, windex, geindex, tmp,
     &                  ncolor, colorpt, coloff, colindex)
c
c  PERFORMS PPM (LAGRANGE REMAP) UPDATE FOR ONE TIMESTEP (WRAPPER)
c
c  NOTE: PPMLR doesnt really work in this amr, because of boundary values
c        problems.  The full algorithm needs at least six ghost zones, but
c        generally we provide only three.
c
c  written by: Greg Bryan
c       (some of this is based on code originally written by Jim Stone)
c  date:       June, 1993
c  modified1:  November, 1994; by GB
c              Modified to work as an entirely callable routine within
c              an AMR packaged.  There is now only one exterior parameter:
c              MAX_ANY_SINGLE_DIRECTION.
c
c  PURPOSE:  Advances the fluid equations by one timestep using the
c    PPM Direct Eulerian algorithm in 3-D.  This routine is a callable
c    wrapper for the 1d sweeps.  It assumes (in classic fortran style)
c    that the arguments are defined to dimension three; however, it does
c    the right thing if called with one or more collapsed dimensions
c    (i.e. with the width set to 1).
c
c  EXTERNALS:
c    x,y,zeuler_sweep - routines to compute the Eulerian step in
c                       one dimension
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     end     - array (of dimension 3) specifying the start of the active
c               region for reach dimension (zero based)
c     eta1    - (dual) selection parameter for gas energy (typically ~0.1)
c     eta2    - (dual) selection parameter for total energy (typically ~0.001)
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     gridvel - bulk grid velocity (vector of length 3)
c     i,j,kn  - dimensions of field arrays
c     idiff   - diffusion flag (0 = off)
c     idual   - dual energy formalism flag (0 = off)
c     ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c     nhy     - cycle number (for better operator splitting)
c     rank    - dimension of problem (not currently used)
c     start   - array (of dimension 3) specifying the start of the active
c               region fo reach dimension (zero based)
c     tmp     - temporary work space (30 * largest_slice)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c     ifloat_size - number of bytes in float (4 or 8)
c
c  LOCALS:
c
c-----------------------------------------------------------------------
      implicit NONE
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer gravity, idiff, idual, iflatten, isteepen, nhy, rank,
     &        ifloat_size
      integer ipresfree, end(3), in, jn, kn, nsubgrids, start(3),
     &        ncolor, coloff(ncolor)
      integer fistart(nsubgrids*3), fiend(nsubgrids*3),
     &        fjstart(nsubgrids*3), fjend(nsubgrids*3), 
     &        lface(nsubgrids*3), rface(nsubgrids*3)
      integer dindex(nsubgrids*6), eindex(nsubgrids*6),
     &        uindex(nsubgrids*6), vindex(nsubgrids*6),
     &        windex(nsubgrids*6),geindex(nsubgrids*6),
     &        colindex(nsubgrids*6,ncolor)
      real d(in,jn,kn), e(in,jn,kn), u(in,jn,kn),
     &     v(in,jn,kn), w(in,jn,kn),ge(in,jn,jn),
     &     gr_xacc(in,jn,kn), gr_yacc(in,jn,kn), gr_zacc(in,jn,kn),
     &     dx(in), dy(jn), dz(kn)
      real dt, eta1, eta2, gamma, gridvel(3), pmin
      real array(1), colorpt(1)
c
c  Parameters
c
c
c  Locals
c
      integer i, ie, is, ixyz, j, je, js, k, ke, ks,
     &        n, nxz, nyz, nzz, ms
      integer i1, i2, j1, j2, k1, k2
      real tmp(1)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Error check
c
      if (max(in,jn,kn) .gt. MAX_ANY_SINGLE_DIRECTION) then
         write(6,*) 'PPM_LR: A grid dimension is too long.'
         write(6,*) '   (increase max_any_single_direction.)'
         ERROR_MESSAGE
      endif
c
      if (ipresfree .eq. 1) then
         write(6,*) 'PPM_LR: does not support pressure-free'
         ERROR_MESSAGE
      endif
c
c  Convert arguments to usable form
c
      is = start(1) + 1
      js = start(2) + 1
      ks = start(3) + 1
      ie = end(1) + 1
      je = end(2) + 1
      ke = end(3) + 1
      ms = max(in*jn, jn*kn, kn*in)
      nxz = ie-is+1
      nyz = je-js+1
      nzz = ke-ks+1

      i1 = 1
      i2 = in
      j1 = 1
      j2 = jn
      k1 = 1
      k2 = kn

!     i1 = is
!     i2 = ie
!     j1 = js
!     j2 = je
!     k1 = ks
!     k2 = ke

c
c  Set minimum pressure (better if it were a parameter)
c
      pmin = tiny
c
c  Loop over directions, using a Stange-type splitting
c
      ixyz = mod(nhy,3)
      do n=ixyz+2,ixyz,-1
c
c  Update in x-direction
c
        if (mod(n,3) .eq. 0 .and. nxz .gt. 1) then
c
!         write (6,*) 'SOLVER: starting x-sweep'
c
!         do k=ks, ke
          do k=k1, k2
            call xlagrange_sweep(k, d, e, u, v, w, ge, in, jn, kn,
     &        gravity, gr_xacc, idual, eta1, eta2,
     &        is, ie, js, je, ks, ke,
     &        gamma, pmin, dt, dx, dy, dz,
     &        idiff, iflatten, isteepen, ipresfree,
     &        nsubgrids, lface, rface,
     &        fistart, fiend, fjstart, fjend,
     &        dindex, eindex, geindex,
     &        uindex, vindex, windex, array,
     &        ncolor, colorpt, coloff, colindex,
     &        tmp(1+ms*0),  tmp(1+ms*1),  tmp(1+ms*2),  tmp(1+ms*3),
     &        tmp(1+ms*4),  tmp(1+ms*5),  tmp(1+ms*6),  tmp(1+ms*7), 
     &        tmp(1+ms*8),  tmp(1+ms*9),  tmp(1+ms*10), tmp(1+ms*11),
     &        tmp(1+ms*12), tmp(1+ms*13), tmp(1+ms*14), tmp(1+ms*15),
     &        tmp(1+ms*16), tmp(1+ms*17), tmp(1+ms*18), tmp(1+ms*19),
     &        tmp(1+ms*20), tmp(1+ms*21), tmp(1+ms*22), tmp(1+ms*23),
     &        tmp(1+ms*24), tmp(1+ms*25), tmp(1+ms*26), tmp(1+ms*27),
     &        tmp(1+ms*28), tmp(1+ms*29),
     &        tmp(1+ms*(30+0*ncolor)), tmp(1+ms*(30+1*ncolor)), 
     &        tmp(1+ms*(30+2*ncolor)), tmp(1+ms*(30+3*ncolor))
     &                       )
          enddo
c
        endif
c
c  Update in y-direction
c
        if (mod(n,3) .eq. 1 .and. nyz .gt. 1) then
c
!         write (6,*) 'SOLVER: starting y-sweep'
c
!         do i=is, ie
          do i=i1, i2
            call ylagrange_sweep(i, d, e, u, v, w, ge, in, jn, kn,
     &        gravity, gr_yacc, idual, eta1, eta2,
     &        is, ie, js, je, ks, ke,
     &        gamma, pmin, dt, dx, dy, dz,
     &        idiff, iflatten, isteepen, ipresfree,
     &        nsubgrids, lface, rface,
     &        fistart, fiend, fjstart, fjend,
     &        dindex, eindex, geindex,
     &        uindex, vindex, windex, array,
     &        ncolor, colorpt, coloff, colindex,
     &        tmp(1+ms*0),  tmp(1+ms*1),  tmp(1+ms*2),  tmp(1+ms*3),
     &        tmp(1+ms*4),  tmp(1+ms*5),  tmp(1+ms*6),  tmp(1+ms*7), 
     &        tmp(1+ms*8),  tmp(1+ms*9),  tmp(1+ms*10), tmp(1+ms*11),
     &        tmp(1+ms*12), tmp(1+ms*13), tmp(1+ms*14), tmp(1+ms*15),
     &        tmp(1+ms*16), tmp(1+ms*17), tmp(1+ms*18), tmp(1+ms*19),
     &        tmp(1+ms*20), tmp(1+ms*21), tmp(1+ms*22), tmp(1+ms*23),
     &        tmp(1+ms*24), tmp(1+ms*25), tmp(1+ms*26), tmp(1+ms*27),
     &        tmp(1+ms*28), tmp(1+ms*29),
     &        tmp(1+ms*(30+0*ncolor)), tmp(1+ms*(30+1*ncolor)), 
     &        tmp(1+ms*(30+2*ncolor)), tmp(1+ms*(30+3*ncolor))
     &                       )
          enddo
c
        endif
c
c  Update in z-direction
c
        if (mod(n,3) .eq. 2 .and. nzz .gt. 1) then
c
!         write (6,*) 'SOLVER: starting z-sweep'
c
!         do j=js, je
          do j=j1, j2
            call zlagrange_sweep(j, d, e, u, v, w, ge, in, jn, kn,
     &        gravity, gr_zacc, idual, eta1, eta2,
     &        is, ie, js, je, ks, ke,
     &        gamma, pmin, dt, dx, dy, dz,
     &        idiff, iflatten, isteepen, ipresfree,
     &        nsubgrids, lface, rface,
     &        fistart, fiend, fjstart, fjend,
     &        dindex, eindex, geindex,
     &        uindex, vindex, windex, array,
     &        ncolor, colorpt, coloff, colindex,
     &        tmp(1+ms*0),  tmp(1+ms*1),  tmp(1+ms*2),  tmp(1+ms*3),
     &        tmp(1+ms*4),  tmp(1+ms*5),  tmp(1+ms*6),  tmp(1+ms*7), 
     &        tmp(1+ms*8),  tmp(1+ms*9),  tmp(1+ms*10), tmp(1+ms*11),
     &        tmp(1+ms*12), tmp(1+ms*13), tmp(1+ms*14), tmp(1+ms*15),
     &        tmp(1+ms*16), tmp(1+ms*17), tmp(1+ms*18), tmp(1+ms*19),
     &        tmp(1+ms*20), tmp(1+ms*21), tmp(1+ms*22), tmp(1+ms*23),
     &        tmp(1+ms*24), tmp(1+ms*25), tmp(1+ms*26), tmp(1+ms*27),
     &        tmp(1+ms*28), tmp(1+ms*29),
     &        tmp(1+ms*(30+0*ncolor)), tmp(1+ms*(30+1*ncolor)), 
     &        tmp(1+ms*(30+2*ncolor)), tmp(1+ms*(30+3*ncolor))
     &                       )
          enddo
c
        endif
c
      enddo
c
      return
      end

