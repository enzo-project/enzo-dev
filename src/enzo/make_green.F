c=======================================================================
c//////////////////////  SUBROUTINE MAKE_GREEN  \\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine make_green(green, nx, ny, nz, in, jn, kn, nalias, 
     &                      S2Part, refinement,
     &                      S2Table, S2Min, S2Step,
     &                      SinTable, SinMin, SinStep)
c
c  COMPUTES OPTIMAL GREENS FUNCTION FOR TSC CLOUDS
c
c     written by: Greg Bryan
c     date:       July, 1995
c
c  PURPOSE: Computes a greens function in k space for the isolated/periodic
c           coloumb potential (on a grid).  The resoluting function can
c           be convolved in real space (i.e. multiplied in k-space) with
c           the produce to produce a potential.  This version uses the
c           Hockney & Eastwood (Computer Simulation using Particles, 1980)
c           optimal influence function technique (see also Eastwood and
c           Brownrigg, J. Comp. Phys. 1979, 32, 24.)  The influence function
c           therefore depends on the other techniques used to solve for
c           the gravitational force.  Here we assume:
c           I. The force is calculated by: Direct K-space solver
c              Choices are:
c               a) Direct k-space solver:  D(k) = -ik
c               b) A four point finite difference:
c                  D(k) = -2i*   alpha *sin(kH/2)*cos(kH/2)/H -
c                           i*(1-alpha)*sin(kH  )*cos(kH  )/H
c               c) A two point finite difference:
c                  D(k) = -2i*sin(kH/2)*cos(kH/2)/H
c           II.  The assignment scheme is: TSC
c              Choices are:
c              i) The NGP(p=1), CIC(p=2) or TSC(p=3) techniques:
c                  W(k) = (sin(kH/2)/(kH/2))^p
c           NOTE: This routine assumes equal-sized zones.
c           Define assignment scheme (NGP, CIC or TSC)
c           Define a finite difference order (TWO_POINT_FD, FOUR_POINT_FD or
c                  DIRECT_KSPACE).
c           The units used are such that when G(k) is multiplied by by
c           d(k), the transformed density field, the result is the potential
c           times the grid spacing (i.e. divide by the grid spacing to get
c           the potential in the correct units).
c
c  INPUTS:
c        i,j,kn      = real dimensions of green
c        nalias      = number of aliases in sum (-nalias -> +nalias)
c        nx,ny,nz    = active dimensions of green
c        refinement  = refinement factor
c        S2Part      = S2 particle size in grid zones
c        S2Table     = Table of precomputed S2 values
c        S2Min       = minimum x value in S2Table
c        S2Step      = x step in S2Table
c        SinTable    = Table of precomputed Sin values
c        SinMin      = minimum x value in SinTable
c        SinStep     = x step in SinTable
c
c  Outputs:
c        green       = Greens function
c
c  LOCALS:
c        num_dim     = number of dimensions to be used for force law
c        nx,y,zmid   = midpoint (+1) of the grid in each axis
c        nx,y,zd2    = number of grid points divided by 2 for each axis
c
c  EXTERNALS: 
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c  arguments
c
      integer in, jn, kn, nalias, nx, ny, nz, refinement
      real*8    S2Min, S2Part, S2Step, SinMin, SinStep
      real*8    S2Table(1), SinTable(1), green(nx/2+1, jn, kn)
c
c  locals
c
      integer bi, bj, bk, index, i, j, k, ii, jj, kk,
     &        naliasx, naliasy, naliasz, nxmid, nymid, nzmid
      real*8    dksqr, keta2, kfact, kmod, rk, sk, sktop, twopi
      real*8    k1, k2, k3, kx, kxbern, ky, kybern, kz, kzbern,
     &        sinx2, siny2, sinz2, sinxb, sinyb, sinzb, 
     &        usqr1, usqr2, usqr3, usqrx, usqry, usqrz
c
c  Define table lookup function
c
      real*8    Table1, Table2, Step, Min, Tablex, TableLookUp
      integer Tablei
      TableLookUp(Table1, Table2, Step, Min, Tablei, Tablex) = 
     &          Table1 + (Tablex - real(Tablei-1)*Step - Min) 
     &          / Step * (Table2 - Table1)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Set constants
c
      twopi  = 8.0*atan(1.0)
c
c     Set maximum wavenumbers in the positive direction
c
      nxmid  = max(nx/2 + 1, 1)
      nymid  = max(ny/2 + 1, 1)
      nzmid  = max(nz/2 + 1, 1)
c
c     Set number of aliases for each dimension
c
      naliasx = nalias
      naliasy = nalias
      naliasz = nalias
      if (nx .eq. 1) naliasx = 0
      if (ny .eq. 1) naliasy = 0
      if (nz .eq. 1) naliasz = 0
c
c     Loop over all wavenumbers and compute optimal influence function.
c     First, compute SUM[U^2] (e.q. 8-45)
c
      do k = 1, nz
         kk = k-1
         if (k .gt. nzmid) kk = kk - nz
         kz    = real(kk)*twopi/real(nz)
         index = int((0.5*kz - SinMin)/SinStep)
         sinz2 = TableLookUp(SinTable(index), SinTable(index+1),
     &                       SinStep, SinMin, index, 0.5*kz)
c         write (6,*) sinz2, sin(0.5*kz), sinz2/sin(0.5*kz)
         usqrz = 1.0 - sinz2**2 + 2.0/15.0*sinz2**4
         if (k .eq. 1) usqrz = 1.0
c
         do j = 1, ny
            jj = j-1
            if (j .gt. nymid) jj = jj - ny
            ky    = real(jj)*twopi/real(ny)
            index = int((0.5*ky - SinMin)/SinStep)
            siny2 = TableLookUp(SinTable(index), SinTable(index+1), 
     &                          SinStep, SinMin, index, 0.5*ky)
            usqry = 1.0 - siny2**2 + 2.0/15.0*siny2**4
            if (j .eq. 1) usqry = 1.0
c
            do i = 1, nxmid
               ii    = i-1
               kx    = real(ii)*twopi/real(nx)
               index = int((0.5*kx - SinMin)/SinStep)
               sinx2 = TableLookUp(SinTable(index), SinTable(index+1),
     &                             SinStep, SinMin, index, 0.5*kx)
               usqrx = 1.0 - sinx2**2 + 2.0/15.0*sinx2**4
               if (i .eq. 1) usqrx = 1.0
c
c              Compute D^2 (eq. 8-46)
c
               dksqr = kx**2 + ky**2 + kz**2
               if (dksqr .eq. 0.0) dksqr = 1.0
c
c              Clear vector sum
c
               k1 = 0.0
               k2 = 0.0
               k3 = 0.0
c
c              Sum over Bernoulli (sp?) zones (c.f. eq. 8-22)
c
               do bk = -naliasz, naliasz
                  kzbern = kz + twopi*real(bk)
                  do bj = -naliasy, naliasy
                     kybern = ky + twopi*real(bj)
                     do bi = -naliasx, naliasx
                        kxbern = kx + twopi*real(bi)
c
c                       Compute the modulus of k
c
                        kmod  = sqrt(kxbern**2 + kybern**2 + kzbern**2)
                        if (kmod .eq. 0.0) kmod = 1.0
c
c                       Compute Sk (eq. 8-41)
c
                        keta2 = 0.5*kmod*S2part
                        index = int((keta2 - S2Min)/S2Step)
                        sk    = TableLookUp(S2Table(index), 
     &                    S2Table(index+1), S2Step, S2Min, index, keta2)
                        sk    = sk**2
c
c                       Compute R(k), the reference force
c
                        rk = 1.0/kmod**2
c
c                       If this is a subgrid, compute Sk for the larger
c                          S2 particle size
c
                        if (refinement .ne. 1) then
                           keta2 = 0.5*kmod*S2Part*real(refinement)
                           index = int((keta2 - S2Min)/S2Step)
                           sktop = TableLookUp(S2Table(index),
     &                                         S2Table(index+1), S2Step,
     &                                         S2Min, index, keta2)
                           sk    = (sk - sktop**2)
                        endif
c
c                       Lookup sin values (for next step).
c
c                        if (min(kxbern,kybern,kzbern) .lt. SinMin) then
c                           write (6,*) 'make_green: out of range:', 
c     &                                  kxbern, kybern, kzbern
c                           stop
c                        endif
                        index = int((0.5*kxbern - SinMin)/SinStep)
                        sinxb = TableLookUp(SinTable(index), 
     &                                      SinTable(index+1), SinStep,
     &                                      SinMin, index, 0.5*kxbern)
                        index = int((0.5*kybern - SinMin)/SinStep)
                        sinyb = TableLookUp(SinTable(index),
     &                                      SinTable(index+1), SinStep, 
     &                                      SinMin, index, 0.5*kybern)
                        index = int((0.5*kzbern - SinMin)/SinStep)
                        sinzb = TableLookUp(SinTable(index),
     &                                      SinTable(index+1), SinStep, 
     &                                      SinMin, index, 0.5*kzbern)
c
c                       Compute U^2 (eq. 8-42)
c
                        usqr1 = sinxb/(0.5*kxbern + tiny)
                        usqr2 = sinyb/(0.5*kybern + tiny)
                        usqr3 = sinzb/(0.5*kzbern + tiny)
                        if (kxbern .eq. 0.0) usqr1 = 1.0
                        if (kybern .eq. 0.0) usqr2 = 1.0
                        if (kzbern .eq. 0.0) usqr3 = 1.0
c
c                       Combine and sum (part of HE, eq. 8-22)
c
                        kfact = sk*rk*(usqr1*usqr2*usqr3)**6
c
                        k1 = k1 + kxbern*kfact
                        k2 = k2 + kybern*kfact
                        k3 = k3 + kzbern*kfact
c
                     enddo
                  enddo
               enddo
c
c        Compute G(k) (eq. 8-22) (k1,2,3 is SUM[U^2 R(k)])
c
               green(i,j,k) =  -(k1*kx + k2*ky + k3*kz)
     &                    /(dksqr*(usqrx*usqry*usqrz)**2)
c
c        Approximate (poor man''s) greens function
c
c#define APPROXIMATE_GREENS_FUNCTION
#ifdef APPROXIMATE_GREENS_FUNCTION
c
              green(i,j,k) = -1.0/(kx**2 + ky**2 + kz**2 + tiny)
              if (real(ii*ii + jj*jj + kk*kk) .gt. 
     &            0.5*max(nx, ny, nz)**2) green(i,j,k) = 0.0
c
#endif /* APPROXIMATE_GREENS_FUNCTION */
c
            enddo
         enddo
      enddo
c
c     Clear the zero wavenumber position
c
      green(1,1,1) = 0.0
c      green(1,1,1) = 1.0
c
      return
      end
