*****************************************************************************
*                                                                           *
* Copyright 2004 Greg Bryan                                                 *
* Copyright 2004 Laboratory for Computational Astrophysics                  *
* Copyright 2004 Board of Trustees of the University of Illinois            *
* Copyright 2004 Regents of the University of California                    *
*                                                                           *
* This software is released under the terms of the "Enzo Public License"    *
* in the accompanying LICENSE file.                                         *
*                                                                           *
*****************************************************************************

#include "fortran.def"
c=======================================================================
c/////////////////////////  SUBROUTINE CALCDIFF  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine calcdiss(
     &            dslice, eslice, uslice, v, w, pslice, dx, dy, dz,
     &            idim, jdim, kdim, i1, i2, j1, j2, k, nzz, idir,
     &            dimx, dimy, dimz, dt, gamma, idiff, iflatten,
     &            diffcoef, flatten
     &                    )
c
c  COMPUTES THE DIFFUSION COEFFICIENTS AND FLATTENING FOR PPM (CW84)
c
c  written by: Greg Bryan
c  date:       July, 1994
c  modified1:  Robert Harkness
c  date:       July, 2003
c              Mods to prevent out-of-bounds memory references with
c              diffusion.  Works with new euler sweep range mod by
c              Alexei Kritsuk to at least allow for flattening.
c              Calculation for idiff = 2 was incorrect in any case
c  modified2:  James Bordner
c  date:       2003-11-21
c              Fixed minor index bug; moved some loop invariants outside
c              inner loops; general housecleaning
c  modified3:  Alexei Kritsuk
c  date:       2004-08-10
c              Fixed iflatten = 1 procedure (A.2), shock capturing (A.1);
c              cleaned and enabled iflatten = 3 (A.7-10)
c
c  PURPOSE:  Given the values of idiff and iflatten, we calculate
c     the diffusion coefficients and amount of local flattening using
c     one of a number of routines.  They are:
c     a) idiff = 0  -> no difussion
c     b) idiff = 1  -> simple diffusion
c       Computes the diffusion coefficient given in Colella &
c       Woodward eq. A.3 for use later when calculating the diffusive
c       flux.  The essential logic is given in 4.4:
c
c            diffcoef =  K max|-Div u(j+1/2), 0|
c
c       where Div is the discrete undivided difference approximation to
c       the multidimensional divergence of u.  This version calculates
c       the diffusion coefficient for a single slice (k).
c     c) idiff = 2  -> alternate diffusion
c       Here we use the scheme described in CW84, eq. A11.  Although,
c       this is intended for local grid motion, here we adapt it for
c       the local amount of explicit diffusion.
c     d) iflatten = 0 -> no flattening
c     e) iflatten = 1 -> The simple flattening scheme described in eqs
c        A1 & A2 of CW84.
c     f) iflatten = 2 -> Flattening needed for Lagrangean hydrodynamics
c        (note: not L+remap).  Stuck in for completeness, but not tested.
c        This is described in equations. A.4 - A.6.
c     g) iflatten = 3 -> Flattening recomended for multidimensional
c        calculations.  Described in CW84 A7-A10.
c     NOTE: this routine is intended to be called for x,y and z directions,
c        cyclically permuting the arguments to obtain the correct result
c        for each direction (i.e. call with dx,dy,dz,... for x-dir, dy,dz,dx
c        for y-dir, and dz,dx,dy for z-dir).  Also, insure idir is set
c        correctly, since a few calculation must know which sweep
c        direction this is.  
c           permute: v,w,dx,dy,dz,i1,i2,j1,j2,k,nzz,idim,jdim,kdim
c           do not permute: dimx,dimy,dimz
c
c  INPUT ARGUMENTS:
c     d      - density field
c     dimx,y,z - dimensions of v and w field arrays
c     dslice - single slice of the density in (sweep, orthodim) direction
c     dx     - distance between Eulerian zone edges in direction 1
c     dy,dz  - distance between Eulerian zone edges in 2,3 directions
c     eslice - single slice of the total energy in (sweep, orthodim) direction
c     i1,i1  - active range in direction-1
c     idiff  - type of explicit diffusion to be used (see above)
c     idir   - direction of sweep (1=x, 2=y, 3=z)
c     idim   - first dimension of slices (direction-1)
c     iflatten - type of local flattening to be used (see above)
c     j2,j2  - active range in direction-2
c     jdim   - second dimension of slices (direction-2)
c     k      - current slice index
c     kdim   - dimension of direction-3
c     nzz    - number of active zones in the z-direction
c     pslice - single slice of the pressure in (sweep, orthodim) direction
c     uslice - single slice of the pressure in (sweep, orthodim) direction
c     v      - 2-velocity field
c     w      - 3-velocity field
c
c  OUTPUT ARGUMENTS:
c     diffcoef - diffusion coefficient in slice k
c     flatten  - amount of local flattening (0-1)
c
c  LOCALS:
c     di     - temporary for 1/dslice
c     wflag  - CW84 eq. A1
c
c  PARAMETERS
c     K      - free parameter in diffusion coefficient calculation [0.1]
c     (for other parameters see CW84)
c
c-----------------------------------------------------------------------
      implicit NONE
c
      integer ijkn
      parameter (ijkn=MAX_ANY_SINGLE_DIRECTION)
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer dimx, dimy, dimz, i1, i2, idiff, idim, idir, iflatten, 
     &        j1, j2, jdim, k, kdim, nzz
      real    gamma, dt, dx(idim), dy(jdim), dz(kdim)
      real    dslice(idim,jdim), eslice(idim,jdim), uslice(idim,jdim),
     &        diffcoef(idim,jdim), flatten(idim,jdim), pslice(idim,jdim)
      real    v(dimx,dimy,dimz), w(dimx,dimy,dimz)
c
c  locals
c
      integer i, j, nyz, is, isp, ism
      real    cj2s, di(ijkn), difftemp(ijkn), flattemp(ijkn),
     &        kappa(ijkn), kappa_tilde, omega(ijkn),
     &        qa, qb, qc, s, sigma(ijkn),
     &        sigma_tilde, vdiff1(ijkn), vdiff2(ijkn), wdiff1(ijkn),
     &        wdiff2(ijkn), wflag(ijkn), Z, ZE
      logical ljx1,ljy1,ljz1,ljx2,ljy2,ljz2,lj1,lj2
      logical lkx1,lky1,lkz1,lkx2,lky2,lkz2,lk1,lk2
      real dp1,dp2,de1,de2,dee,dpp,omega_tilde
c
c  parameters
c
      real    epsilon, kappa1, kappa2, Kparam, nu1, nu2, nu3,
     &        omega1, omega2, sigma1, sigma2
      parameter (epsilon = 0.33  , kappa1 = 2.0 , kappa2 =  0.01,
     &            Kparam = 0.1   ,    nu1 = 2.0 , nu2    =  0.1 ,
     &               nu3 = 0.2333, omega1 = 0.75, omega2 = 10.0 ,
     &            sigma1 = 0.5   , sigma2 = 1.0)
c
c Notice, for 1D calculations nu2 = 0 and nu3 = 0.3333, see CW84, p. 200.
c
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c  Compute number of active zones in transverse direction
c
      nyz = j2 - j1 + 1
c
c  Precompute logical masks for k loop bounds
c
      lk1  = (1.lt.k .and. k.lt.kdim)
      lkx1 = (1.lt.k .and. k.lt.dimx)
      lky1 = (1.lt.k .and. k.lt.dimy)
      lkz1 = (1.lt.k .and. k.lt.dimz)

      lk2  = (2.lt.k .and. k.lt.kdim-1)
      lkx2 = (2.lt.k .and. k.lt.dimx-1)
      lky2 = (2.lt.k .and. k.lt.dimy-1)
      lkz2 = (2.lt.k .and. k.lt.dimz-1)
c
c  Loop across the k slice
c
      do j=j1, j2
c
c  Precompute logical masks for j loop bounds
c
         lj1  = (1.lt.j .and. j.lt.jdim)
         ljx1 = (1.lt.j .and. j.lt.dimx)
         ljy1 = (1.lt.j .and. j.lt.dimy)
         ljz1 = (1.lt.j .and. j.lt.dimz)

         lj2  = (2.lt.j .and. j.lt.jdim-1)
         ljx2 = (2.lt.j .and. j.lt.dimx-1)
         ljy2 = (2.lt.j .and. j.lt.dimy-1)
         ljz2 = (2.lt.j .and. j.lt.dimz-1)
c
c  All (well, almost all) routines below need this quantity (CW84, eq. a1)
c
         do i=i1-2, i2+2
            qb =   abs(pslice(i+1,j) - pslice(i-1,j))  ! AK abs()
     &           / min(pslice(i+1,j),  pslice(i-1,j))
            if ( qb.gt.epsilon .and.
     &          uslice(i-1,j).gt.uslice(i+1,j) ) then
               wflag(i) = 1.0
            else
               wflag(i) = 0.0
            endif
         end do
c
c  Diffusion routine (B) needs this calculation, which depends on v & w,
c     so it depends on the sweep direction.
c     
         if (idiff .eq. 1) then
c     
            if (idir .eq. 1) then

               if (dimy.gt.1 .and. ljy1) then
                  do i = i1-2, i2+2
                     vdiff1(i) = (v(i,j-1,k) + v(i-1,j-1,k)) 
     &                    -      (v(i,j+1,k) + v(i-1,j+1,k))
                  end do
               else
                  do i = i1-2, i2+2
                     vdiff1(i) = 0.0
                  end do
               end if
c     
               if (dimz.gt.1 .and. lkz1) then
                  do i = i1-2, i2+2
                     wdiff1(i) = (w(i,j,k-1) + w(i-1,j,k-1))
     &                    -      (w(i,j,k+1) + w(i-1,j,k+1))
                  end do
               else 
                  do i = i1-2, i2+2
                     wdiff1(i) = 0.0
                  end do
               end if
c     
c     
            elseif (idir .eq. 2) then
c     
c     
               if (dimz.gt.1 .and. ljz1) then
                  do i = i1-2, i2+2
                     vdiff1(i) = (v(k,i,j-1) + v(k,i-1,j-1))
     &                    -      (v(k,i,j+1) + v(k,i-1,j+1))
                  end do
               else
                  do i = i1-2, i2+2
                     vdiff1(i) = 0.0
                  end do
               endif
c     
               if (dimx.gt.1 .and. lkx1) then
                  do i = i1-2, i2+2
                     wdiff1(i) = (w(k-1,i,j) + w(k-1,i-1,j))
     &                    -      (w(k+1,i,j) + w(k+1,i-1,j))
                  end do
               else
                  do i = i1-2, i2+2
                     wdiff1(i) = 0.0
                  end do
               endif
c
c     
            elseif (idir .eq. 3) then
c
c     
c     
               if (dimx.gt.1 .and. ljx1) then
                  do i = i1-2, i2+2
                     vdiff1(i) = (v(j-1,k,i) + v(j-1,k,i-1)) 
     &                    -      (v(j+1,k,i) + v(j+1,k,i-1))
                  end do
               else
                  do i = i1-2, i2+2
                     vdiff1(i) = 0.0
                  end do
               end if
c     
               if (dimy.gt.1 .and. lky1) then
                  do i = i1-2, i2+2
                     wdiff1(i) = (w(j,k-1,i) + w(j,k-1,i-1)) 
     &                    -      (w(j,k+1,i) + w(j,k+1,i-1))
                  end do
               else
                  do i = i1-2, i2+2
                     wdiff1(i) = 0.0
                  end do
               endif
c     
            end if
c
         endif
c
c  Diffusion routine (C) needs this calculation, which depends on v & w,
c    so it depends on the sweep direction.
c     
         if (idiff .eq. 2) then
c     
            if (idir .eq. 1) then
c     
               if (dimy.gt.1 .and. ljy2) then
                  do i=i1-2, i2+2
                     vdiff2(i) = (v(i,j-2,k) + v(i-1,j-2,k))
     &                    -      (v(i,j+2,k) + v(i-1,j+2,k))
                  end do
               else
                  do i=i1-2, i2+2
                     vdiff2(i) = 0.0
                  end do
               endif
c     
               if (dimz.gt.1 .and. lkz2) then
                  do i=i1-2, i2+2
                     wdiff2(i) = (w(i,j,k-2) + w(i-1,j,k-2)) 
     &                    -      (w(i,j,k+2) + w(i-1,j,k+2))
                  end do
               else
                  do i=i1-2, i2+2
                     wdiff2(i) = 0.0
                  end do
               endif
c     
c     
            elseif (idir .eq. 2) then
c     
               if (dimz.gt.1 .and. ljz2) then
                  do i=i1-2, i2+2
                     vdiff2(i) = (v(k,i,j-2) + v(k,i-1,j-2)) 
     &                    -      (v(k,i,j+2) + v(k,i-1,j+2))
                  end do
               else
                  do i=i1-2, i2+2
                     vdiff2(i) = 0.0
                  end do
               endif
c     
               if (dimx.gt.1 .and. lkx2) then
                  do i=i1-2, i2+2
                     wdiff2(i) = (w(k-2,i,j) + w(k-2,i-1,j)) 
     &                    -      (w(k+2,i,j) + w(k+2,i-1,j))
                  end do
               else
                  do i=i1-2, i2+2
                     wdiff2(i) = 0.0
                  end do
               endif
c     
c
            elseif (idir .eq. 3) then
c     
               if (dimx.gt.1 .and. ljx2) then
                  do i = i1-2, i2+2
                     vdiff2(i) = (v(j-2,k,i)+v(j-2,k,i-1)) 
     &                    -      (v(j+2,k,i)+v(j+2,k,i-1))
                  end do
               else
                  do i = i1-2, i2+2
                     vdiff2(i) = 0.0
                  end do
               end if
c     
               if (dimy.gt.1 .and. lky2) then
                  do i = i1-2, i2+2
                     wdiff2(i) = (w(j,k-2,i)+w(j,k-2,i-1)) 
     &                    -      (w(j,k+2,i)+w(j,k+2,i-1))
                  end do
               else
                  do i = i1-2, i2+2
                     wdiff2(i) = 0.0
                  end do
               end if
c     
c                
            endif
c
         endif
c  
c-----------------------------------------------------------------------
c (B)  Compute diffusion coefficient
c
         if (idiff .eq. 1) then
c
            do i=i1, i2+1
               diffcoef(i,j) = uslice(i-1,j) - uslice(i,j)
            end do
c     
c
            if (nyz.gt.1 .and. lj1) then
               do i=i1, i2+1
                  diffcoef(i,j) = diffcoef(i,j) + (0.25*(dx(i)+dx(i-1))
     &                 /   (0.5*(dy(j+1)+dy(j-1)) + dy(j)))*vdiff1(i)
               end do
            endif
c     
            if (nzz.gt.1 .and. lk1) then
               do i=i1, i2+1
                  diffcoef(i,j) = diffcoef(i,j) + (0.25*(dx(i)+dx(i-1))
     &                 /   (0.5*(dz(k+1)+dz(k-1)) + dz(k)))*wdiff1(i)
               end do
            endif
c     
            do i=i1, i2+1
               diffcoef(i,j) = Kparam*max(0.0, diffcoef(i,j))
            end do
c
         endif
c
c-----------------------------------------------------------------------
c (E)  Construct flattening parameter (eqns A1 and A2), also 68-74
c
         if (iflatten .eq. 1) then

            do i=i1-1, i2+1
               if (abs(pslice(i+2,j) - pslice(i-2,j))/
     &            min(pslice(i+2,j),pslice(i-2,j)) .lt. epsilon) then
                 qa = 1.0
               else
                 qa =     ( pslice(i+1,j) - pslice(i-1,j) )
     &                /   ( pslice(i+2,j) - pslice(i-2,j) )
               endif
               flattemp(i) = min( 1.0, (qa-omega1)*omega2*wflag(i) )
c               flattemp(i) = wflag(i) !AK: Svetsov 2001, ShockWaves 11, 229
               flattemp(i) = max(0.0, flattemp(i))
c              flattemp(i) = 1.0-flattemp(i)  !AK: typo in CW84 eq. (A.2)
             end do

            flattemp(i1-2) = flattemp(i1-1)
            flattemp(i2+2) = flattemp(i2+1)

c     Now, choose the maximum (eq. A2, first part)

            do i=i1-1, i2+1
               if (pslice(i+1,j) - pslice(i-1,j) .lt. 0.0) then
                  flatten(i,j) = max(flattemp(i),flattemp(i+1))
               else
                  flatten(i,j) = max(flattemp(i),flattemp(i-1))
               endif
c        There are some additional recipes used by others  !AK
c               flatten(i,j) = 0.025 !Runacres+1 astro-ph/0405315v1, p.7;
c               flatten(i,j) = 1.0 !Godunov interpolation (constant)
c               flatten(i,j) = max(flatten(i,j), 0.1)
            end do

         endif

c-----------------------------------------------------------------------
c (F) Construct second type flattening parameter (eq. A4-A6)
c
         if (iflatten .eq. 2) then

            do i=i1-3, i2+3
               di(i) = 1.0/dslice(i,j)
            end do

            do i=i1-1, i2+1

               is = i + sign(2.0, pslice(i+1,j) - pslice(i-1,j))

               omega(i) = max(0.0, omega1 * (omega2 
     &              - (pslice(i+1,j) - pslice(i-1,j))
     &              / (pslice(i+2,j) - pslice(i-2,j)) ))

               Z = sqrt(( max(pslice(i+2,j),pslice(i-2,j)) +
     &              0.5*(pslice(i+2,j)+pslice(i-2,j))*(gamma-1.0))
     &              / max(di(i+2),di(i-2)))

               kappa_tilde = 
     &              (Z + sqrt(gamma*pslice(is,j)*dslice(is,j))) / Z

               kappa(i)   = max(0.0, (kappa_tilde - kappa1)
     &              /                (kappa_tilde + kappa2))

               flattemp(i) = min(wflag(i)*omega(i), kappa(i))

            end do

            flattemp(i1-2) = flattemp(i1-1)
            flattemp(i2+2) = flattemp(i2+1)

c     Choose minimum (eq. A6)

            do i=i1-1, i2+1
               flatten(i,j) = 
     &              max(flattemp(i-1),flattemp(i),flattemp(i+1))
            end do

         endif
c-----------------------------------------------------------------------
c (G) Construct third type flattening parameter (eq. A7-A9)
c
         if (iflatten .eq. 3 .or. idiff .eq. 2) then
c     
            do i=i1-3, i2+3
               di(i) = 1.0/dslice(i,j)
            end do
c     
            do i=i1-1, i2+1
c     
               dp1 = pslice(i+1,j) - pslice(i-1,j)
               dp2 = pslice(i+2,j) - pslice(i-2,j)
               de1 = eslice(i+1,j) - eslice(i-1,j)
               de2 = eslice(i+2,j) - eslice(i-2,j)
c
               dpp = 0.0
               dee = 0.0
               if (dp2 .ne. 0.0) dpp = dp1/dp2
               if (de2 .ne. 0.0) dee = de1/de2
               omega_tilde = max( dpp , dee )
c
               ism = i + sign(2.0, dp1)        ! post-shock
               isp = i - sign(2.0, dp1)        ! upstream
               s   =   - sign(1.0, dp1)        ! i+s - zone upstream from i
               if (dp1 .eq. 0.0) s = 0.0
c
c              Compute sigma(i): strength of a shock near zone i
c     
               sigma_tilde = 
     &              wflag(i)*abs(dp2)/min(pslice(i+2,j),pslice(i-2,j))
               sigma(i) = max( 0.0, (sigma_tilde - sigma1)
     &              /               (sigma_tilde + sigma2) )
c
c              Compute omega(i): steepness of a shock near zone i 
c
               omega(i) = 
     &              max(0.0, omega2*(omega_tilde - omega1)) ! as in Prometheus
c     &              max(0.0, omega1*(omega2 - omega_tilde)) !AK: a CW84 typo?
c     
c              Compute an estimate of the Lagrangean shock speed W_j
c
               Z = sqrt(( max(pslice(i+2,j),pslice(i-2,j)) +
     &              0.5*(pslice(i+2,j)+pslice(i-2,j))*(gamma-1.0))
     &              / max(di(i+2),di(i-2)))
c
c              Compute W_j^E; Possible typos in kappa_tilda definition in 
c              CW84 (A.9). !AK
c
               ZE    = s*Z/dslice(ism,j) + uslice(ism,j) + tiny
               cj2s = sqrt(gamma*pslice(isp,j)/dslice(isp,j))
c
c              Compute kappa(i): the wavelength of the noise
c
               kappa_tilde = abs((ZE - uslice(isp,j) + s*cj2s)/ZE)
               kappa(i) = max(0.0, (kappa_tilde - kappa1)
     &              /              (kappa_tilde + kappa2))
c     
               flattemp(i) =
     &              min(kappa(i),wflag(i)*omega(i),wflag(i)*sigma(i))
c     
            end do
c     
         endif
c     
         if (iflatten .eq. 2 .or. iflatten .eq. 3) then
c     
            flattemp(i1-2) = flattemp(i1-1)
            flattemp(i2+2) = flattemp(i2+1)
c     
c     Now, choose the maximum (eq. A.2, first part; also A.10)
c     
            do i=i1-1, i2+1
               if (pslice(i+1,j) - pslice(i-1,j) .lt. 0.0) then
                  flatten(i,j) = max(flattemp(i),flattemp(i+1))
               else
                  flatten(i,j) = max(flattemp(i),flattemp(i-1))
               endif
            end do
c     
         endif
c     
c-----------------------------------------------------------------------
c (C)  Compute complex diffusion coefficient
c     

         if (idiff .eq. 2) then

            do i=i1-1, i2+1
               difftemp(i) = uslice(i-2,j) - uslice(i+1,j)
            end do
            
            if (nyz.ne.1 .and. lj2) then
               do i=i1-1, i2+1
                  qa = 0.5*(0.5*(dx(i+1)+dx(i-2))+dx(i)+dx(i-1))
     &                 / (0.5*(dy(j+2)+dy(j-2)) + dy(j+1)+dy(j)+dy(j-1))
                  difftemp(i) = difftemp(i) + qa*vdiff2(i)
               end do
            end if

            if (nzz.ne.1 .and. lk2) then
               do i=i1-1, i2+1
                  qb = 0.5*(0.5*(dx(i+1)+dx(i-2))+dx(i)+dx(i-1))
     &                 / (0.5*(dz(k+2)+dz(k-2)) + dz(k+1)+dz(k)+dz(k-1))
                  difftemp(i) = difftemp(i) + qb*wdiff2(i)
               end do
            end if

            do i=i1-1, i2+1
               difftemp(i) = max(0.0, nu1*difftemp(i))
            end do

            difftemp(i1-2) = difftemp(i1-1)
            difftemp(i2+2) = difftemp(i2+1)

            do i=i1-1, i2+1

               s  = - sign(1.0, pslice(i-1,j) - pslice(i,j))
               diffcoef(i,j) = s*min(
     &              max(difftemp(i-1),difftemp(i),difftemp(i+1)),
     &              (nu2 + nu3*max(sigma(i-1)*kappa(i-1)**2
     &              ,              sigma(i  )*kappa(i  )**2 ))
     &              * dx(i+1)*s/dt)

            end do

         endif
c
c  Done
c
      end do
c
      return
      end
