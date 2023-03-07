#include "fortran.def"
#include "phys_const.def"
#include "error.def"
c
c=======================================================================
c/////////////////////  SUBROUTINE STAR_FEEDBACK \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine star_feedback_emf(nx, ny, nz,
     &                      d, te, ge, u, v, w,
     &                      idual, imethod, dt, dx, t,
     &                      d1, x1, v1, t1,
     &                      npart, xstart, ystart, zstart, ibuff,
     &                      xp, yp, zp, up, vp, wp,
     &                      mp, iminit, minit, tcp, type,
     &                      alpha, tfb, p0)

c
c  RELEASES "STELLAR" WINDS
c
c  Calculates early, pre-SNe momentum injection following 
c  Keller, Kruijssen, & Chevance 2022 and deposits it as in kinetic
c  feedback algorithm from Simpson et al 2015 (star_maker3mom.F)
c
c  INPUTS:
c
c    d     - density field
c    te,ge - total energy and gas energy fields
c    u,v,w - velocity fields
c    dt    - current timestep
c    dx    - zone size (code units)
c    t     - current time
c    d1,x1,v1,t1 - factors to convert d,dx,v,t to physical units
c    nx,ny,nz - dimensions of field arrays
c    ibuff    - number of buffer zones at each end of grid
c    idual    - dual energy flag
c    imethod  - hydro method (0 - PPMDE, 1 - PPMLR, 2 - ZEUS)
c    npart    - particle array size specified by calling routine
c    x/y/z start - starting position of grid origin
c    xp,yp,zp - positions of created particles
c    up,vp,wp - velocities of created particles
c    mp       - mass of new particles
c    tcp      - creation time of particle (-1 if not a star particle)
c    type     - particle type
c    iminit   - whether initial particle mass is stored separately.
c               This is *required* to be true
c    minit    - initial particle mass
c    alpha    - free parameter for wind driving mechanism
c    tfb      - feedback timescale parameter
c    p0       - momentum scale parameter
c
c  OUTPUTS:
c    d,u,v,w,ge,e - modified field
c
c
c-----------------------------------------------------------------------
       implicit none
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC nx, ny, nz, ibuff, npart, idual, imethod, iminit
      R_PREC    d(nx,ny,nz), dm(nx,ny,nz), te(nx,ny,nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    ge(nx,ny,nz)
      R_PREC    dt, dx, t
      R_PREC    d1, x1, v1, t1, justburn
      P_PREC    xstart, ystart, zstart
      P_PREC    xp(npart), yp(npart), zp(npart)
      R_PREC    up(npart), vp(npart), wp(npart)
      R_PREC    mp(npart), tcp(npart), minit(npart)
      R_PREC    alpha, tfb, p0
      INTG_PREC type(npart)
c
c  Locals
c
      INTG_PREC n, ic, jc, kc, ip, jp, kp,
     &     iface, jface, kface, i, j, k
      R_PREC dxf, dyf, dzf, dxc, dyc, dzc, xfc, yfc, zfc,
     &     xfcshift, yfcshift, zfcshift, face_shift,
     &     xface, yface, zface, xpos, ypos, zpos, fbuff
      R_PREC delta_p, delta_ke, ke_thermalized
      R_PREC ke_before(4,4,4), ke_after
c
c-----------------------------------------------------------------------
c
c     Loop over particles
c
      do n=1, npart
         if (tcp(n) .gt. 0 .and. mp(n) .gt. 0 .and. type(n) .eq. 2 
     &       .and. iminit .eq. 1) then
c         
c         Calculate the amount of momentum deposited per cell
c           Use Equation 10 from Keller et al. 2022
c           Note that particle mass is actually a density
c
            delta_p = alpha * p0 * minit(n)*dx**3 
     &         * ( ((t+dt)/tfb)**(4*alpha-1) - (t/tfb)**(4*alpha-1) )
     &         / 27
c
c          skip if very little momentum deposited.
c           TODO adjust
            if (delta_p/d(i,j,k) .lt. 1.e-10_RKIND) goto 10
c
c         Compute index of the cell that the star particlet
c           resides in.
c 
            ip = int((xp(n) - xstart)/dx,IKIND) + 1_IKIND
            jp = int((yp(n) - ystart)/dx,IKIND) + 1_IKIND
            kp = int((zp(n) - zstart)/dx,IKIND) + 1_IKIND
c
c	        Set center of feedback zone
c
            xfc = xp(n)
            yfc = yp(n)
            zfc = zp(n)
            fbuff = ibuff + 2._RKIND
c
c         Check bounds - if star particle is near grid edge
c           then shift center of feedback region
c
            if (xfc .lt. xstart+fbuff*dx .or. 
     &          xfc .gt. xstart+dx*nx-fbuff*dx .or.
     &          yfc .lt. ystart+fbuff*dx .or. 
     &          yfc .gt. ystart+dx*ny-fbuff*dx .or.
     &          zfc .lt. zstart+fbuff*dx .or. 
     &          zfc .gt. zstart+dx*nz-fbuff*dx) then
c               
               xfcshift = xfc
               yfcshift = yfc 
               zfcshift	= zfc
c
               xfc = max(xfc,xstart+fbuff*dx) 
               yfc = max(yfc,ystart+fbuff*dx) 
               zfc = max(zfc,zstart+fbuff*dx) 
c              
               xfc = min(xfc,xstart+dx*nx-fbuff*dx) 
               yfc = min(yfc,ystart+dx*ny-fbuff*dx) 
               zfc = min(zfc,zstart+dx*nz-fbuff*dx)
c
               xfcshift = xfcshift - xfc
               yfcshift = yfcshift - yfc 
               zfcshift	= zfcshift - zfc
            endif
c
c         If using zeus, then velocities are face-centered so shift
c
            face_shift = 0._RKIND
            if (imethod .eq. 2) face_shift = 0.5_RKIND
c
c         Compute index of the first cell to add momentum,
c           accounting for possible face-centering
c
            xface = (xfc - xstart)/dx - 0.5_RKIND - face_shift
            yface = (yfc - ystart)/dx - 0.5_RKIND - face_shift
            zface = (zfc - zstart)/dx - 0.5_RKIND - face_shift
c
            iface = int(xface + 0.5_RKIND)
            jface = int(yface + 0.5_RKIND)
            kface = int(zface + 0.5_RKIND)
c
            dxf = real(iface) + 0.5_RKIND - xface
            dyf = real(jface) + 0.5_RKIND - yface
            dzf = real(kface) + 0.5_RKIND - zface
c
c         Compute index of the first cell to add mass, assuming cell-centering
c           This is still used for CIC algorithm even though we don't add mass
c 
            xpos = (xfc - xstart)/dx - 0.5_RKIND
            ypos = (yfc - ystart)/dx - 0.5_RKIND
            zpos = (zfc - zstart)/dx - 0.5_RKIND
c
            ic = int(xpos + 0.5_RKIND)
            jc = int(ypos + 0.5_RKIND)
            kc = int(zpos + 0.5_RKIND)
c
            dxc = real(ic) + 0.5_RKIND - xpos
            dyc = real(jc) + 0.5_RKIND - ypos
            dzc = real(kc) + 0.5_RKIND - zpos
c
c         Compute kinetic energy of affected cells before adding
c           feedback. Energy lost to momentum cancellation will be
c           added to the particle's host cell as thermal energy.
c
            do k = -1_IKIND, +2_IKIND
               do j = -1_IKIND, +2_IKIND
                  do i = -1_IKIND, +2_IKIND
                     ke_before(i+2, j+2, k+2) = 
     &                       0.5_RKIND*d(ic+i, jc+j, kc+k)*
     &                               (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                                w(ic+i ,jc+j ,kc+k)**2_IKIND)
                  enddo
               enddo
            enddo
c
c         Convert velocities to momenta and transform
c           into frame comoving with particle
c
            call momentum_emf(u, v, w, d, up(n), vp(n), wp(n),
     &                        nx, ny, nz, ic, jc, kc, 
     &                        iface, jface, kface, imethod,
     &                        +1_IKIND)
c
c         Add momentum using 3-point CIC stencil. Kinetic energy
c           will be updated after we convert back to the lab frame
c           for non-ZEUS methods.
c
            call add_feedback_emf(u, v, w, nx, ny, nz, ic, jc, kc,
     &                            iface, jface, kface, dxf, dyf, dzf,
     &                            dxc, dyc, dzc, delta_p)
c
c         Convert momenta back to velocities and transform
c           back to lab frame
c
            call momentum_emf(u, v, w, d, up(n), vp(n), wp(n),
     &                        nx, ny, nz, ic, jc, kc, 
     &                        iface, jface, kface, imethod,
     &                        -1_IKIND)
c
c         Calculate new kinetic energy. Energy lost to momentum 
c           cancellation will be added to the particle's host 
c           cell as thermal energy.
c
            ke_thermalized = 0._RKIND
            do k = -1_IKIND, +2_IKIND
               do j = -1_IKIND, +2_IKIND
                  do i = -1_IKIND, +2_IKIND
                     ke_after = 0.5_RKIND*
     &                          d(ic+i, jc+j, kc+k) *
     &                         (u(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                          v(ic+i ,jc+j ,kc+k)**2_IKIND + 
     &                          w(ic+i ,jc+j ,kc+k)**2_IKIND) 
                     delta_ke = ke_after - ke_before(i+2,j+2,k+2)
                     if (delta_ke .lt. 0) then
                        ke_thermalized =  ke_thermalized - delta_ke
                     endif
                     if (imethod .ne. 2_IKIND) then
                        te(ic+i ,jc+j ,kc+k) = 
     &                       te(ic+i ,jc+j ,kc+k) + 
     &                       delta_ke/d(ic+i, jc+j, kc+k)
                     endif
                  enddo
               enddo
            enddo
            te(ip,jp,kp) = te(ip,jp,kp) + ke_thermalized/d(ip,jp,kp)
            if (idual .eq. 1) then
               ge(ip,jp,kp) = ge(ip,jp,kp) + ke_thermalized/d(ip,jp,kp)
            endif
c
c
 10         continue
         endif
c
 100     continue
c
      enddo
c
      return
      end
c
c ==========================================================
c
c     Convert velocities to momentum and back
c
      subroutine momentum_emf(u, v, w, d, up, vp, wp, 
     &                        nx, ny, nz, ic, jc, kc, 
     &                        iface, jface, kface, imethod, idir)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      INTG_PREC iface, jface, kface, imethod, idir
      R_PREC    d(nx, ny, nz)
      R_PREC    u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz)
      R_PREC    up, vp, wp
c
c     Locals
c
      INTG_PREC i, j, k
c
c     Error check
c
      if (idir .ne. -1_IKIND .and. idir .ne. 1_IKIND) then
         write(6,*) 'incorrect idir value in momentum call'
         ERROR_MESSAGE
      endif
c
c     Loop over velocities, multiplying by densities (or dividing if
c       converting back)
c
      do k = -1_IKIND, +2_IKIND
         do j = -1_IKIND, +2_IKIND
            do i = -1_IKIND, +2_IKIND
c
c              idir = +1: convert vel -> mom
c               
               if (idir. eq. +1_IKIND) then
                  if (imethod .eq. 2_IKIND) then
                     u(iface+i ,jc+j ,kc+k) = 
     &                    (u(iface+i ,jc+j ,kc+k) - up)
     &                    * 0.5_RKIND * (d(iface+i  , jc+j, kc+k) +
     &                    d(iface+i+1_IKIND, jc+j, kc+k))
                     v(ic+i ,jface+j ,kc+k) = 
     &                    (v(ic+i ,jface+j ,kc+k) - vp)
     &                    * 0.5_RKIND * (d(ic+i, jface+j  , kc+k) + 
     &                    d(ic+i, jface+j+1_IKIND, kc+k))
                     w(ic+i ,jc+j ,kface+k) = 
     &                    (w(ic+i ,jc+j ,kface+k) - wp)
     &                    * 0.5_RKIND * (d(ic+i, jc+j, kface+k  ) + 
     &                    d(ic+i, jc+j, kface+k+1_IKIND))
                  else
                     u(ic+i ,jc+j ,kc+k) = (u(ic+i ,jc+j ,kc+k)-up) *
     &                                      d(ic+i, jc+j, kc+k)
                     v(ic+i ,jc+j ,kc+k) = (v(ic+i ,jc+j ,kc+k)-vp) *
     &                                      d(ic+i, jc+j, kc+k)
                     w(ic+i ,jc+j ,kc+k) = (w(ic+i ,jc+j ,kc+k)-wp) *
     &                                      d(ic+i, jc+j, kc+k)
                  endif
c
c              if idir = -1: convert mom -> vel
c               
               else
                  if (imethod .eq. 2_IKIND) then
                     u(iface+i ,jc+j ,kc+k) = u(iface+i ,jc+j ,kc+k)
     &                    /( 0.5_RKIND * (d(iface+i  , jc+j, kc+k) +
     &                             d(iface+i+1_IKIND, jc+j, kc+k))) + up
                     v(ic+i ,jface+j ,kc+k) = v(ic+i ,jface+j ,kc+k) 
     &                    /( 0.5_RKIND * (d(ic+i, jface+j  , kc+k) + 
     &                             d(ic+i, jface+j+1_IKIND, kc+k))) + vp
                     w(ic+i ,jc+j ,kface+k) = w(ic+i ,jc+j ,kface+k) 
     &                    /( 0.5_RKIND * (d(ic+i, jc+j, kface+k  ) + 
     &                             d(ic+i, jc+j, kface+k+1_IKIND))) + wp
                  else
                     u(ic+i ,jc+j ,kc+k) = u(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + up
                     v(ic+i ,jc+j ,kc+k) = v(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + vp
                     w(ic+i ,jc+j ,kc+k) = w(ic+i ,jc+j ,kc+k) /
     &                                     d(ic+i, jc+j, kc+k) + wp
                  endif
               endif
c
            enddo
         enddo
      enddo
c
      return
      end
c
c ==========================================================
c
c     Sum mass, momentum, and energy, and add mass and momentum, if requested
c       Note that pu, pv, pw are momenta
c
      subroutine add_feedback_emf(pu, pv, pw, nx, ny, nz, 
     &                            ic, jc, kc, iface, jface, kface,
     &                            dxf, dyf, dzf, dxc, dyc, dzc,
     &                            mom_per_cell)
c
      implicit none
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC nx, ny, nz, ic, jc, kc
      INTG_PREC iface, jface, kface
      R_PREC    pu(nx,ny,nz), pv(nx,ny,nz), pw(nx,ny,nz)
      R_PREC    mom_per_cell
      R_PREC    dxf, dyf, dzf, dxc, dyc, dzc,
     &        dxf1, dyf1, dzf1, dxc1, dyc1, dzc1
c
c     Locals
c
      INTG_PREC i, j, k, i1, j1, k1
      R_PREC    delta_pu, delta_pv, delta_pw
c
c     Loop "cells" in particle-frame
c
      do k = -1_IKIND, +1_IKIND
         do j = -1_IKIND, +1_IKIND
            do i = -1_IKIND, +1_IKIND
c
c     For each particle "cell", do CIC-like deposit
c     compute zone and face centered weight factors
c
               do i1 = i, i+1_IKIND
                  dxf1 = dxf
                  dxc1 = dxc
                  if (i1 .eq. i+1_IKIND) dxf1 = 1._RKIND - dxf
                  if (i1 .eq. i+1_IKIND) dxc1 = 1._RKIND - dxc
                  do j1 = j, j+1_IKIND
                     dyf1 = dyf
                     dyc1 = dyc
                     if (j1 .eq. j+1_IKIND) dyf1 = 1._RKIND - dyf
                     if (j1 .eq. j+1_IKIND) dyc1 = 1._RKIND - dyc
                     do k1 = k, k+1_IKIND
                        dzf1 = dzf
                        dzc1 = dzc
                        if (k1 .eq. k+1_IKIND) dzf1 = 1._RKIND - dzf
                        if (k1 .eq. k+1_IKIND) dzc1 = 1._RKIND - dzc

                        delta_pu   =  i*mom_per_cell*dxf1*dyc1*dzc1
                        delta_pv   =  j*mom_per_cell*dxc1*dyf1*dzc1
                        delta_pw   =  k*mom_per_cell*dxc1*dyc1*dzf1
c
c     Add momentum
c
                        pu(iface+i1 ,jc+j1 ,kc+k1) = 
     &                       pu(iface+i1 ,jc+j1 ,kc+k1)+ delta_pu
                        pv(ic+i1 ,jface+j1 ,kc+k1) = 
     &                       pv(ic+i1 ,jface+j1 ,kc+k1) + delta_pv
                        pw(ic+i1 ,jc+j1 ,kface+k1) = 
     &                       pw(ic+i1 ,jc+j1 ,kface+k1)+ delta_pw
c     
c     End loop over CIC-deposit
c     
                     enddo
                  enddo
               enddo
c     
c     End loop over "cells" in particle-frame
c     
            enddo
         enddo
      enddo
c
c     write(6,*) 'tot_mass:',tot_mass
c     
      return
      end