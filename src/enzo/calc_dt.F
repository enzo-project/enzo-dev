#include "fortran.def"
c=======================================================================
c////////////////////////  SUBROUTINE CALC_DT  \\\\\\\\\\\\\\\\\\\\\\\\\
c
c
      subroutine calc_dt(rank, idim, jdim, kdim, 
     &                   i1, i2, j1, j2, k1, k2, ihydro, C2,
     &                   dx, dy, dz, vgx, vgy, vgz, gamma, ipfree, aye,
     &                   d, p, u, v, w, dt, dtviscous)
#ifndef CONFIG_PFLOAT_16
c
c  COMPUTES TIME STEP FOR NEXT CYCLE
c
c     written by: Greg Bryan
c     date:       February, 1996
c     modified1:  Alexei Kritsuk, Jan. 2001; changed the expression for
c                 hydro time step for PPM_DE. Now it follows the linear 
c                 stability condition for multidimensional Godunov scheme 
c                 (Godunov 1959). It is now safe to use CourantSafetyNumber =
c                 0.8 for PPM_DE runs.
c
c  PURPOSE:  Computes the new timestep using the Courant condition.
c            (For rank < 3, the unused fields and cell widths may be
c             null)
c
c  INPUTS:
c    rank    - rank of fields
c    i,j,dim - declared dimensions of fields
c    i,j,k1  - start index of active region in fields (0 based)
c    i,j,k2  - end index of active region in fields (0 based)
c    ihydro  - Hydro method (2 - Zeus), used for viscosity computation
c    C2      - coefficient of quadratic artificial viscosity
c    dx,y,z  - cell widths along each dimension
c    vgx,y,z - grid bulk velocity
c    gamma   - ratio of specific heats
c    ipfree  - pressure free flag (1 = on, 0 = off)
c    aye     - expansion factor (or 1 if not using comoving coordinates)
c    d,p     - density and pressure fields
c    u,v,w   - velocity fields (x,y,z)
c    dtviscous - viscous time for stability (if used)
c
c  OUTPUTS:
c    dt      - minimum allowed dt (without Courant safety factor)
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c     Arguments
c
      INTG_PREC idim, jdim, kdim, i1, i2, j1, j2, k1, k2, rank, ipfree,
     &        ihydro
      P_PREC dx(idim), dy(jdim), dz(kdim)
      R_PREC    dt, vgx, vgy, vgz, gamma, aye, C2, dtviscous
      R_PREC    d(idim,jdim,kdim), p(idim,jdim,kdim), u(idim,jdim,kdim),
     &        v(idim,jdim,kdim), w(idim,jdim,kdim)
c
c     Locals
c
      INTG_PREC i,j,k
      R_PREC    cs, dt1
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c     Set initial timestep to a large number
c
      dt = huge
c
c     one-dimensional version
c
      if (rank .eq. 1) then
         do i = i1+1, i2+1
            cs = max(sqrt(gamma*p(i,1,1)/d(i,1,1)), tiny)
            if (ipfree .eq. 1) cs = tiny
            dt = min(dt, REAL(dx(i)*aye/(cs + abs(u(i,1,1)-vgx)),RKIND))
         enddo
         if (ihydro .eq. 2) then
            do i = i1+1, i2+1
               dtviscous = min(dtviscous,
     &              REAL(dx(i)*aye/(4._RKIND*C2*
     &              max(-u(i+1,1,1)+u(i,1,1),tiny)),
     &              RKIND))
            enddo
         endif
      endif
c
c     two-dimensional version
c
      if (rank .eq. 2) then
         do j = j1+1, j2+1
            do i = i1+1, i2+1
               cs = max(sqrt(gamma*p(i,j,1)/d(i,j,1)), tiny)
               if (ipfree .eq. 1) cs = tiny
               if (ihydro .eq. 0) then
c
c 	Here is Godunov's formula (to make sure ppm works with 0.8)
c
 	       dt = min(dt, REAL(aye/((cs + abs(u(i,j,1)-vgx))/dx(i)+
     &                       (cs + abs(v(i,j,1)-vgy))/dy(j)), RKIND))
               else
c
c       The way it was originally in enzo
c
               dt = min(dt, REAL(dx(i)*aye/(cs + abs(u(i,j,1)-vgx)),
     &                           RKIND),
     &                      REAL(dy(j)*aye/(cs + abs(v(i,j,1)-vgy)),
     &                           RKIND))
               endif
            enddo
c
            if (ihydro .eq. 2) then
               do i = i1+1, i2+1
                  dtviscous = min(dtviscous, 
     &                 REAL(dx(i)*aye/(4._RKIND*C2*max(
     &                 -u(i+1,j,1)+u(i,j,1),tiny)), RKIND),
     &                 REAL(dy(j)*aye/(4._RKIND*C2*max(
     &                 -v(i,j+1,1)+v(i,j,1),tiny)), RKIND))
               enddo
            endif
         enddo
      endif
c
c     three-dimensional version
c
      if (rank .eq. 3) then
         do k = k1+1, k2+1
            do j = j1+1, j2+1
               do i = i1+1, i2+1
                  if (d(i,j,k) .ne. d(i,j,k) .or.
     &                p(i,j,k) .ne. p(i,j,k))
     &               write(6,*) 'calc_dt',d(i,j,k),p(i,j,k),i,j,k
                  cs = max(sqrt(gamma*p(i,j,k)/d(i,j,k)), tiny)
                  if (ipfree .eq. 1) cs = tiny
                  if (ihydro .eq. 0) then
c     
c     Godunov's formula.
c     
                     dt1 = aye/((cs + abs(u(i,j,k)-vgx))/dx(i) +
     &                    (cs + abs(v(i,j,k)-vgy))/dy(j) +
     &                    (cs + abs(w(i,j,k)-vgz))/dz(k))
                  else
c     
c     The way it was originally in enzo
c     
                     dt1 = min(dx(i)*aye/(cs + abs(u(i,j,k)-vgx)),
     &                    dy(j)*aye/(cs + abs(v(i,j,k)-vgy)),
     &                    dz(k)*aye/(cs + abs(w(i,j,k)-vgz)))
                  endif
                  dt = min(dt, dt1)
C                   if (dt1 .lt. 1.0e-5) then
C                      write(6,*) d(i-1,j,k), d(i+1,j,k), d(i,j-1,k),
C      $                    d(i,j+1,k), d(i,j,k-1), d(i,j,k+1)
C                      write(6,1000) dt1,d(i,j,k),
C      &                    p(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k)
C                   end if
 1000             format('calc_dt (small dt): dt,d,p,uvw=',1p,6e12.3)
               enddo
c
               if (ihydro .eq. 2) then
                  do i = i1+1, i2+1
	            dt1 = min(
     &                    dx(i)*aye/(4._RKIND*C2*max(
     &                        -u(i+1,j,k)+u(i,j,k), tiny)),
     &                    dy(j)*aye/(4._RKIND*C2*max(
     &                        -v(i,j+1,k)+v(i,j,k), tiny)),
     &                    dz(k)*aye/(4._RKIND*C2*max(
     &                        -w(i,j,k+1)+w(i,j,k), tiny)))
                     dtviscous = min(dtviscous, dt1)
                  enddo
               endif
            enddo
         enddo
      endif
c
      return
#endif
      end
