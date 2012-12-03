#include "fortran.def"
#include "phys_const.def"
#include "error.def"

!=======================================================================
!///////////////////////  SUBROUTINE CALC_TDUST  \\\\\\\\\\\\\\\\\\\\\\\

      subroutine calc_tdust(
     &     tdust, tgas, nh, gasgr, itmask, 
     &     trad, in, is, ie, j, k)

!  CALCULATE EQUILIBRIUM DUST TEMPERATURE
!
!  written by: Britton Smith
!  date:       February, 2011
!  modified1: 
!
!  PURPOSE:
!    Calculate dust temperature.
!
!  INPUTS:
!     in       - dimension of 1D slice
!
!     tdust    - dust temperature
!
!     tgas     - gas temperature
!     nh       - H number density
!     gasgr    - gas/grain heat transfer rate
!
!     trad     - CMB temperature
!
!     is,ie    - start and end indices of active region (zero based)
!     j,k      - indices of 1D slice
!
!     itmask   - iteration mask
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!  Arguments

      INTG_PREC in, is, ie, j, k

      R_PREC tdust(in), tgas(in), nh(in), gasgr(in)
      R_PREC trad

!  Iteration mask

      LOGIC_PREC itmask(in)

!  Parameters

      R_PREC radf, t_subl
      parameter(radf = 4._RKIND * sigma_sb)
      parameter(t_subl = 1.5e3_RKIND) ! grain sublimation temperature
      R_PREC tol, bi_tol, minpert
      parameter(tol = 1.e-5_RKIND, bi_tol = 1.e-3_RKIND, 
     &     minpert = 1.e-10_RKIND)
      INTG_PREC itmax, bi_itmax
      parameter(itmax = 50, bi_itmax = 30)

!  Locals

      INTG_PREC i, iter, c_done, c_total, nm_done

      real*8 pert_i, trad4

!  Slice Locals

      real*8 kgr(in), kgrplus(in), sol(in), solplus(in), 
     &     slope(in), tdplus(in), tdustnow(in), tdustold(in), pert(in),
     &     bi_t_mid(in), bi_t_high(in)
      LOGIC_PREC nm_itmask(in), bi_itmask(in)

      real*8 sola, solb, solpa, solpb

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      pert_i = 1.e-3_RKIND

      trad  = max(1._RKIND, trad)
      trad4 = trad**4

!     Set total cells for calculation

      c_done = 0
      nm_done = 0
      c_total = ie - is + 1

!     Set local iteration mask and initial guess

      do i = is+1, ie+1
         nm_itmask(i) = itmask(i)
         bi_itmask(i) = itmask(i)
         if ( nm_itmask(i) ) then

            if (trad .ge. tgas(i)) then
               tdustnow(i) = trad
               nm_itmask(i) = .false.
               bi_itmask(i) = .false.
               c_done = c_done + 1
               nm_done = nm_done + 1
            else if (tgas(i) .gt. t_subl) then
!     Use bisection if T_gas > grain sublimation temperature.
               nm_itmask(i) = .false.
               nm_done = nm_done + 1
            else
               tdustnow(i) = 0.5_RKIND * (trad + tgas(i))
               pert(i) = pert_i
            endif

         else
            c_done = c_done + 1
            nm_done = nm_done + 1
         endif
      enddo

!     Iterate to convergence with Newton's method

      do iter = 1, itmax

!     Loop over slice

         do i = is+1, ie+1
            if ( nm_itmask(i) ) then

               tdplus(i) = max(1.e-3_RKIND, ((1._RKIND+pert(i)) 
     &              * tdustnow(i)))

            endif
         enddo

!     Calculate grain opacities

         call calc_kappa_gr(tdustnow, kgr, nm_itmask, 
     &        in, is, ie, t_subl)

         call calc_kappa_gr(tdplus, kgrplus, nm_itmask, 
     &        in, is, ie, t_subl)

         do i = is+1, ie+1
            if ( nm_itmask(i) ) then

!     Use Newton's method to solve for Tdust

               sol(i) = radf * kgr(i) * 
     &              (trad4 - tdustnow(i)**4) +
     &              (gasgr(i) * nh(i) * 
     &              (tgas(i) - tdustnow(i)))

               sola = radf * kgr(i) * 
     &              (trad4 - tdustnow(i)**4)
               solb = (gasgr(i) * nh(i) * 
     &              (tgas(i) - tdustnow(i)))

               solplus(i) = radf * kgrplus(i) * 
     &              (trad4 - tdplus(i)**4) +
     &              (gasgr(i) * nh(i) * 
     &              (tgas(i) - tdplus(i)))

               solpa = radf * kgrplus(i) * 
     &              (trad4 - tdplus(i)**4)
               solpb = (gasgr(i) * nh(i) * 
     &              (tgas(i) - tdplus(i)))

               slope(i) = (solplus(i) - sol(i)) / 
     &              (pert(i) * tdustnow(i))

               tdustold(i) = tdustnow(i)
               tdustnow(i) = tdustnow(i) - (sol(i) / slope(i))

               pert(i) = max(min(pert(i), 
     &              (0.5_RKIND * abs(tdustnow(i) - tdustold(i)) / 
     &              tdustnow(i))), minpert)

!     If negative solution calculated, give up and wait for bisection step.
               if (tdustnow(i) .lt. trad) then
                  nm_itmask(i) = .false.
                  nm_done = nm_done + 1
!     Check for convergence of solution
               else if (abs(sol(i) / solplus(i)) .lt. tol) then
                  nm_itmask(i) = .false.
                  c_done = c_done + 1
                  bi_itmask(i) = .false.
                  nm_done = nm_done + 1
               endif

!     if ( nm_itmask(i) )
            endif

!     End loop over slice
         enddo

!     Check for all cells converged
         if (c_done .ge. c_total) go to 666

!     Check for all cells done with Newton method
!     This includes attempts where a negative solution was found
         if (nm_done .ge. c_total) go to 555

!     End iteration loop for Newton's method
      enddo

 555  continue

!     If iteration count exceeded, try once more with bisection
      if (c_done .lt. c_total) then
         do i = is+1, ie+1
            if ( bi_itmask(i) ) then
               tdustnow(i)  = trad
               bi_t_high(i) = tgas(i)
            endif
         enddo

         do iter = 1, bi_itmax

            do i = is+1, ie+1
               if ( bi_itmask(i) ) then

                  bi_t_mid(i) = 0.5_RKIND * (tdustnow(i) + bi_t_high(i))
                  if (iter .eq. 1) then
                     bi_t_mid(i) = min(bi_t_mid(i), t_subl)
                  endif

                  call calc_kappa_gr(bi_t_mid, kgr, bi_itmask, 
     &                 in, is, ie, t_subl)
                  sol(i) = radf * kgr(i) * 
     &                 (trad4 - bi_t_mid(i)**4) +
     &                 (gasgr(i) * nh(i) * 
     &                 (tgas(i) - bi_t_mid(i)))
                  if (sol(i) .gt. 0._RKIND) then
                     tdustnow(i) = bi_t_mid(i)
                  else
                     bi_t_high(i) = bi_t_mid(i)
                  endif

                  if ((abs(bi_t_high(i) - tdustnow(i)) / tdustnow(i)) 
     &                 .le. bi_tol) then
                     bi_itmask(i) = .false.
                     c_done = c_done + 1
                  endif

!     Check for all cells converged
                  if (c_done .ge. c_total) go to 666

!     if ( nm_itmask(i) )
               endif

!     End loop over slice
            enddo

!     End iteration loop for bisection
         enddo

!     If iteration count exceeded with bisection, end of the line.
         if (iter .gt. itmax) then
            write(6,*) 'CALC_TDUST failed using both methods for ',
     &           (c_total - c_done), 'cells.'
            ERROR_MESSAGE
         endif

!     if (iter .gt. itmax) then
      endif

 666  continue

!     Copy values back to thrown slice
      do i = is+1, ie+1
         if ( itmask(i) ) then

!     Check for bad solutions
            if (tdustnow(i) .lt. 0._RKIND) then
               write(6, *) 'CALC_TDUST Newton method - ',
     &              'T_dust < 0: i = ', i, 'j = ', j,
     &              'k = ', k, 'nh = ', nh(i), 
     &              't_gas = ', tgas(i), 't_rad = ', trad,
     &              't_dust = ', tdustnow(i)
               ERROR_MESSAGE
            endif

            tdust(i) = tdustnow(i)
         endif
      enddo

      return
      end

!=======================================================================
!///////////////////////  SUBROUTINE CALC_TDUST  \\\\\\\\\\\\\\\\\\\\\\\

      subroutine calc_kappa_gr(
     &     tdust, kgr, itmask, in, is, ie, t_subl)

!  CALCULATE GRAIN PLANK MEAN OPACITY
!
!  written by: Britton Smith
!  date:       September, 2011
!  modified1: 
!
!  PURPOSE:
!    Calculate grain plank mean opacity
!
!  INPUTS:
!     in       - i dimension of 3D fields
!
!     tdust    - dust temperature
!
!     is,ie    - start and end indices of active region (zero based)
!
!     itmask   - iteration mask
!
!     t_subl   - grain sublimation temperature
!
!  OUTPUTS:
!     kgr      - opacities
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!  Arguments

      INTG_PREC in, is, ie
      R_PREC t_subl
      real*8 tdust(in)

!  Iteration mask

      LOGIC_PREC itmask(in)

!  Parameters

      R_PREC kgr1, kgr200
      parameter(kgr1 = 4.0e-4_RKIND, kgr200 = 16.0_RKIND)

!  Locals

      INTG_PREC i

!  Slice Locals

      real*8 kgr(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

      do i = is+1, ie+1
         if ( itmask(i) ) then

!     Temperature dependence from Dopcke et al. (2011).
!     Normalized to Omukai (2000).

            if (tdust(i) .lt. 200._RKIND) then
               kgr(i) = kgr1 * tdust(i)**2
            else if (tdust(i) .lt. t_subl) then
               kgr(i) = kgr200
            else
               kgr(i) = max(tiny, 
     &              (kgr200 * (tdust(i) / 1.5e3_RKIND)**(-12)))
            endif

         endif
      enddo

      return
      end
