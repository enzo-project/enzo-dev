#include "fortran.def"
#include "phys_const.def"
!=======================================================================
!\////////////////////  SUBROUTINE CALC_TDUST_3D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine calc_tdust_3d(
     &     d, de, HI, HII, 
     &     HeI, HeII, HeIII,
     &     HM, H2I, H2II, 
     &     in, jn, kn, 
     &     nratec, iexpand,
     &     ispecies, idim,
     &     is, js, ks, 
     &     ie, je, ke, 
     &     aye, temstart, temend,
     &     gasgra,
     &     utem, uxyz, uaye,
     &     urho, utim,
     &     gas_temp, dust_temp)

!  COMPUTE THE DUST TEMPERATURE
!
!  written by: Britton Smith
!  date: July, 2011
!  modified1: 
!
!  PURPOSE:
!    Calculate dust heat balance to get the dust temperature.
!
!  INPUTS:
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"
!  Arguments

      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, nratec,
     &        iexpand, ispecies, idim
      R_PREC    aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim
      R_PREC    d(in,jn,kn),
     &     de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &     HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &     HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &     gas_temp(in,jn,kn), dust_temp(in,jn,kn)
      R_PREC    gasgra(nratec)

!  Parameters

      real*8 mh
      parameter (mh = mass_h)      !DPC

!  Locals

      INTG_PREC i, j, k
      R_PREC trad, zr, logtem0, logtem9, dlogtem
      real*8 coolunit, dbase1, tbase1, xbase1

!  Slice locals
 
      INTG_PREC indixe(in)
      R_PREC t1(in), t2(in), logtem(in), tdef(in), 
     &     tgas(in), tdust(in), nh(in), gasgr(in)
          
!  Iteration mask for multi_cool

      LOGIC_PREC itmask(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set log values of start and end of lookup tables

      logtem0  = log(temstart)
      logtem9  = log(temend)
      dlogtem  = (log(temend) - log(temstart))/real(nratec-1)

!     Set units

      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      zr       = 1._RKIND/(aye*uaye) - 1._RKIND

!     Set compton cooling coefficients (and temperature)

      if (iexpand .eq. 1) then
         trad = 2.73_RKIND * (1._RKIND + zr)
      else
         trad = tiny
      endif

!     Loop over slices (in the k-direction)

      do k = ks+1, ke+1
         do j = js+1, je+1

            do i = is+1, ie+1

!     Set itmask to all true

               itmask(i) = .true.

!     Compute hydrogen number density

               nh(i) = HI(i,j,k) + HII(i,j,k)
               if (ispecies .gt. 1) then
                  nh(i) = nh(i) + H2I(i,j,k) + H2II(i,j,k)
               endif

!     We have not converted to proper, so use urho and not dom

               nh(i) = nh(i) * urho / mh

!     Compute log temperature and truncate if above/below table max/min

               tgas(i)   = gas_temp(i,j,k)
               logtem(i) = log(tgas(i))
               logtem(i) = max(logtem(i), logtem0)
               logtem(i) = min(logtem(i), logtem9)

!     Compute index into the table and precompute parts of linear interp

               indixe(i) = min(nratec-1,
     &              max(1,int((logtem(i)-logtem0)/dlogtem,IKIND)+1))
               t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
               t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
               tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

!     Lookup values and do a linear temperature in log(T)
!     Convert back to cgs

               gasgr(i) = gasgra(indixe(i)) + tdef(i)
     &              *(gasgra(indixe(i)+1) -gasgra(indixe(i)))
               gasgr(i) = gasgr(i) * coolunit / mh

            enddo

!     --- Compute dust temperature in a slice ---

            call calc_tdust(tdust, tgas, nh, gasgr, 
     &           itmask, trad, in, is, ie, j, k)

!     Copy slice values back to grid

            do i = is+1, ie+1
               dust_temp(i,j,k) = tdust(i)
            enddo

         enddo
      enddo

      return
      end

