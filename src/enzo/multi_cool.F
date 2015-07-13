#include "fortran.def"
#include "phys_const.def"
#include "error.def"

!=======================================================================
!//////////////////////  SUBROUTINE MULTI_COOL  \\\\\\\\\\\\\\\\\\\\\\\\


      subroutine multi_cool(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, imethod, igammah,
     &                idual, ispecies, imetal, imcool, idust, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht,
     &                dt, aye, temstart, temend,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma, z_solar,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
     &                gpldla, gphdla, hdltea, hdlowa,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                gasgra, metala, n_xe, xe_start, xe_end,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma, 
     &                ih2optical, iciecool, ciecoa, 
     &                icmbTfloor, iClHeat,
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating)


!  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
!  modified2: October, 1996 by GB; moved to AMR
!  modified3: February, 2003 by Robert Harkness; iteration mask
!  modified4: November, 2003 by Robert Harkness; tighten convergence
!
!  PURPOSE:
!    Solve the energy cooling equations.
!
!  INPUTS:
!    is,ie   - start and end indicies of active region (zero-based!)
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "fortran_types.def"

!  Arguments

      INTG_PREC in, jn, kn, is, js, ks, ie, je, ke, nratec, imethod,
     &        idual, iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        iradtype, nfreq, imetalregen, iradshield, n_xe, iradtrans,
     &        imcool, idust, igammah
      R_PREC    dt, aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim, z_solar,
     &        eta1, eta2, gamma, xe_start, xe_end
      R_PREC    d(in,jn,kn),   ge(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      R_PREC    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)
      R_PREC    photogamma(in,jn,kn)
      INTG_PREC ih2optical, iciecool
      R_PREC    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec)
      R_PREC    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec), gasgra(nratec),
     &        ciecoa(nratec)
      R_PREC    metala(nratec, n_xe)
      R_PREC    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      R_PREC    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        inutot(nfreq), avgsighp, avgsighep, avgsighe2p
      R_PREC    gammaha

!  Cloudy cooling data

      INTG_PREC icmbTfloor, iClHeat, clGridRank, clDataSize
      INTG_PREC clGridDim(clGridRank)
      R_PREC clEleFra
      R_PREC clPar1(clGridDim(1)), clPar2(clGridDim(2)), 
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)), 
     &     clPar5(clGridDim(5))
      R_PREC clCooling(clDataSize), clHeating(clDataSize)

!  Parameters

      INTG_PREC itmax, ijk
      parameter (itmax = 10000, ijk = MAX_ANY_SINGLE_DIRECTION)

#ifdef CONFIG_BFLOAT_4
      R_PREC tolerance
      parameter (tolerance = 1.e-5_RKIND)
#endif

#ifdef CONFIG_BFLOAT_8
      R_PREC tolerance
      parameter (tolerance = 1.e-10_RKIND)
#endif

      real*8 mh
      parameter (mh = mass_h)

!  Locals

      INTG_PREC i, j, k, n, iter
      R_PREC dom, energy
      R_PREC dt2, ttmin, comp1, comp2

!  Row locals
 
      INTG_PREC indixe(ijk)
      R_PREC t1(ijk), t2(ijk), logtem(ijk), tdef(ijk), p2d(ijk),
     &     dtit(ijk), ttot(ijk), tgas(ijk), tgasold(ijk),
     &     tdust(ijk), metallicity(ijk), rhoH(ijk)
      real*8 edot(ijk)

!  Cooling/heating row locals

      real*8 ceHI(ijk), ceHeI(ijk), ceHeII(ijk),
     &     ciHI(ijk), ciHeI(ijk), ciHeIS(ijk), ciHeII(ijk),
     &     reHII(ijk), reHeII1(ijk), reHeII2(ijk), reHeIII(ijk),
     &     brem(ijk), cieco(ijk)
      R_PREC hyd01k(ijk), h2k01(ijk), vibh(ijk), roth(ijk), rotl(ijk),
     &     gpldl(ijk), gphdl(ijk), hdlte(ijk), hdlow(ijk)

!  Iteration mask

      LOGIC_PREC itmask(ijk)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set units

      dom      = urho*(aye**3)/mh

!     Convert densities from comoving to proper

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)/aye**3
               de(i,j,k)    = de(i,j,k)/aye**3
               HI(i,j,k)    = HI(i,j,k)/aye**3
               HII(i,j,k)   = HII(i,j,k)/aye**3
               HeI(i,j,k)   = HeI(i,j,k)/aye**3
               HeII(i,j,k)  = HeII(i,j,k)/aye**3
               HeIII(i,j,k) = HeIII(i,j,k)/aye**3
            enddo
            if (ispecies > 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)/aye**3
                  H2I(i,j,k)  = H2I(i,j,k)/aye**3
                  H2II(i,j,k) = H2II(i,j,k)/aye**3
               enddo
            endif
            if (ispecies > 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)/aye**3
                  DII(i,j,k) = DII(i,j,k)/aye**3
                  HDI(i,j,k) = HDI(i,j,k)/aye**3
               enddo
            endif
            if (imetal == 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)/aye**3
               enddo
            endif
         enddo
      enddo
 
!     Solve energy cooling (subcycle)

!     Loop over rows of cells

      do k = ks+1, ke+1
       do j = js+1, je+1

!       tolerance = 1.0e-06 * dt

        do i = is+1, ie+1
           itmask(i) = .true.
        end do

!       Set time elapsed to zero for each cell

        do i = is+1, ie+1
           ttot(i) = 0._RKIND
        enddo

!       Loop over cooling subcycles
     
        do iter = 1, itmax

!       Compute the cooling rate on this row

          call cool1d_multi(
     &                d, e, ge, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, idual, imethod,
     &                iexpand, ispecies, imetal, imcool, idust, idim,
     &                is, ie, j, k, ih2co, ipiht, iter, igammah,
     &                aye, temstart, temend, z_solar,
     &                utem, uxyz, uaye, urho, utim,
     &                eta1, eta2, gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                gasgra, metala, n_xe, xe_start, xe_end,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d, tdust, metallicity, rhoH,
     &                inutot, iradtype, nfreq, imetalregen,
     &                iradshield, avgsighp, avgsighep, avgsighe2p,
     &                iradtrans, photogamma, 
     &                ih2optical, iciecool, ciecoa, cieco,
     &                icmbTfloor, iClHeat,
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating,
     &                itmask
     &                     )

!         Compute maximum timstep that keeps any fractional change to 10%
!         (maximum timestep is 1/2 hydro step).
!         Then, use each cell~s individual timestep tp update it~s energy
!         Find minimum elapsed timestep (ttot)

          dt2 = dt/2._RKIND
          ttmin = huge
	  write(0,*) 'is=',is, ie
	  write(81,*) 'is=',is, ie
          do i = is+1, ie+1

!            Set energy of this cell (the gamma used here is the right
!            one even for H2 since p2d is calculated with this gamma).

             if (tgas(i) <= temstart .and. edot(i) < 0.0) 
     &              edot(i) = tiny*1.e-3_RKIND

             if (abs(edot(i)) < tiny) edot(i) = tiny

             energy = max(p2d(i)/(gamma-1._RKIND), tiny)

c            energy = max(ge(i,j,k)*d(i,j,k), p2d(i)/(gamma-1._RKIND), 
c    &                    tiny)
c            if (energy < tiny) energy = d(i,j,k)*(e(i,j,k) - 
c    &              0.5_RKIND*(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))
c            energy = p(i,j,k)/(gamma-1._RKIND)


!            Compute timestep for 10% change

c            if (iter > 100) then
c               dtit(i) = min(REAL(abs(0.01_RKIND*energy/edot(i)),RKIND), 
c    &                        dt-ttot(i), dt2)
c            else
c               dtit(i) = min(REAL(abs(0.1_RKIND*energy/edot(i)),RKIND), 
c    &                        dt-ttot(i), dt2)
c            endif

             dtit(i) = min(REAL(abs(0.1_RKIND*energy/edot(i)),RKIND),
     &                     dt-ttot(i), dt2)

             if ( dt-ttot(i) <= tolerance*dt ) then
                itmask(i) = .false.
             end if

             if ( itmask(i) ) then

!            call open_mpi_error_file( 'F30', 30, 'unknown' )
!            if (iter == 1 .and. j == 4 .and. k == 4)
!    &       write(30,1000) i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
!    &       energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k),dtit(i)
!            call close_mpi_error_file( 30 )

#define NO_FORTRAN_DEBUG
#ifdef FORTRAN_DEBUG
             if (ge(i,j,k) <= 0._RKIND .and. idual == 1)
     &         write(6,*) 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),iter

             if (e(i,j,k) <= 0._RKIND .or. e(i,j,k) 
     &            .ne. e(i+ih2co-1,j,k))
     &         write(6,*) 'mc:e',ge(i,j,k),energy,d(i,j,k),e(i,j,k),
     &                        iter, edot(i), tgas(i)

             if (edot(i) .ne. edot(i+ih2co-1) .or. 
     &           d(i,j,k) .ne. d(i+ih2co-1,j,k) .or.
     &           e(i,j,k) .ne. e(i+ih2co-1,j,k) .or.
     &           dtit(i) .ne. dtit(i+ih2co-1) .or. 
     &           d(i,j,k) == 0._RKIND) 
     &           then 
                write(6,*) 'multi_cool_a:',
     &          edot(i),d(i,j,k),dtit(i),i,j,k,iter,p2d(i),
     &          energy,dom,tgas(i),e(i,j,k),energy,
     &          t1(i),t2(i),tdef(i),indixe(i),logtem(i)
                write(6,*) '2',de(i,j,k),HI(i,j,k),HII(i,j,k),
     &              HeI(i,j,k),HeII(i,j,k),HeIII(i,j,k)
                write(6,*) '3',ceHI(i),ceHeI(i),ceHeII(i),ciHI(i),
     &              ciHeI(i),ciHeII(i),ciHeIS(i),reHII(i),
     &              reHeII1(i),reHeII2(i),reHeIII(i),comp1,comp2,
     &              brem(i)
              if (ispecies > 1)
     &          write(6,*) 'H2:',rotl(i),roth(i),vibh(i),HM(i,j,k),
     &              H2I(i,j,k),H2II(i,j,k),h2k01(i),hyd01k(i)
                write(0,*) 'FATAL error (1) in MULTI_COOL'
                ERROR_MESSAGE
             endif

             if (ge(i,j,k) .ne. ge(i,j,k) .or.
     &            e(i,j,k) .ne.  e(i,j,k)) write(6,*)
     &           'multi_cool_b:',ge(i,j,k),i,j,k,iter,tgas(i),e(i,j,k)

             if (idual == 1 .and.
     &           ge(i,j,k)+edot(i)/d(i,j,k)*dtit(i) <= 0._RKIND)
     &         write(6,*) i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
     &              energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k)
#endif /* FORTRAN_DEBUG */

             if (((dtit(i)/dt < 1.e-2_RKIND) .and. 
     &            (iter > 1000) .and.
     &            (abs((dt-ttot(i))/dt) > 1.e-3_RKIND)) .or. 
     &            (iter > 1800)) then
!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!               call open_mpi_error_file( 'F3', 3, 'unknown' )
               write(6,1000) i,j,k,iter,ge(i,j,k),edot(i),tgas(i),
     &         energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k),dtit(i)
!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!               call close_mpi_error_file( 3 )
             end if

 1000        format(4(i4,1x),1p,10(e14.3))


!            Update total and gas energy

             e(i,j,k)  = e(i,j,k) + edot(i)/d(i,j,k)*dtit(i)

             if (e(i,j,k) < 0.0) then
!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!                call open_mpi_error_file( 'F3', 3, 'unknown' )
                write(6,*) 'Eijk < 0', e(i,j,k)
!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!                call close_mpi_error_file( 3 )
             end if

             if (idual == 1) then
                ge(i,j,k) = ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i)

c               ge(i,j,k) = max(ge(i,j,k) + edot(i)/d(i,j,k)*dtit(i),
c     &                      0.5*ge(i,j,k))
c            if (ge(i,j,k) <= tiny) ge(i,j,k) = (energy + 
c     &           edot(i)*dtit(i))/d(i,j,k)

                if (ge(i,j,k) <= 0._RKIND) then
!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!                   call open_mpi_error_file( 'F3', 3, 'unknown' )
                   write(6,*) 'a',ge(i,j,k),energy,d(i,j,k),e(i,j,k),
     &                        iter
!!! JB 2008-05-20 commented to avoid portability issues with Fortran MPI
!!!                call close_mpi_error_file( 3 )
                end if
             endif

!            Update time-stepping

             ttot(i) = ttot(i) + dtit(i)
             ttmin = min(ttot(i), ttmin)

             end if  ! test of itmask(i)

          enddo  ! end of loop over cells(i)

!         call open_mpi_error_file( 'F37', 37, 'unknown' )
!         write(37,'(i4,128l1)') iter, (itmask(i),i=is+1,ie+1)
!         call close_mpi_error_file( 37 )

!         If the all cells are done then skip out of loop

          if (abs(dt-ttmin) < tolerance*dt) go to 8888

         enddo  ! end of iteration loop

 8888    continue
      
!       Abort if iteration count exceeds maximum
	  write(0,*) 'is=',is, ie
         if (iter > itmax) then
	    write(0,*) 'inside if statement is=',is, ie
            write(6,*) 'MULTI_COOL iter > ',itmax,' at j,k =',j,k
            write(0,*) 'FATAL error (2) in MULTI_COOL'
            write(0,'(" dt = ",1pe10.3," ttmin = ",1pe10.3)') dt, ttmin
            write(0,'((16(1pe8.1)))') (dtit(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (ttot(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (edot(i),i=is+1,ie+1)
            ERROR_MESSAGE
         endif

         if (iter > itmax/2) then
            write(6,*) 'MULTI_COOL iter,j,k =',iter,j,k
         end if

!        call open_mpi_error_file( 'F30', 30, 'unknown' )
!        write(30,'("J=",i4,4x,"K=",i4)') j,k
!        write(30,'((12(1pd10.3)))') (ge(i,j,k),i=is+1,ie+1)
!        call close_mpi_error_file( 30 )


!      Next j,k row

       enddo
      enddo

!     Convert densities back to comoving from proper

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*aye**3
               de(i,j,k)    = de(i,j,k)*aye**3
               HI(i,j,k)    = HI(i,j,k)*aye**3
               HII(i,j,k)   = HII(i,j,k)*aye**3
               HeI(i,j,k)   = HeI(i,j,k)*aye**3
               HeII(i,j,k)  = HeII(i,j,k)*aye**3
               HeIII(i,j,k) = HeIII(i,j,k)*aye**3
            enddo
            if (ispecies > 1) then
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*aye**3
                  H2I(i,j,k)  = H2I(i,j,k)*aye**3
                  H2II(i,j,k) = H2II(i,j,k)*aye**3
               enddo
            endif
            if (ispecies > 2) then
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*aye**3
                  DII(i,j,k) = DII(i,j,k)*aye**3
                  HDI(i,j,k) = HDI(i,j,k)*aye**3
               enddo
            endif
            if (imetal == 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*aye**3
               enddo
            endif
         enddo
      enddo

      return
      end
