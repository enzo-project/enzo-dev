#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c///////////////////  SUBROUTINE CALC_PHOTO_RATES  \\\\\\\\\\\\\\\\\\\\\
c
      subroutine calc_photo_rates(NFREQ, FREQDEL, iradshield, aye,
     &                    SIGH, SIGHE, SIGHE2, INUTOT,
     &                    PHTH, PHTHE2, PHTHE, EXRY, TXRY,
     &                    PHTLAMH, PHTLAMHE, PHTLAMHE2,
     &                    AVGSIGH, AVGSIGHE, AVGSIGHE2,
     &                    AVGSIGHP, AVGSIGHEP, AVGSIGHE2P,
     &                    utim, uxyz, urho, uaye)
c
c  CALCULATES THE PHOTO-HEATING AND IONIZATION RATES GIVEN THE
c     RADIATION FIELD
c
c  written by: Renyue Cen
c  date:       
c  modified1:  September, 1999 by Greg Bryan; converted to AMR
c
c  PURPOSE:
c    Given a radiation field, this routine calcultes the photo-
c      ionization rates and photo-heating rates (including
c      Compton X-ray).
c
c  INPUTS:
c    NFREQ    - Number of frequency bins
c    FREQDEL  - space between frequency bins, in log10(eV)
c    SIGH     - HI photo-ionization heating cross-section
c    SIGHE    - HeI photo-ionization heating cross-section
c    SIGHE2   - HeII photo-ionization heating cross-section
c    iradshield - integer flag indicating if approximate radiative
c                 shielding should be used (0 - no, 1 - yes)
c    INUTOT   - total radiation intensity
c    AYE      - expansion factor in code units
c
c  OUTPUTS:
c    PHTH     - HI photoionization rate, in s^-1
c    PHTHE    - HeI photoionization rate
c    PHTHE2   - HeII photoionization rate
c    PHTLAMTH     - HI photo-heating rate, in 1e-30 erg/s
c    PHTLAMTHE    - HeI photo-heating rate
c    PHTLAMTHE2   - HeII photo-heating rate
c    EXRY     - X-ray backgroun radiation energy density, in erg/cm^3
c    TXRY     - X-ray background temperature, in K
c    AVGSIGH  - intensity weighted average cross section
c
c  PARAMETERS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
c
c  Arguments
c
      integer NFREQ, iradshield
      real    FREQDEL, SIGH(NFREQ), SIGHE(NFREQ), SIGHE2(NFREQ),
     &        PHTH, PHTHE, PHTHE2, EXRY, TXRY, INUTOT(NFREQ),
     &        PHTLAMH, PHTLAMHE, PHTLAMHE2,
     &        AVGSIGH, AVGSIGHE, AVGSIGHE2,
     &        AVGSIGHP, AVGSIGHEP, AVGSIGHE2P,
     &        utim, uxyz, urho, uaye, aye
c  Parameters
c
      real    EV2HZ, PI
c
c  Locals
c
      integer N11
      real    FNUDEL, AVGNU, FNUMNUH0, FNUMNUHE0, FNUMNUHE20,
     &        QNICK, CROSSS, FJXRY11, FJXRY21,
     &        FJXRY12, FJXRY22, FJXRY13, FJXRY23
      double precision tbase1, xbase1, dbase1, coolunit, mh
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     FIRST, DEFINE SOME CONSTANTS
c
      PI        = pi_val
      EV2HZ     = 2.415d14       ! from eV to Hz

      mh        = mass_h         ! DPC

c
c     Clear sums used for X-ray compton temperature
c
      TXRY      = 0.d0
      EXRY      = 0.d0
c
c     Clear sums used in approximate radiative shield
c
      FJXRY11   = 0.d0
      FJXRY12   = 0.d0
      FJXRY13   = 0.d0
      FJXRY21   = 0.d0
      FJXRY22   = 0.d0
      FJXRY23   = 0.d0
c
c     Clear photo-heating sums
c
      PHTH      = 0.d0
      PHTHE     = 0.d0
      PHTHE2    = 0.d0
      PHTLAMH   = 0.d0
      PHTLAMHE  = 0.d0
      PHTLAMHE2 = 0.d0
c
c     Integrate over frequencies to get rates which depend on the
c         shape of the spectrum.
c
      DO N11=2,NFREQ
c
c        BIN WIDTH IN Hz
c
         FNUDEL    = (10.d0**((N11-0.5d0)*FREQDEL)
     .               -10.d0**((N11-1.5d0)*FREQDEL)
     .               )*EV2HZ
c
c        AVGERAGE ENERGY OF THE BIN IN QUESTION IN ergs
c
         AVGNU     = 10.d0**((N11-1.d0)*FREQDEL)*ev2erg
c
c        energy difference (in ergs) between mean frequency and cutoffs
c
         FNUMNUH0  = AVGNU - 13.6d0*ev2erg
         FNUMNUHE0 = AVGNU - 24.6d0*ev2erg
         FNUMNUHE20= AVGNU - 54.4d0*ev2erg
c
c        This is bin energy in units of the electron rest mass
c
         QNICK  = (10.d0**((N11-1.d0)*FREQDEL))/5.12d5
c
c        This is the Klein-Nishina cross-section
c
         IF(QNICK.LE.0.03d0) THEN
           CROSSS = 1.d0 - 2.d0*QNICK
         ELSE
           CROSSS = 3.d0/8.d0/QNICK*((1.d0-2.d0/QNICK-2.d0/QNICK**2)
     .             *LOG(1.d0+2.d0*QNICK) + 1.d0/2.d0
     .             -1.d0/2.d0/(1.d0+2.d0*QNICK)**2 + 4.d0/QNICK)
         ENDIF
C
C        1) FOR COMPTON HEATING
C
C        X-RAY BACKGROUND ENERGY DENSITY MULTIPLIED BY \SIGMA_{e,photon},
C        IN 1.E-21 erg/cm^3*cm^2, NOTE THAT INUTOT() IS IN
C        1.E-21 erg/cm^2/sec/hz/sr
C
         EXRY   = EXRY  + 4.d0*PI/c_light*INUTOT(N11)
     .                       *FNUDEL*CROSSS
C
C        X-RAY BACKGROUND EFFECTIVE TEMPERATURE
C
         TXRY   = TXRY  + 4.d0*PI/c_light*INUTOT(N11)*AVGNU
     .                       *FNUDEL*CROSSS
C
C        2) FOR PHOTOIONIZATION RATES
C
C        PHOTON-IONIZATION RATE DEFINED AS
C        PHTH   = \int_nu_0(H)^\infty 4*PI*sigma_H(nu)*I_nu/(h nu) d nu
C        IN 1.E-21/sec
C
C        HYDROGEN PHOTON IONIZATION RATE
C
         PHTH   = PHTH   + SIGH(N11)  *4.d0*PI*INUTOT(N11)/AVGNU*FNUDEL
C
C        HELIUM I PHOTON IONIZATION RATE
C
         PHTHE  = PHTHE  + SIGHE(N11) *4.d0*PI*INUTOT(N11)/AVGNU*FNUDEL
C
C        HELIUM II PHOTON IONIZATION RATE
C
         PHTHE2 = PHTHE2 + SIGHE2(N11)*4.d0*PI*INUTOT(N11)/AVGNU*FNUDEL
c
c        If using approximate radiation-shielding, compute the
c          mean optical depth
c
         if (iradshield .eq. 1) then
c     
            IF(FNUMNUH0.GE.0.0) THEN
               FJXRY11  = FJXRY11 + 4.d0*PI*INUTOT(N11)/AVGNU*FNUDEL
               FJXRY21  = FJXRY21 + 4.d0*PI*INUTOT(N11)/AVGNU*FNUMNUH0
     .                            *FNUDEL
            ENDIF
C
            IF(FNUMNUHE0.GE.0.0) THEN
               FJXRY12  = FJXRY12 + 4.d0*PI*INUTOT(N11)/AVGNU*FNUDEL
               FJXRY22  = FJXRY22 + 4.d0*PI*INUTOT(N11)/AVGNU*FNUMNUHE0
     .                            *FNUDEL
            ENDIF
C
            IF(FNUMNUHE20.GE.0.0) THEN
               FJXRY13  = FJXRY13 + 4.d0*PI*INUTOT(N11)/AVGNU*FNUDEL
               FJXRY23  = FJXRY23 + 4.d0*PI*INUTOT(N11)/AVGNU*FNUMNUHE20
     .                            *FNUDEL
            ENDIF
c
         endif
C
C        3) FOR PHOTOIONIZATION HEATING RATES
C           PHOTON IONIZATION HEATING RATE IN 1.E-21 erg/sec
C
         PHTLAMH   = PHTLAMH   + SIGH(N11)
     .                   *4.d0*PI*INUTOT(N11)/AVGNU*FNUMNUH0*FNUDEL
         PHTLAMHE  = PHTLAMHE  + SIGHE(N11)
     .                   *4.d0*PI*INUTOT(N11)/AVGNU*FNUMNUHE0*FNUDEL
         PHTLAMHE2 = PHTLAMHE2 + SIGHE2(N11)
     .                   *4.d0*PI*INUTOT(N11)/AVGNU*FNUMNUHE20*FNUDEL
C
      ENDDO
C
C
C     X-RAY BACKGROUND RADIATION FIELD TEMPERATURE CONVERTED TO IN KELVIN
C
      TXRY      = MIN(1.d-7,TXRY/MAX(1.d-10,EXRY))/(4.d0*kboltz)
c
c     For the approximate radiative shield, compute the
c        AVERAGE CROSS SECTION
c
      if (iradshield .eq. 1) then
c
         AVGSIGH   = PHTH/(1.d-20+FJXRY11)
         AVGSIGHE  = PHTHE/(1.d-20+FJXRY12)
         AVGSIGHE2 = PHTHE2/(1.d-20+FJXRY13)
         AVGSIGHP  = PHTLAMH/(1.d-20+FJXRY21)
         AVGSIGHEP = PHTLAMHE/(1.d-20+FJXRY22)
         AVGSIGHE2P= PHTLAMHE2/(1.d-20+FJXRY23)
c         WRITE(40,*)'NCYC,ZR,AVGSIGH,AVGSIGHE,AVGSIGHE2=',
c     .               NCYC,ZR,AVGSIGH,AVGSIGHE,AVGSIGHE2,
c     .               AVGSIGHP,AVGSIGHEP,AVGSIGHE2P
c
      endif
C
C     X-RAY BACKGROUND RADIATION ENERGY DENSITY IN erg/cm^3
C
      EXRY      = EXRY*1.d-21
      EXRY      = MAX(0.d0,EXRY)
C
C     PHOTO-IONIZATION RATES IN 1/sec
C       (now in AMR code units -- see calc_rates from more details)
C
      PHTH      = PHTH     *1.D-21 * utim
      PHTH      = MAX(0.d0,PHTH)
      PHTHE     = PHTHE    *1.D-21 * utim
      PHTHE     = MAX(0.d0,PHTHE)
      PHTHE2    = PHTHE2   *1.D-21 * utim
      PHTHE2    = MAX(0.d0,PHTHE2)
C
C     PHOTO-IONIZATION HEATING RATES IN 1.D-30 erg/sec
C       (now in AMR code units -- see calc_rates from more details)
C
      tbase1 = utim
      xbase1 = uxyz/(aye*uaye)    ! uxyz is [x]*a     
      dbase1 = urho*(aye*uaye)**3 ! urho is [dens]/a^3
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      PHTLAMH   = PHTLAMH  *1.D-21 / coolunit
      PHTLAMH   = MAX(0.d0,PHTLAMH)
      PHTLAMHE  = PHTLAMHE *1.D-21 / coolunit
      PHTLAMHE  = MAX(0.d0,PHTLAMHE)
      PHTLAMHE2 = PHTLAMHE2*1.D-21 / coolunit
      PHTLAMHE2 = MAX(0.d0,PHTLAMHE2)
C
      RETURN
      END
