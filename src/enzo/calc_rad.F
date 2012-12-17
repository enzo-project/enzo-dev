#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c/////////////////////  SUBROUTINE CALC_RAD  \\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine calc_rad(NFREQ, FREQDEL,
     &                    AAA, ANEW, OMEGA0, H, DT1, UTIM, 
     &                    NFREQT, FREQTDEL, TMIN,
     &                    JUSTBURN, EUVSTAR, EUVQUASAR,
     &                    SPSTAR1, SPSTAR2,
     &                    INU, INUG, INUS, INUQ,
     &                    JNU, JNUG, JNUS, JNUQ,
     &                    DENSQ, DSH2, DSHE2, DSHE3,
     &                    H1DEN, HE1DEN, HE2DEN,
     &                    SIGH, SIGHE, SIGHE2)
c
c  CALCULATES THE NEW, HOMOGENEOUS RADIATION FIELD
c
c  written by: Renyue Cen
c  date:       
c  modified1:  September, 1999 by Greg Bryan; converted to AMR
c
c  PURPOSE:
c    Given information about the H and He density fields and the
c      the star formation rate, calculates the new radiation field
c      more-or-less as described in Cen & Ostriker.
c
c  INPUTS:
c    NFREQ    - Number of frequency bins
c    FREQDEL  - space between frequency bins, in log10(eV)
c    NFREQT   - Number of temperature bins for emission calculation
c    FREQTDEL - space between temperature bins, in log10(K)
c    TMIN     - minimum temperature bin, in log10(K)
c    iradshield - INTG_PREC flag indicating if approximate radiative
c                 shielding should be used (0 - no, 1 - yes)
c    AAA      - original expansion factor (were AAA=1 means z=0)
c    ANEW     - new expansion factor
c    H1DEN    - mean HI density, in particles/cm^3
c    HE1DEN   - mean HeI density, in particles/cm^3
c    HE2DEN   - mean HeII density, in particles/cm^3
c    DT1      - time step since last time this was done (in code units)
c    UTIM     - time conversion factor from code units to seconds
c    DENSQ    - mean value of n_e * (sum(Z_i * n_i)) as a function of T
c               (units are in number^2/cm^6)
c    DENSH2   - mean value of n_e * n_HII as a function of T
c    DENSHE2  - mean value of n_e * n_HII as a function of T
c    DENSHE3  - mean value of n_e * n_HII as a function of T
c    OMEGA0   - Omega matter at z=0
c    H        - Hubble constant in 100 km/s/Mpc
c    INU      - radiation intensity in 10^-12 erg/cm^2/sec/hz/sr
c    SPSTAR1  - normalized shape of stellar radiation field
c    SPSTAR2  - normalized shape of quasar radiation field
c    SIGH     - HI photo-ionization heating cross-section
c    SIGHE    - HeI photo-ionization heating cross-section
c    SIGHE2   - HeII photo-ionization heating cross-section
c    JUSTBURN - integrated mean star formation rate (units of mean density)
c    EUVSTAR  - parameter (rest mass energy fraction converted to star UV)
c    EUVQUASAR - parameter (rest mass energy fraction converted to quasar UV)
c
c  OUTPUTS:
c    JNU      - total emissivity in units of 10^-21 erg/cm^2/sec^2/hz/sr
c    JNUG     - gas emissivity in units of 10^-21 erg/cm^2/sec^2/hz/sr
c    JNUS     - stellar emissivity in units of 10^-21 erg/cm^2/sec^2/hz/sr
c    JNUQ     - quasar emissivity in units of 10^-21 erg/cm^2/sec^2/hz/sr
c    INU      - total radiation intensity, updated
c    INUG     - gas radiation intensity, updated
c    INUS     - stellar radiation intensity, updated
c    INUQ     - quasar radiation intensity, updated
c
c  PARAMETERS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c  Arguments
c
      INTG_PREC NFREQ, NFREQT
      R_PREC    FREQDEL, FREQTDEL, TMIN, DT1, JNU(NFREQ), INU(NFREQ),
     &        JNUG(NFREQ), JNUS(NFREQ), JNUQ(NFREQ),
     &        INUG(NFREQ), INUS(NFREQ), INUQ(NFREQ),
     &        DENSQ(NFREQT), DSH2(NFREQT), DSHE2(NFREQT), DSHE3(NFREQT),
     &        AAA, ANEW, OMEGA0, H, UTIM, H1DEN, HE1DEN, HE2DEN,
     &        SPSTAR1(NFREQ), SPSTAR2(NFREQ), 
     &        SIGH(NFREQ), SIGHE(NFREQ), SIGHE2(NFREQ), JUSTBURN,
     &        EUVSTAR, EUVQUASAR
c
c  Parameters
c
      R_PREC    EV2HZ, PI
c
c  Locals
c
      INTG_PREC N11, N22, NSHFT, NUPP, NLOW
      R_PREC    FNU, TEMP, SUB2, SUB3, SUB4, UFF, UFB, DELTAOA, ETA,
     &        ZR, COEFUV1, COEFUV2, DTST, FKAPA
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     FIRST, DEFINE SOME CONSTANTS
c
      PI        = pi_val
      EV2HZ     = 2.415e14_RKIND       ! from eV to Hz
c
c     Set the timestep in seconds
c
      DTST      = DT1*UTIM
C
C     FREE-FREE EMISSION COEFFICIENT IN UNITS OF
C     1E-21 erg/cm^3/sec/hz/sr * cm/sec
C     SEE SPITZER "Physical Processes in the Interstellar Medium" PAGE 57
C
      UFF       = 1.3_RKIND*5.44e-18_RKIND*c_light
C
C     FREE-BOUND EMISSION COEFFICIENT IN UNITS OF
C     1E-21 erg/cm^3/sec/hz/sr * cm/sec
C
      UFB       = 3.66d-4*c_light
C     
C     COEFICIENT FOR THE MASS BURNING STARS PRODUCING UV 
C     RADIATION RATE IN UNITS 1.0E-21 erg/cm^2/sec/sr
C
      ZR = 1._RKIND/ANEW - 1._RKIND
      COEFUV1= 1.88e-8_RKIND*OMEGA0*H**2*(1._RKIND+ZR)**3
     .                *2.7e31_RKIND*EUVSTAR/4._RKIND/PI
C
      IF(ZR.LE.2._RKIND) THEN
         COEFUV2= 1.88e-8_RKIND*OMEGA0*H**2*(1._RKIND+ZR)**3
     .        *2.7e31_RKIND*EUVQUASAR/4._RKIND/PI
      ELSE
         COEFUV2= 1.88e-8_RKIND*OMEGA0*H**2*(1._RKIND+ZR)**3
     .        *2.7e31_RKIND*EUVQUASAR/4._RKIND/PI
     .        *EXP(-(ZR-2._RKIND)**2)
C     NOTE THAT /hz IS HANDLED IN coef.src
      ENDIF
c
c-----------------------------------------------------------------
c     Calculate emissivity
c
      DO N11=1,NFREQ
c
c        Set emissivity to zero and compute frequency of this bin, in eV
c
         JNU(N11) = 0._RKIND
         FNU      = 10._RKIND**(FREQDEL*(N11-1._RKIND)) 
c
c        Loop over bins in temperature space
c
         DO N22=2,NFREQT
c
c           Compute temperature of this bin
c
            TEMP     = 10._RKIND**((N22-1._RKIND)*FREQTDEL+TMIN)
c
c           HI recombination emissivity
c
            IF(FNU.GE.13.6_RKIND) THEN
              SUB2 = EXP(-MIN(50._RKIND,
     .              1.16e4_RKIND*(FNU-13.6_RKIND)/TEMP))/TEMP/SQRT(TEMP)
     .              * SQRT((FNU-13.6_RKIND)*ev2erg)
            ELSE 
              SUB2=0._RKIND
            ENDIF
c
c           HeI recombination emissivity
c
            IF(FNU.GE.24.6_RKIND) THEN
               SUB3 = EXP(-MIN(50._RKIND,
     .              1.16e4_RKIND*(FNU-24.6_RKIND)/TEMP))/TEMP/SQRT(TEMP)
     .              * SQRT((FNU-24.6_RKIND)*ev2erg)
            ELSE 
              SUB3=0._RKIND
            ENDIF
c
c           HeII recombination emissivity
c
            IF(FNU.GE.54.4_RKIND) THEN
               SUB4 = EXP(-MIN(50._RKIND,
     .              1.16e4_RKIND*(FNU-54.4_RKIND)/TEMP))/TEMP/SQRT(TEMP)
     .              * SQRT((FNU-54.4_RKIND)*ev2erg)
            ELSE 
              SUB4=0._RKIND
            ENDIF
c
c           Add emissivity to this frequency bin (from this temp bin)
c
            JNU(N11) = JNU(N11) 
c
c                          1) free-free
c
     .           + UFF*DENSQ(N22)/SQRT(TEMP)
     .           *EXP(-MIN(50._RKIND,FNU*1.16e4_RKIND/TEMP))
c
c                          2) free-bound (HI, HeI, HeII)
c
     .           + UFB*DSH2(N22)*SUB2
     .           *FNU*ev2erg
     .           *6.63e-11_RKIND/SQRT(TEMP) 
     .           /(TEMP/1.e3_RKIND)**0.2_RKIND
     .           /(1._RKIND+(TEMP/1.e6_RKIND)**0.7_RKIND)
C
     .           + UFB*DSHE2(N22)*SUB3
     .           *FNU*ev2erg
     .           *1.5e-10_RKIND/TEMP**0.6353_RKIND
C
     .           + UFB*DSHE3(N22)*SUB4
     .           *FNU*ev2erg
     .           *2.65e-10_RKIND/SQRT(TEMP)
     .           /(TEMP/1.e3_RKIND)**0.2_RKIND
     .           /(1._RKIND+(TEMP/1.e6_RKIND)**0.7_RKIND)
C
         ENDDO
         JNU(N11) = MAX(0._RKIND,JNU(N11))
      ENDDO
c
c-----------------------------------------------------------------
c     Compute new radiation field
c
      DELTAOA = AAA/ANEW
      ETA     = -LOG10(DELTAOA)/FREQDEL
      NSHFT   = ETA
c
      DO N11=1,NFREQ
c
c        Compute line kapa at this frequency
c           (H1DEN is mean HI density in particles/cm^3)
c
         FKAPA =(H1DEN *SIGH  (N11)
     .         + HE1DEN*SIGHE (N11)
     .         + HE2DEN*SIGHE2(N11))
c
c        Compute the two bins from which to read the radiation
c          (due to redshifting)
c
         NUPP = MIN(NFREQ,N11+1+NSHFT)
         NLOW = MIN(NFREQ,N11+NSHFT)
c
c        Compute new intensity by interpolating from old field
c          (with redshifting factor) and adding new emissivity
c
         INU(N11) =((INU(NUPP)*(ETA-NSHFT)
     .              +INU(NLOW)*(1._RKIND-ETA+NSHFT))*DELTAOA**3
c
c                  Emissivity from gas (free-free and free-bound)
c
     .             +(JNU(N11)
     .              )*DTST
c
c                  Emissivity from stars
c
     .             +SPSTAR1(N11)*JUSTBURN*COEFUV1
c
c                  Emissivity from quasars
c
     .             +SPSTAR2(N11)*JUSTBURN*COEFUV2

     .             )/(1._RKIND+c_light*FKAPA*DTST)
c
c        Save the various emissivitivies
c
         JNUG(N11) = JNU(N11)          
         JNUS(N11) = SPSTAR1(N11)*JUSTBURN*COEFUV1/DTST
         JNUQ(N11) = SPSTAR2(N11)*JUSTBURN*COEFUV2/DTST
c
c        Compute the intensities of the components (gas, stars, quasars)
c
         INUG(N11) = ((INUG(NUPP)*(ETA-NSHFT)
     .              +INUG(NLOW)*(1._RKIND-ETA+NSHFT))*DELTAOA**3
     .              +JNUG(N11)*DTST
     .               )/(1._RKIND+c_light*FKAPA*DTST)
         INUS(N11) =((INUS(NUPP)*(ETA-NSHFT)
     .              +INUS(NLOW)*(1._RKIND-ETA+NSHFT))*DELTAOA**3
     .              +SPSTAR1(N11)*JUSTBURN*COEFUV1
     .              )/(1._RKIND+c_light*FKAPA*DTST)
C
         INUQ(N11) =((INUQ(NUPP)*(ETA-NSHFT)
     .              +INUQ(NLOW)*(1._RKIND-ETA+NSHFT))*DELTAOA**3
     .              +SPSTAR2(N11)*JUSTBURN*COEFUV2
     .              )/(1._RKIND+c_light*FKAPA*DTST)
C
      ENDDO
c
      RETURN
      END
