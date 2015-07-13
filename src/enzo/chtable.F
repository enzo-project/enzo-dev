#include "fortran.def"
#include "phys_const.def"
#include "error.def"

c=======================================================================
c//////////////////////  SUBROUTINE CHTABLE  \\\\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine chtable(NFREQ, FREQDEL, SIGH, SIGHE, SIGHE2,
     &                   SPSTAR1, SPSTAR2)

c
c  GENERATES THE LOOK-UP TABLE USED FOR MANY ATOMIC RATES
c       (CURRENTLY ONLY RADIATIVE)
c
c  written by: Renyue Cen
c  date:       
c  modified1:  September, 1999 by Greg Bryan; converted to AMR
c  modified2:  February, 2004 by Robert Harkness
c              Remove obsolete syntax
c
c  PURPOSE:
C     THIS ROUTINE CREATES A LOOK UP TABLE FOR COLLISIONAL IONIZATION
C     AND RECOMBINATION COEFFICIENTS, AND COOLING AND HEATING TERMS,
C     EXCLUDING PHOTOIONIZATION AND PHOTOHEATING, COMPTON HEATING/COOLING
C   
C     THE LOOKUP TABLE IS A FUNCTION OF TEMPERATURE ONLY
c
c  INPUTS:
c    NFREQ    - Number of frequency bins
c    FREQDEL  - space between frequency bins, in log10(eV)
c
c  OUTPUTS:
c    SIGH     - HI photo-ionization heating cross-section
c    SIGHE    - HeI photo-ionization heating cross-section
c    SIGHE2   - HeII photo-ionization heating cross-section
c    SPSTAR1  - normalized shape of stellar radiation field
c    SPSTAR2  - normalized shape of quasar radiation field
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
      INTG_PREC NFREQ
      R_PREC    FREQDEL, SIGH(NFREQ), SIGHE(NFREQ), SIGHE2(NFREQ),
     &        SPSTAR1(NFREQ), SPSTAR2(NFREQ)
c
c  Parameters
c
      R_PREC    EV2HZ, PI
c
c  Locals
c
      INTG_PREC N, NTI
      R_PREC    E0, sigma0, ya, P, yw, y0, y1, x, y, Fy, FNU, FNU2,
     &        SPTMO, SUM, TBB, RAT, FLAMRAT, FNURAT, ALA, A,
     &        FNUPOWBK1, FNUPOWBK2, SPITP(9), FNUITP(9), WORK(9),
     &        FNU1
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c     Set some cosntants
c
      PI        = pi_val
      EV2HZ     = 2.415e14_RKIND
C
c-----------------------------------------------------------------------
C     CREATE THE CROSS SECTIONS
C
      E0 = 13.61_RKIND 
      sigma0 = 9.492e-16_RKIND
      ya = 1.469_RKIND
      P  = 3.188_RKIND
      yw = 2.039_RKIND
      y0 = 0.4434_RKIND
      y1 = 2.136_RKIND
      DO N=1,NFREQ
         FNU  = 10._RKIND**((N-1._RKIND)*FREQDEL)
C  
C        i) HYDROGEN I PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU .GT. 13.6_RKIND) THEN
           FNU1    = SQRT(FNU/13.6_RKIND -1._RKIND)
           SIGH(N) = 1.18d-11*FNU**(-4)*EXP(-4._RKIND*ATAN(FNU1)/FNU1)
     .                       /(1._RKIND -EXP(-2._RKIND*PI/FNU1))
         ELSE 
           IF(FNU .EQ. 13.6_RKIND) THEN
             SIGH(N) = 2.16e-13_RKIND*FNU**(-4)
           ELSE
             SIGH(N) = 0._RKIND
           ENDIF
         ENDIF
C
C
C        ii) HELIUM I PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU.GE.24.59_RKIND) THEN
c          SIGHE(N) = 1.13E-14*(FNU**(-2.05) -9.775*FNU**(-3.05))
           x       = FNU/E0-y0
           y       = sqrt(x*x+y1*y1)
           Fy      = ((x-1._RKIND)**2+yw**2)*y**(0.5_RKIND*P-5.5_RKIND)
     .             *(1._RKIND+sqrt(y/ya))**(-P)
           SIGHE(N) = sigma0*Fy
         ELSE
           SIGHE(N) = 0._RKIND
         ENDIF
C
C
C        iii) HELIUM II PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU.GT.54.4_RKIND) THEN
           FNU2      = SQRT(FNU/54.4_RKIND -1._RKIND)
           SIGHE2(N) = 7.55e-10_RKIND*FNU**(-4)
     .          * EXP(-4._RKIND*ATAN(FNU2)/FNU2)
     .          / (1._RKIND -EXP(-2._RKIND*PI/FNU2))
         ELSE
           IF(FNU.EQ.54.4_RKIND) THEN
             SIGHE2(N) = 1.38e-11_RKIND*FNU**(-4)
           ELSE
             SIGHE2(N) = 0._RKIND
           ENDIF
         ENDIF
      ENDDO
C
c-----------------------------------------------------------------------
C     GENERATE PHOTOIONIZATION SPECTRUM FOR MASSIVE STARS
C
      FNUITP(1)= 10._RKIND
      SPITP(1) = 0.291_RKIND

      FNUITP(2)= 12.6_RKIND
      SPITP(2) = 0.287_RKIND

      FNUITP(3)= 15.85_RKIND
      SPITP(3) = 0.283_RKIND

      FNUITP(4)= 19.95_RKIND
      SPITP(4) = 0.28_RKIND

      FNUITP(5)= 25.12_RKIND
      SPITP(5) = 0.21_RKIND

      FNUITP(6)= 30.63_RKIND
      SPITP(6) = 0.13_RKIND

      FNUITP(7)= 39.81_RKIND
      SPITP(7) = 0.081_RKIND

      FNUITP(8)= 50.12_RKIND
      SPITP(8) = 0.018_RKIND

      FNUITP(9)= 60._RKIND
      SPITP(9) = 0.01_RKIND
C
      DO N=1,NFREQ
         FNU = 10._RKIND**((N-1._RKIND)*FREQDEL)
         IF(FNU.LT.50.12_RKIND) THEN
           CALL IN_EXTP(FNUITP,SPITP,9_IKIND,FNU,SPTMO,WORK)
           SPSTAR1(N) = SPTMO
         ELSE
           SPSTAR1(N) = 10._RKIND**(-(FNU-60._RKIND) / 
     .           (60._RKIND-50.12_RKIND)
     .           *(LOG10(0.018_RKIND)-LOG10(0.01_RKIND))
     .           + LOG10(0.01_RKIND))
         ENDIF
      ENDDO
C
      SUM = 0._RKIND
      DO N=1,NFREQ
         FNU = 10._RKIND**((N-1._RKIND)*FREQDEL)
         IF(FNU.GE.13.6_RKIND) THEN
           SUM = SUM + SPSTAR1(N)*EV2HZ*FNU
     .           *(10._RKIND**(0.5_RKIND*FREQDEL)
     .           -10._RKIND**(-0.5_RKIND*FREQDEL))
         ENDIF
      ENDDO
C
      OPEN(15,FILE='STAR.sp')
      DO N=1,NFREQ
         FNU = 10._RKIND**((N-1._RKIND)*FREQDEL)
C
C        SPSTAR1(N) IS NORMALIZED SUCH THAT 
C        \int_(13.6 eV)^{\infty} SPSTAR1(N) dNU = 1 erg
C        WHERE dNU is in hz
         SPSTAR1(N) = SPSTAR1(N)/SUM
         WRITE(15,505)N,FNU,SPSTAR1(N),SPSTAR1(N)*FNU
     .                  ,SIGH(N),SIGHE(N),SIGHE2(N)
      ENDDO
      CLOSE(15)
C
  505 FORMAT(I4,6(1X,E12.4))
C
c-----------------------------------------------------------------------
C     GENERATE PHOTOIONIZATION SPECTRUM FOR QUASAR
C
C     BLACK BODY TEMPERATURE IN KELVIN
c     TBB     = 2.6e4_RKIND
      TBB     = 3.0e4_RKIND
C     BLACK BODY TEMPERATURE IN eV
      TBB     = Tbb/1.16e4_RKIND
      RAT     = 0.3_RKIND
      FLAMRAT = 5450._RKIND
c     h\nu IN eV
      FNURAT  = hplanck*c_light/(FLAMRAT*1.e-8_RKIND)/
     &     ev2erg
c
      ALA     = -1.3_RKIND
c
      A       = RAT/FNURAT**(3._RKIND-ALA)*(EXP(FNURAT/TBB)-1._RKIND) 
c
c
      FNUPOWBK1= 1.e3_RKIND
      FNUPOWBK2= 1.e5_RKIND
c
      DO N=1,NFREQ
         FNU = 10._RKIND**((N-1._RKIND)*FREQDEL)
         IF(FNU.LT.FNUPOWBK1) THEN
           SPSTAR2(N) = A*FNU**3/(exp(min(50._RKIND,FNU/Tbb))-1._RKIND) 
     .                 + FNU**ala
         ELSE
           IF(FNU.LT.FNUPOWBK2) THEN
             SPSTAR2(N) = FNU**ala*(1._RKIND+(FNU/FNUPOWBK1)**0.7_RKIND)
     .             / 2._RKIND
           ELSE
             SPSTAR2(N) = FNU**ala*(1._RKIND+(FNU/FNUPOWBK1)**0.7_RKIND)
     .             / 2._RKIND / (1._RKIND+(FNU/FNUPOWBK2)**1._RKIND)
     .             * 2._RKIND
             ENDIF
         ENDIF
      ENDDO
C
      SUM = 0._RKIND
      DO N=1,NFREQ
         FNU = 10._RKIND**((N-1._RKIND)*FREQDEL)
         IF(FNU.GE.13.6_RKIND) THEN
           SUM = SUM + SPSTAR2(N)*EV2HZ*FNU
     .           *(10._RKIND**(0.5_RKIND*FREQDEL)
     .           -10._RKIND**(-0.5_RKIND*FREQDEL))
         ENDIF
      ENDDO
C
      OPEN(15,FILE='QUASAR.sp')
      DO N=1,NFREQ
         FNU        = 10._RKIND**((N-1._RKIND)*FREQDEL)
         SPSTAR2(N) = SPSTAR2(N)/SUM
C
C        SPSTAR2(N) IS NORMALIZED SUCH THAT 
C        \int_(13.6 eV)^{\infty} SPSTAR2(N) dNU = 1 erg
C        WHERE dNU is in hz
         WRITE(15,505)N,FNU,SPSTAR2(N),SPSTAR2(N)*FNU
     .                 ,SIGH(N),SIGHE(N),SIGHE2(N)
      ENDDO
      CLOSE(15)
C
C
c-----------------------------------------------------------------------
C     NOW GENERATE TABLE AS A FUNCTION OF TEMPERATURE 
C     FOR IONIZATION COEFFICIENTS AND COOLING/HEATING 
C     FOR A PLASMA OF HYDROGEN AND HELIUM WITH PRIMORDIAL COMPOSITION 
C
#ifdef USE_CEN_RATES
C      
      DO NTI=1,NT
         TEMP    = 10._RKIND**(TMIN+(NTI-1._RKIND)*TDEL)
         TLIM    = 1._RKIND/(1._RKIND+SQRT(TEMP*1.e-5_RKIND))
         TSQRT   = SQRT(TEMP)
         TEMP4   = TEMP*1.e-4_RKIND
         TEMP4LOG = LOG10(TEMP4)
         TEMP4LOG4= LOG10(TEMP4/4._RKIND)
C
C        IONIZATION COEFFICIENTS IN cm^3/sec
C        RECOMBINATION COEFFICIENTS
C        HYDROGEN RECOMBINATION COEFFICIENT
         CHTERMS(1,NTI) = 1.567e-1_RKIND/(SQRT(TEMP/1.023_RKIND)
     .        *(1._RKIND+SQRT(TEMP/1.023_RKIND))**(
     .                                    1._RKIND-0.7545_RKIND)
     .        *(1._RKIND+SQRT(TEMP/6.8628e5_RKIND))**(
     .                                    1._RKIND+0.7545_RKIND))
C
C        HELIUMN II RECOMBINATION COEFFICIENT
         CHTERMS(2,NTI) = 3.294e-11_RKIND/(SQRT(TEMP/15.54_RKIND)
     .        *(1._RKIND+SQRT(TEMP/15.54_RKIND))**(1._RKIND-0.691_RKIND)
     .        *(1._RKIND+SQRT(TEMP/3.676d7))**(1._RKIND+0.691_RKIND))
C
C        HELIUMN III RECOMBINATION COEFFICIENT
         CHTERMS(3,NTI) = 2._RKIND*1.567e-10_RKIND / 
     .        (SQRT(TEMP/4._RKIND/1.023_RKIND)*(1._RKIND +
     .        SQRT(TEMP/4._RKIND/1.023_RKIND))**(1._RKIND-0.7545_RKIND)
     .        *(1._RKIND+SQRT(TEMP/4._RKIND/6.8628d5))**(1._RKIND0 +
     .        0.7545_RKIND))
C
C
C        COLLISIONAL IONIZATION COEFFICIENTS
C        HYDROGEN COLLISIONAL IONIZATION COEFFICIENT
         CHTERMS(4,NTI) = 2.91e-8_RKIND *
     .        (157809.1_RKIND/TEMP)**(0.39_RKIND)
     .        *EXP(-MIN(50.D0,157809.1_RKIND/TEMP))
     .        /(0.232_RKIND+157809.1_RKIND/TEMP)
C
C        HELIUM I COLLISIONAL IONIZATION COEFFICIENT
         CHTERMS(5,NTI) = 1.75e-8_RKIND * 
     .        (285335.4_RKIND/TEMP)**(0.35_RKIND)
     .        *EXP(-MIN(50.D0,285335.4_RKIND/TEMP))
     .        /(0.18_RKIND+285335.4_RKIND/TEMP)
C
C        HELIUM II COLLISIONAL IONIZATION COEFFICIENT
         CHTERMS(6,NTI) = 2.05e-9_RKIND *
     .        (1._RKIND+SQRT(631515._RKIND/TEMP))
     .        *(631515._RKIND/TEMP)**(0.25_RKIND)
     .        *EXP(-MIN(50._RKIND,631515._RKIND/TEMP))
     .        /(0.265_RKIND+631515._RKIND/TEMP)
C
C        HELIUMN II DIELECTRIC RECOMBINATION COEFFICIENT
         CHTERMS(7,NTI) = 1.9e-3_RKIND/TEMP**1.5_RKIND * 
     .        EXP(-MIN(50._RKIND,470000._RKIND/TEMP)) * 
     .        (1._RKIND+0.3_RKIND *
     .        EXP(-MIN(10._RKIND,94000._RKIND/TEMP)))
C
C
C
C
C        HEATING/COOLING RATES IN 1.E-30 erg cm^3/sec
C        HYDROGEN I COLLISIONAL IONIZATION COOLING TERM
         CHTERMS(8,NTI) = 6.33e11_RKIND
     .        *(157809.1_RKIND/TEMP)**(0.39_RKIND)
     .        *EXP(-MIN(50._RKIND,157809.1_RKIND/TEMP))
     .        /(0.232_RKIND+157809.1_RKIND/TEMP)
CCCC .                    *ELTR*HRNI 
C
C        HELIUM I COLLISIONAL IONIZATION COOLING TERM
         CHTERMS(9,NTI) = 6.89e11_RKIND
     .        *(285335.4_RKIND/TEMP)**(0.35_RKIND)
     .        *EXP(-MIN(50._RKIND,285335.4_RKIND/TEMP))
     .        /(0.18_RKIND+285335.4_RKIND/TEMP)
CCCC .                    *ELTR*HEMI
C
C        HELIUM I HE(2^3S) THREE-BODY COLLISIONAL IONIZATION COOLING TERM
         CHTERMS(10,NTI)= 5.01e3_RKIND/TEMP**0.1687_RKIND*TLIM
     .                    *EXP(-MIN(40._RKIND,5.53_RKIND/TEMP4))
CCCC .                    *ELTR**2*HEMII
C
C        HELIUM II COLLISIONAL IONIZATION COOLING TERM
         CHTERMS(11,NTI) = 1.78e11_RKIND*(1._RKIND
     .        +SQRT(631515._RKIND/TEMP))
     .        *(631515._RKIND/TEMP)**(0.25_RKIND)
     .        *EXP(-MIN(50._RKIND,631515._RKIND/TEMP))
     .        /(0.265_RKIND+631515._RKIND/TEMP)
CCCC .                    *ELTR*HEMII
C
C
C        HYDROGEN II RECOMBINATION COOLING TERM
         CHTERMS(12,NTI)= 10._RKIND**(5.6553_RKIND 
     /        + (0.2418_RKIND*TEMP4LOG
     .        -0.0745_RKIND*TEMP4LOG*TEMP4LOG
     .        -0.0316_RKIND*TEMP4LOG*TEMP4LOG*TEMP4LOG)
     .        /(1._RKIND+0.01685_RKIND*TEMP4LOG*TEMP4LOG*TEMP4LOG
     .        -0.00186_RKIND*TEMP4LOG*TEMP4LOG*TEMP4LOG*TEMP4LOG))
CCCC                *ELTR*HRNII
C
C        HELIUM II RECOMBINATION COOLING TERM
C        JORDI, I AM NOT CERTAIN THAT THIS EXPRESSION IS RIGHT??
C        THERE WAS A TEMP FACTOR MISSING
         CHTERMS(13,NTI) = 3.393e3_RKIND*TEMP/(SQRT(TEMP/15.54_RKIND)
     .        *(1._RKIND+SQRT(TEMP/15.54_RKIND))**(
     .                                      1._RKIND-0.691_RKIND)
     .        *(1._RKIND+SQRT(TEMP/3.676e7_RKIND))**(
     .                                      1._RKIND+0.691_RKIND))
CCCC .                    *ELTR*HEMII
C
C        HELIUM II DIELECTRONIC RECOMBINATION COOLING TERM
         CHTERMS(14,NTI)= 1.24e17_RKIND/TEMP**1.5_RKIND
     .        *EXP(-MIN(40._RKIND,47._RKIND/TEMP4))*(1._RKIND
     .        +0.3_RKIND*EXP(-MIN(40._RKIND,9.4_RKIND/TEMP4)))
CCCC .                    *ELTR*HEMII
C
C        HELIUM III RECOMBINATION COOLING TERM
         CHTERMS(15,NTI)= 8._RKIND*10._RKIND**(5.6553_RKIND 
     .        + (0.2418_RKIND*TEMP4LOG4
     .        -0.0745_RKIND*TEMP4LOG4*TEMP4LOG4
     .        -0.0316_RKIND*TEMP4LOG4*TEMP4LOG4*TEMP4LOG4)
     .        /(1._RKIND+0.01685_RKIND*TEMP4LOG4*TEMP4LOG4*TEMP4LOG4
     .        -0.00186_RKIND*TEMP4LOG4*TEMP4LOG4*TEMP4LOG4*TEMP4LOG4))
CCCC .                    *ELTR*HEMIII
C
C
C        HYDROGEN I COLLISIONAL EXCITATION COOLING TERM
         CHTERMS(16,NTI)= 7.5e11_RKIND*TLIM
     .        *EXP(-MIN(40._RKIND,11.83_RKIND/TEMP4))
CCCC .                    *ELTR*HRNI
C
C        HELIUM I COLLISIONAL EXCITATION COOLING TERM, HE(n=2,3,4 TRIPLETS) 
         CHTERMS(17,NTI)= 9.10e3_RKIND/TEMP**0.1687_RKIND*TLIM
     .        *EXP(-MIN(40._RKIND,1.32_RKIND/TEMP4))
CCCC .                    *ELTR**2*HEMII
C
C        HELIUM II COLLISIONAL EXCITATION COOLING TERM, He^+ (n=2)
         CHTERMS(18,NTI)= 5.54e13_RKIND/TEMP**0.397_RKIND*TLIM
     .                    *EXP(-MIN(40._RKIND,47.36_RKIND/TEMP4))
CCCC .                    *ELTR*HEMII
C
C        LITHIUM COLLISIONAL EXCITATION COOLING TERM, FIRST EXCITED (1.85eV)
C        STATE, IT DECAYS BY EMISSION OF A 6708A PHOTON
         CHTERMS(19,NTI)= 1.e1_RKIND*SQRT(10._RKIND*TEMP4)
     .                    *EXP(-MIN(40._RKIND,2.144_RKIND/TEMP4))
CCCC .                    *ELTR*(HRNI+HRNII)*RATLI2H
C
C        BREMSSTRAHLUNG RADIATION COOLING TERM, Gff=1.30 IS USED
         CHTERMS(20,NTI)= 1.3_RKIND*1.42e3_RKIND*TSQRT
CCCC .                    *ELTR*(HRNII+HEMII+4._RKIND*HEMIII)
C
      ENDDO 
C
#endif /* USE_CEN_RATES */
C
C
      RETURN
      END
C
C
C
C
C     NAME IN_EXTP
      SUBROUTINE IN_EXTP(XN,YN,N,X,Y,XX)
C
C     INPUT  : XN(N),YN(N),N,XX(N),X
C              WHERE XX(N) IS A WORKING ARRAY
C     OUTPUT : Y
C
C     PURPOSE: GIVEN N VALUES AT N POINTS, THIS ROUTINE INTERPOLATES OT
C              EXTRAPOLATES TO GIVE THE VALUE Y AT X
C
C     WRITTEN: BY RENYUE CEN, 11/7/92
C
C
      implicit none
#include "fortran_types.def"
C
      INTG_PREC N, NPOLY, NSHIFT, I
      R_PREC    XN, YN, X, Y, XX, XNTMP, YNTMP, DIRV1, DIRVN, DY
      DIMENSION XN(N),YN(N),XX(N)
      DIMENSION XNTMP(10),YNTMP(10)
C
C
      IF(X.GE.XN(2) .AND. X.LE.XN(N-1)) THEN
        DIRV1 = (YN(2)-YN(1))/(XN(2)-XN(1))
        DIRVN = (YN(N)-YN(N-1))/(XN(N)-XN(N-1))
        CALL SPLINE(XN,YN,N,DIRV1,DIRVN,XX)
        CALL SPLINT(XN,YN,XX,N,X,Y)
      ELSE
        IF(X.LT.XN(2)) THEN
          IF(N.GE.10) THEN
            NPOLY = 10
          ELSE
            NPOLY = N
          ENDIF
          DO I=1,NPOLY
             XNTMP(I) = XN(I)
             YNTMP(I) = YN(I)
          ENDDO
          CALL RATINT(XNTMP,YNTMP,NPOLY,X,Y,DY)
        ELSE
          IF(N.GE.10) THEN
            NPOLY  = 10
            NSHIFT = N - 10
            DO I=1,NPOLY
               XNTMP(I) = XN(NSHIFT+I)
               YNTMP(I) = YN(NSHIFT+I)
            ENDDO
          ELSE
            NPOLY  = N
            DO I=1,NPOLY
               XNTMP(I) = XN(I)
               YNTMP(I) = YN(I)
            ENDDO
          ENDIF
          CALL RATINT(XNTMP,YNTMP,NPOLY,X,Y,DY)
        ENDIF
      ENDIF
c
c
      RETURN
      END
c
c
c
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      implicit none
#include "fortran_types.def"
      INTG_PREC N, NMAX, I, K
      R_PREC    X, Y, YP1, YPN, Y2, U, SIG, P, QN, UN
      PARAMETER (NMAX=1000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99e3_RKIND) THEN
        Y2(1)=0._RKIND
        U(1)=0._RKIND
      ELSE
        Y2(1)=-0.5_RKIND
        U(1)=(3._RKIND/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2._RKIND
        Y2(I)=(SIG-1._RKIND)/P
        U(I)=(6._RKIND*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     &      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99e3_RKIND) THEN
        QN=0._RKIND
        UN=0._RKIND
      ELSE
        QN=0.5_RKIND
        UN=(3._RKIND/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1._RKIND)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
c
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      implicit none
#include "fortran_types.def"
      INTG_PREC N, KLO, KHI, K
      R_PREC    XA, YA, Y2A, X, Y, H, A, B
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)

      if (h.eq.0._RKIND) then
        write(0,'("Bad XA input in SPLINT")')
        ERROR_MESSAGE
      end if

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6._RKIND
      RETURN
      END
C
C
      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      implicit none
#include "fortran_types.def"
      INTG_PREC N, NMAX, NS, I, M
      R_PREC    XA, YA, X, Y, DY, TINY, C, D, HH, H, W, T, DD
      PARAMETER (NMAX=10,TINY=1.E-25)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0._RKIND)THEN
          Y=YA(I)
          DY=0._RKIND
          RETURN
        ELSE IF (H.LT.HH) THEN
          NS=I
          HH=H
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)+TINY
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          W=C(I+1)-D(I)
          H=XA(I+M)-X
          T=(XA(I)-X)*D(I)/H
          DD=T-C(I+1)

          if (dd.eq.0._RKIND) then
            write(0,'("DD = 0 in RATINT")')
            ERROR_MESSAGE
          end if

          DD=W/DD
          D(I)=C(I+1)*DD
          C(I)=T*DD
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END
