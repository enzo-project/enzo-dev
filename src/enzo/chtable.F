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
c
c  Arguments
c
      integer NFREQ
      real    FREQDEL, SIGH(NFREQ), SIGHE(NFREQ), SIGHE2(NFREQ),
     &        SPSTAR1(NFREQ), SPSTAR2(NFREQ)
c
c  Parameters
c
      real    EV2HZ, PI
c
c  Locals
c
      integer N, NTI
      real    E0, sigma0, ya, P, yw, y0, y1, x, y, Fy, FNU, FNU2,
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
      EV2HZ     = 2.415d14 
C
c-----------------------------------------------------------------------
C     CREATE THE CROSS SECTIONS
C
      E0 = 13.61d0 
      sigma0 = 9.492d-16
      ya = 1.469d0
      P  = 3.188d0
      yw = 2.039d0
      y0 = 0.4434d0
      y1 = 2.136d0
      DO N=1,NFREQ
         FNU  = 10.d0**((N-1.d0)*FREQDEL)
C  
C        i) HYDROGEN I PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU .GT. 13.6d0) THEN
           FNU1    = SQRT(FNU/13.6d0 -1.d0)
           SIGH(N) = 1.18d-11*FNU**(-4)*EXP(-4.d0*ATAN(FNU1)/FNU1)
     .                       /(1.d0 -EXP(-2.d0*PI/FNU1))
         ELSE 
           IF(FNU .EQ. 13.6d0) THEN
             SIGH(N) = 2.16d-13*FNU**(-4)
           ELSE
             SIGH(N) = 0.d0
           ENDIF
         ENDIF
C
C
C        ii) HELIUM I PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU.GE.24.59d0) THEN
c          SIGHE(N) = 1.13E-14*(FNU**(-2.05) -9.775*FNU**(-3.05))
           x       = FNU/E0-y0
           y       = sqrt(x*x+y1*y1)
           Fy      = ((x-1.d0)**2+yw**2)*y**(0.5d0*P-5.5d0)
     .             *(1.d0+sqrt(y/ya))**(-P)
           SIGHE(N) = sigma0*Fy
         ELSE
           SIGHE(N) = 0.d0
         ENDIF
C
C
C        iii) HELIUM II PHOTON-IONIZATION CROSS SECTION 
C
         IF(FNU.GT.54.4d0) THEN
           FNU2      = SQRT(FNU/54.4d0 -1.d0)
           SIGHE2(N) = 7.55d-10*FNU**(-4)*EXP(-4.d0*ATAN(FNU2)/FNU2)
     .                         /(1.d0 -EXP(-2.d0*PI/FNU2))
         ELSE
           IF(FNU.EQ.54.4d0) THEN
             SIGHE2(N) = 1.38d-11*FNU**(-4)
           ELSE
             SIGHE2(N) = 0.d0
           ENDIF
         ENDIF
      ENDDO
C
c-----------------------------------------------------------------------
C     GENERATE PHOTOIONIZATION SPECTRUM FOR MASSIVE STARS
C
      FNUITP(1)= 10.00d0
      SPITP(1) = 0.291d0

      FNUITP(2)= 12.60d0
      SPITP(2) = 0.287d0

      FNUITP(3)= 15.85d0
      SPITP(3) = 0.283d0

      FNUITP(4)= 19.95d0
      SPITP(4) = 0.280d0

      FNUITP(5)= 25.12d0
      SPITP(5) = 0.210d0

      FNUITP(6)= 30.63d0
      SPITP(6) = 0.130d0

      FNUITP(7)= 39.81d0
      SPITP(7) = 0.081d0

      FNUITP(8)= 50.12d0
      SPITP(8) = 0.018d0

      FNUITP(9)= 60.00d0
      SPITP(9) = 0.010d0
C
      DO N=1,NFREQ
         FNU = 10.d0**((N-1.d0)*FREQDEL)
         IF(FNU.LT.50.12d0) THEN
           CALL IN_EXTP(FNUITP,SPITP,9,FNU,SPTMO,WORK)
           SPSTAR1(N) = SPTMO
         ELSE
           SPSTAR1(N) = 10.d0**(-(FNU-60.d0)/(60.d0-50.12d0)
     .                        *(LOG10(0.018d0)-LOG10(0.010d0))
     .                        + LOG10(0.010d0))
         ENDIF
      ENDDO
C
      SUM = 0.d0
      DO N=1,NFREQ
         FNU = 10.d0**((N-1.d0)*FREQDEL)
         IF(FNU.GE.13.6d0) THEN
           SUM = SUM + SPSTAR1(N)*EV2HZ*FNU
     .                *(10.d0**(0.5d0*FREQDEL)-10.d0**(-0.5d0*FREQDEL))
         ENDIF
      ENDDO
C
      OPEN(15,FILE='STAR.sp')
      DO N=1,NFREQ
         FNU = 10.d0**((N-1.d0)*FREQDEL)
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
c     TBB     = 2.6e4
      TBB     = 3.0d4
C     BLACK BODY TEMPERATURE IN eV
      TBB     = Tbb/1.16d4
      RAT     = 0.3d0
      FLAMRAT = 5450.d0
c     h\nu IN eV
      FNURAT  = hplanck*c_light/(FLAMRAT*1.d-8)/
     &     ev2erg
c
      ALA     = -1.3d0
c
      A       = RAT/FNURAT**(3.d0-ALA)*(EXP(FNURAT/TBB)-1.d0) 
c
c
      FNUPOWBK1= 1.d3
      FNUPOWBK2= 1.d5
c
      DO N=1,NFREQ
         FNU = 10.d0**((N-1.d0)*FREQDEL)
         IF(FNU.LT.FNUPOWBK1) THEN
           SPSTAR2(N) = A*FNU**3/(exp(min(50.d0,FNU/Tbb))-1.d0) 
     .                 + FNU**ala
         ELSE
           IF(FNU.LT.FNUPOWBK2) THEN
             SPSTAR2(N) = FNU**ala*(1.d0+(FNU/FNUPOWBK1)**0.7d0)/2.d0
           ELSE
             SPSTAR2(N) = FNU**ala*(1.d0+(FNU/FNUPOWBK1)**0.7d0)/2.d0
     .                /(1.d0+(FNU/FNUPOWBK2)**1.d0)*2.d0
             ENDIF
         ENDIF
      ENDDO
C
      SUM = 0.d0
      DO N=1,NFREQ
         FNU = 10.d0**((N-1.d0)*FREQDEL)
         IF(FNU.GE.13.6d0) THEN
           SUM = SUM + SPSTAR2(N)*EV2HZ*FNU
     .                *(10.d0**(0.5d0*FREQDEL)-10.d0**(-0.5d0*FREQDEL))
         ENDIF
      ENDDO
C
      OPEN(15,FILE='QUASAR.sp')
      DO N=1,NFREQ
         FNU        = 10.d0**((N-1.d0)*FREQDEL)
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
         TEMP    = 10.d0**(TMIN+(NTI-1.d0)*TDEL)
         TLIM    = 1.d0/(1.d0+SQRT(TEMP*1.d-5))
         TSQRT   = SQRT(TEMP)
         TEMP4   = TEMP*1.d-4
         TEMP4LOG = LOG10(TEMP4)
         TEMP4LOG4= LOG10(TEMP4/4.d0)
C
C        IONIZATION COEFFICIENTS IN cm^3/sec
C        RECOMBINATION COEFFICIENTS
C        HYDROGEN RECOMBINATION COEFFICIENT
         CHTERMS(1,NTI) = 1.567d-10/(SQRT(TEMP/1.023d0)
     .                   *(1.d0+SQRT(TEMP/1.023d0))**(1.d0-0.7545d0)
     .                   *(1.d0+SQRT(TEMP/6.8628d5))**(1.d0+0.7545d0))
C
C        HELIUMN II RECOMBINATION COEFFICIENT
         CHTERMS(2,NTI) = 3.294d-11/(SQRT(TEMP/15.54d0)
     .                   *(1.d0+SQRT(TEMP/15.54d0))**(1.d0-0.691d0)
     .                   *(1.d0+SQRT(TEMP/3.676d7))**(1.d0+0.691d0))
C
C        HELIUMN III RECOMBINATION COEFFICIENT
         CHTERMS(3,NTI) = 2.d0*1.567d-10/(SQRT(TEMP/4.d0/1.023d0)
     .        *(1.d0+SQRT(TEMP/4.d0/1.023d0))**(1.d0-0.7545d0)
     .        *(1.d0+SQRT(TEMP/4.d0/6.8628d5))**(1.d00+0.7545d0))
C
C
C        COLLISIONAL IONIZATION COEFFICIENTS
C        HYDROGEN COLLISIONAL IONIZATION COEFFICIENT
         CHTERMS(4,NTI) = 2.91d-8*(157809.1d0/TEMP)**(0.39d0)
     .                    *EXP(-MIN(50.D0,157809.1d0/TEMP))
     .                    /(0.232d0+157809.1d0/TEMP)
C
C        HELIUM I COLLISIONAL IONIZATION COEFFICIENT
         CHTERMS(5,NTI) = 1.75d-8*(285335.4d0/TEMP)**(0.35d0)
     .                    *EXP(-MIN(50.D0,285335.4d0/TEMP))
     .                    /(0.18d0+285335.4d0/TEMP)
C
C        HELIUM II COLLISIONAL IONIZATION COEFFICIENT
         CHTERMS(6,NTI) = 2.05d-9*(1.d0+SQRT(631515.d0/TEMP))
     .                    *(631515.d0/TEMP)**(0.25d0)
     .                    *EXP(-MIN(50.D0,631515.d0/TEMP))
     .                    /(0.265d0+631515.d0/TEMP)
C
C        HELIUMN II DIELECTRIC RECOMBINATION COEFFICIENT
         CHTERMS(7,NTI) = 1.90d-3/TEMP**1.5d0 * 
     .        EXP(-MIN(50.D0,470000.d0/TEMP)) * 
     .        (1.d0+0.3d0*EXP(-MIN(10.D0,94000.d0/TEMP)))
C
C
C
C
C        HEATING/COOLING RATES IN 1.E-30 erg cm^3/sec
C        HYDROGEN I COLLISIONAL IONIZATION COOLING TERM
         CHTERMS(8,NTI) = 6.33d11*(157809.1d0/TEMP)**(0.39d0)
     .                    *EXP(-MIN(50.D0,157809.1d0/TEMP))
     .                    /(0.232d0+157809.1d0/TEMP)
CCCC .                    *ELTR*HRNI 
C
C        HELIUM I COLLISIONAL IONIZATION COOLING TERM
         CHTERMS(9,NTI) = 6.89d11*(285335.4d0/TEMP)**(0.35d0)
     .                    *EXP(-MIN(50.D0,285335.4d0/TEMP))
     .                    /(0.18d0+285335.4d0/TEMP)
CCCC .                    *ELTR*HEMI
C
C        HELIUM I HE(2^3S) THREE-BODY COLLISIONAL IONIZATION COOLING TERM
         CHTERMS(10,NTI)= 5.01d3/TEMP**0.1687d0*TLIM
     .                    *EXP(-MIN(40.d0,5.53d0/TEMP4))
CCCC .                    *ELTR**2*HEMII
C
C        HELIUM II COLLISIONAL IONIZATION COOLING TERM
         CHTERMS(11,NTI) = 1.78d11*(1.d0+SQRT(631515.d0/TEMP))
     .                    *(631515.d0/TEMP)**(0.25d0)
     .                    *EXP(-MIN(50.D0,631515.d0/TEMP))
     .                    /(0.265d0+631515.d0/TEMP)
CCCC .                    *ELTR*HEMII
C
C
C        HYDROGEN II RECOMBINATION COOLING TERM
         CHTERMS(12,NTI)= 10.d0**(5.6553d0 + (0.2418d0*TEMP4LOG
     .               -0.0745d0*TEMP4LOG*TEMP4LOG
     .               -0.0316d0*TEMP4LOG*TEMP4LOG*TEMP4LOG)
     .               /(1.d0+0.01685d0*TEMP4LOG*TEMP4LOG*TEMP4LOG
     .               -0.00186d0*TEMP4LOG*TEMP4LOG*TEMP4LOG*TEMP4LOG))
CCCC                *ELTR*HRNII
C
C        HELIUM II RECOMBINATION COOLING TERM
C        JORDI, I AM NOT CERTAIN THAT THIS EXPRESSION IS RIGHT??
C        THERE WAS A TEMP FACTOR MISSING
         CHTERMS(13,NTI) = 3.393d3*TEMP/(SQRT(TEMP/15.54d0)
     .                   *(1.d0+SQRT(TEMP/15.54d0))**(1.d0-0.691d0)
     .                   *(1.d0+SQRT(TEMP/3.676d7))**(1.d0+0.691d0))
CCCC .                    *ELTR*HEMII
C
C        HELIUM II DIELECTRONIC RECOMBINATION COOLING TERM
         CHTERMS(14,NTI)= 1.24d17/TEMP**1.5d0
     .                    *EXP(-MIN(40.d0,47.d0/TEMP4))*(1.d0
     .                      +0.3d0*EXP(-MIN(40.d0,9.4d0/TEMP4)))
CCCC .                    *ELTR*HEMII
C
C        HELIUM III RECOMBINATION COOLING TERM
         CHTERMS(15,NTI)= 8.d0*10.d0**(5.6553d0 + (0.2418d0*TEMP4LOG4
     .               -0.0745d0*TEMP4LOG4*TEMP4LOG4
     .               -0.0316d0*TEMP4LOG4*TEMP4LOG4*TEMP4LOG4)
     .               /(1.d0+0.01685d0*TEMP4LOG4*TEMP4LOG4*TEMP4LOG4
     .          -0.00186d0*TEMP4LOG4*TEMP4LOG4*TEMP4LOG4*TEMP4LOG4))
CCCC .                    *ELTR*HEMIII
C
C
C        HYDROGEN I COLLISIONAL EXCITATION COOLING TERM
         CHTERMS(16,NTI)= 7.5d11*TLIM*EXP(-MIN(40.d0,11.83d0/TEMP4))
CCCC .                    *ELTR*HRNI
C
C        HELIUM I COLLISIONAL EXCITATION COOLING TERM, HE(n=2,3,4 TRIPLETS) 
         CHTERMS(17,NTI)= 9.10d3/TEMP**0.1687d0*TLIM
     .                    *EXP(-MIN(40.d0,1.32d0/TEMP4))
CCCC .                    *ELTR**2*HEMII
C
C        HELIUM II COLLISIONAL EXCITATION COOLING TERM, He^+ (n=2)
         CHTERMS(18,NTI)= 5.54d13/TEMP**0.397d0*TLIM
     .                    *EXP(-MIN(40.d0,47.36d0/TEMP4))
CCCC .                    *ELTR*HEMII
C
C        LITHIUM COLLISIONAL EXCITATION COOLING TERM, FIRST EXCITED (1.85eV)
C        STATE, IT DECAYS BY EMISSION OF A 6708A PHOTON
         CHTERMS(19,NTI)= 1.00d10*SQRT(10.d0*TEMP4)
     .                    *EXP(-MIN(40.d0,2.144d0/TEMP4))
CCCC .                    *ELTR*(HRNI+HRNII)*RATLI2H
C
C        BREMSSTRAHLUNG RADIATION COOLING TERM, Gff=1.30 IS USED
         CHTERMS(20,NTI)= 1.30d0*1.42d3*TSQRT
CCCC .                    *ELTR*(HRNII+HEMII+4.0*HEMIII)
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
      PARAMETER (NMAX=1000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99d30) THEN
        Y2(1)=0.d0
        U(1)=0.d0
      ELSE
        Y2(1)=-0.5d0
        U(1)=(3.d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.d0
        Y2(I)=(SIG-1.)/P
        U(I)=(6.d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     &      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99d30) THEN
        QN=0.d0
        UN=0.d0
      ELSE
        QN=0.5d0
        UN=(3.d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.d0)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
c
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
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

      if (h.eq.0.d0) then
        write(0,'("Bad XA input in SPLINT")')
        ERROR_MESSAGE
      end if

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.d0
      RETURN
      END
C
C
      SUBROUTINE RATINT(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10,TINY=1.E-25)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      HH=ABS(X-XA(1))
      DO 11 I=1,N
        H=ABS(X-XA(I))
        IF (H.EQ.0.)THEN
          Y=YA(I)
          DY=0.d0
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

          if (dd.eq.0.) then
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
