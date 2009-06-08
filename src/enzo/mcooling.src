#include "fortran.def"
#include "error.def"

C     NAME MTLCLHT
      SUBROUTINE MTLCLHT(DEN,RADT,IPHOTT,CBOVCOOL,CBOVHEAT)
#ifdef CEN_METALS
C
C     THIS PROGRAM CALCULATES A COOLING/HEATING TABLE FOR METALS FROM
C     C TO NI EXCLUDING H AND HE, GIVEN AN IONIZING RADIATION FIELD
C     THE ASSUMPTIONS ARE: 1) ELECTRON NUMBER DENSITY IS DOMINATED BY
C                             THOSE FROM H AND HE, I.E., ELECTRON
C                             NUMBER DENSITY FROM METALS CAN BE IGNORED.
C                             THIS ASSUMPTION IS LIKELY TO BE VALID FOR
C                             COSMOLOGICAL SIMULATIONS WHERE METALLICITY
C                             IS NOT LARGE COMPARED TO THE SOLAR ONE.
C     
C                          2) METALS ARE IN IONIZATION EQUILIBRIUM.
C
C
C     *******************************************************************
C     11/20/96, R. Cen
C
C
C     INPUT:
C     NID    : NUMBER OF ENTRIES IN THE ELECTRON NUMBER DENSITY DIMENSION 
C     DENMIN : MINIMUM DENSITY     
C     DELD   : DENSITY INCREMENT IN LOG SCALE
C     NIT    : NUMBER OF ENTRIES IN THE TEMPERATURE DIMENSION 
C     TEMMIN : MINIMUM TEMPERATURE
C     DELT   : TEMPERATURE INCREMENT IN LOG SCALE
C
C     NFREQ  : NUMBER OF ENTRIES IN THE SPECTRUM 
C     FREQDEL: FREQUENCY INCREMENT IN LOG SCALE, MINIMUM=1eV
C     JNU    : THE INTENSITY OF THE SPECTRUM
C     IPHOT  : = 0 TO IGNORE  PHOTOIONIZATION 
C              > 0 TO INCLUDE PHOTOIONIZATION 
C     NOTE THAT IPHOT HAS TO BE SET BEFORE CALLING MTLCLHT
C
C
C
C     OUTPUT:
C     MTLCL  : COOLING RATE DUE TO METALS EXCLUDING H AND HE
C              KI* CBOVCOOL * N_H       IN 10^{-23} ERG/CM^3/S
C     MTLHT  : HEATING RATE DUE TO METALS EXCLUDING H AND HE
C              KI* CBOVHEAT * N_H       IN 10^{-23} ERG/CM^3/S
C
C     *******************************************************************
C
      INCLUDE "MTLPARAM.h"
C
      REAL DEN(NID)
      REAL RADT(NIB)
      REAL CBOVCOOL(NIT,NID)
      REAL CBOVHEAT(NIT,NID)
C
C
C
C
C
      IPHOT = IPHOTT
c     WRITE(*,*)'IPHOT = ',IPHOT
      IF(IPHOT.GT.0) THEN
        DO I=1,NIB
           ABIN(I) = RADT(I)
        ENDDO
      ENDIF
C
      DO ID=1,NID
         AN(ID) = DEN(ID)
      ENDDO
C
C
C     INITIALIZATION
      CALL MTLINIT
C
C     GENERATE COOLING/HEATING TABLE
      CALL METAL(CBOVCOOL,CBOVHEAT)
C 
c     WRITE(93,*)CBOVHEAT(1,1),CBOVHEAT(10,10),CBOVHEAT(100,100)
c     CALL FLUSH(93)
C
      RETURN
      END
C
C
C     NAME = MTLINIT
      SUBROUTINE MTLINIT
C
C     THIS ROUTINE CALCULATES THE PHOTOIONIZATION CROSS SECTION AND HEATING
C     TERMS THAT ARE INDEPENDENT OF TEMPERATURE AND IONIZATION FRACTION.
C     IT ALSO CALCULATES THE EMISSION SPECTRUM OF H AND HE
C     THIS SHOULD BE CALLED ONLY ONCE IN THE BEGINNING OF THE CALCULATION.
C
C     13 ELEMENTS.  H AND HE ARE ALWAYS CALCULATED SEPARATELY; 
C     THE REMAINING ELEMENTS ARE: C,N,O,NE,MG,SI,S,AR,CA,FE,NI.
C     ALL SUBROUTINES OTHER THAN THE MAIN ROUTINE AND ROUTINES IN ISODEN.F
C     CAME COURTESY OF JOHN C. RAYMOND IN 1988.  
C     SEVERAL MODIFICATIONS TO THEM HAVE BEEN MADE BY BOB BENJAMIN.  
C
C
C     SOME USEFUL IDENTIFICATIONS FOR FIGURING OUT THE PROGRAM:
C     RHY      : H+/H0
C     CNC(I,J) : CONCENTRATION OF ION J OF ELEMENT I
C     NJ(I)    : ATOMIC NUMBER OF ELEMENT I
C     C1J(I,J) : IONIZATION COEFFICIENT
C     C2J(I,J) : RECOMBINATION COEFFICIENT
C     CAJ(I,J) : AUGER IONIZATION COEFFICIENT
C
C     WJ(I,J)  : WAVE(I,J) : WAVELENGTH OF ATOMIC TRANSITION
C     E3J(I,J) : LINE EXCITATION ENERGY
C     FJ(I,J)  : F-VALUES, OSCILLATOR STRENGTH        
C
C     E(I)     : IONIZATION POTENTIAL
C     EA(I)    : IONIZATION POTENTIAL FOR INNERSHELL IONIZATION 
C     S2J(),S3J(),S4J(),S5J() : COEFFICIENTS FOR PHOTOIONIZ. CROSS SECTION
C
C
      INCLUDE "MTLPARAM.h"
C
      DIMENSION WEIGHT(12), SOLAR(12)
C
C     ATOMIC WEIGHT IN ATOMIC MASS UNITS FOR THE ELEMENTS IN THE PROGRAM
      DATA WEIGHT/4.00,12.01,14.01,16.00,20.18,24.31,28.09,32.06,
     .  39.95,40.08,55.85,58.71/
C
C     NEW SOLAR ABUNDACES FROM CLOUDY (7/28/94) 
      DATA SOLAR/10.991, 8.56, 8.049, 8.93, 8.09, 7.58, 7.55,
     .  7.21, 6.56, 6.36, 7.67, 6.246/
C
C
C
      DO I=1,NUME
C        SOLAR ABUNDANCES
         ABUNJ(I)=SOLAR(I)
      ENDDO
C
C
C
C     READ IN ATOMIC DATA 
      CALL ATREAD
C
C
      RETURN
      END
C
C
C
C     NAME METAL
      SUBROUTINE METAL(CBOVCOOL,CBOVHEAT)
C 
C
      INCLUDE "MTLPARAM.h"
C
      REAL CBOVCOOL(NIT,NID)
      REAL CBOVHEAT(NIT,NID)
C
      REAL ABNN(12)
      REAL PHOTR(12,30),AP1(12,30),AP2(12,30)
      REAL CNC(12,30,NID)
      REAL HEATR(12,30)
      REAL PCOOL(NID)
      REAL HEATING(NID)
      REAL DELNE(NID)
C
      DO NO=1,NUME
C        ABUNDANCE OF EACH ELEMENT FROM HE TO FE IN NUMBER 
C        RELATIVE TO HYDROGEEN
         ABNN(NO)=10.0**(ABUNJ(NO)-12.0)
      ENDDO
C
C
C
C PREGION 0
C$DOACROSS LOCAL(I,J,ID)
C$&       ,SHARE(CNC)
C     INITIALIZATION OF CNC
      DO I=1,12
         DO J=1,30
            DO ID=1,NID
               CNC(I,J,ID)=0.0
            ENDDO
         ENDDO
      ENDDO
C
C
C     GENERATE THE PHOTOIONIZATION CROSS SECTION FOR EACH IONS
C     CALCULATE THE PHOTOIONZ. RATE AND HEATING FOR A GIVEN FIELD ABIN(IB) 
      IF(IPHOT.GT.0) THEN
C PREGION 1
C$DOACROSS LOCAL(NO,NN)
C$&       ,SHARE(NJ,PHOTR,AP1,AP2,HEATR)
        DO NO=2,NUME
           NN=NJ(NO)
           CALL PHTR(NO,NN,PHOTR,HEATR)
           CALL APHTR(NO,NN,AP1,AP2,HEATR)
        ENDDO
      ELSE
        DO NO=2,NUME
           NN = NJ(NO)
           DO J = 1,NN
              HEATR(NO,J) = 0.0
           ENDDO
        ENDDO
      ENDIF
C
C
C
C     GENERATE ELECTRON DENSITY, COOLING TABLE
C     LOG(TEMP) = 3.0 + DELT*(IT-1),  IT=1,..NIT FOR NE AND COOLING TABLE
C     REAL METAL COOL = KI* CBOVCOOL * N_H       IN 10^{-23} ERG/CM^3/S
C     REAL METAL HEAT = KI* CBOVHEAT * N_H       IN 10^{-23} ERG/CM^3/S
C
C
C PREGION 2
C$DOACROSS LOCAL(IT,BREM,TEMPK,CNC,PCOOL,ID,HEATING,I,J
C$&             ,NO,NN,DELNE,FRACT,BRMZ,FKI,DENH)
C$&       ,SHARE(ABNN,NJ,AN,PHOTR,AP1,AP2,CBOVCOOL,CBOVHEAT)
C
      DO IT=1,NIT
         TEMPK  = 10.0**(TEMMIN+(IT-0.5)*DELT)
         BREM   = SQRT(TEMPK)*1.5*1.42E-4
C
         CALL FELINE(TEMPK,PHOTR,AP1,AP2,CNC,PCOOL)
C
         DO ID=1,NID
            HEATING(ID) = 0.0
            DO NO = 2,NUME
               NN = NJ(NO)
               DO J = 1,NN
                  HEATING(ID) = HEATING(ID) 
     .                        + HEATR(NO,J)*CNC(NO,J,ID)*ABNN(NO)
               ENDDO
            ENDDO
         ENDDO
C
C
C        COMPUTES THE ELECTRON NUMBER DNESITY FROM METALS AS A CHECK
C        TO ASSUMPTION 1
         DO ID=1,NID
            DELNE(ID) = 0.0 
            DO NO=2,NUME
               NN=NJ(NO)
               DO J=1,NN+1
                  FRACT=ABNN(NO)*CNC(NO,J,ID)
                  DELNE(ID)=DELNE(ID) + (J-1)*FRACT
               ENDDO
            ENDDO
C
C
            BRMZ=0.0
            DO NO=2,NUME
               DO J=1,NJ(NO)+1
                  BRMZ = BRMZ+(J-1)**2*ABNN(NO)*CNC(NO,J,ID)
               ENDDO
            ENDDO
C           NOTE THAT FKI: METALLICITY, DENH: HYDROGEN ATOM DENSITY NEEDS TO
C           BE CHANGED
            CBOVCOOL(IT,ID) = (PCOOL(ID) + BREM*BRMZ)*AN(ID)
            CBOVHEAT(IT,ID) = HEATING(ID)
         ENDDO
      ENDDO
C
C
      RETURN 
  776 FORMAT(I5,6E13.5)
      END
C     
C
C     NAME=APHTR
      SUBROUTINE APHTR(NO,N,AP1,AP2,HEATR)
C
      INCLUDE "MTLPARAM.h"
C
      REAL AP1(12,30),AP2(12,30),HEATR(12,30)
C
C     AUGER PHOTOIONIZATION : N=1 AND N=2 SHELLS
C     USES COX & DALTABUIT AND REILMAN & MANSON CROSS SECTIONS
      DIMENSION A2P(20,20),B2P(20,20)
      DATA A2P/20*0,12.60,20.58, 38*0., 4.520,3.945,5.546,7.213,
     .  36*0., 1.286,1.316,1.843,2.556,3.363,3.654, 34*0.,
     .  .6109,.6256,.8755,1.111,1.254,1.578,1.919,2.307, 32*0.,
     .  .4075,.4438,.4801,.5164,.6798,.8431,.9660,1.089,1.288,1.488,
     .  110*0., .1339,.1531,.1722,.1909,.2095,.2320,.2545,.3022,.3499,
     .  .3628,.3756,.4630,.5504,.5820,.6136,.6679, 24*0.,
     .  .1145,.1294,.1409,.1525,.1644,.1763,.1797,.1831,.2005,.2179,
     .  .2425,.2670,.3393,.4116,.4266,.4415,.5339,.6263, 42*0./
      DATA B2P/20*0.,.1887,-7.098, 38*0., 1.349,5.109,4.291,6.750,
     .  36*0., 2.408,2.546,3.116,3.060,2.981,5.273, 34*0.,
     .  1.377,1.456,1.782,1.930,2.384,2.572,3.023,3.451, 32*0.,
     .  .7788,.9714,1.164,1.356,1.488,1.620,1.869,2.118,2.261,2.405,
     .  110*0.,.3454,.3781,.4107,.4565,.5023,.5590,.6157,.6283,.6409,
     .  .7392,.8375,.8258,.8142,.9381,1.062,1.207, 24*0.,
     .  .2541,.2787,.3087,.3387,.3770,.4152,.4907,.5661,.6055,.6499,
     .  .6729,.7010,.6410,.5809,.6970,.8130,.7768,.7406, 42*0./
C
C
      HE = 0.0
      RN = N
      NN = N-2
      IF(N .LE. 2) RETURN
C
C
C     FACT1= 1.E-21/1.6E-12*1.E-18
C     FIRST FACTOR: 1.E-21 DUE TO THE ADOPTED UNITS OF THE PHOTOIONIZATION
C                   FIELD: 1.E-21 ERG/CM^2/SEC
C     SECOND FACTOR: CONVERTING eV to ERG
C     THIRD FACTOR: UNITS OF THE PHOTOIONIZATION CROSS SECTION
C
      FACT1= 1.E-9/1.6*1.E-18
C
C
C
C     FACT = 1.6E-12*1.E23 
C     FIRST FACTOR: CONVERTING BACK ERG to eV      
C     SECOND FACTOR: THE UNITS OF THE PHOTOHEATING IS OUTPUTED IN 1.E-23 ERG/SEC
      FACT = 1.6E11
C
C
C
      DO J=1,NN
         ISO = N-J+1
         RISO = ISO
         PA1 = 0.
         PA2 = 0.
         HE = 0.
         EKS = 13.6 * RN*RN * RISO ** (.24-.43/LOG10(RN))
C
         IK = max(1.0,log10(EKS/freqmin)/freqdel)
         IMAX = nib
         IF(ISO .GE. 11) IMAX = MIN(IK,nib)
         IMIN = max(1.0,log10(EAJ(NO,J)/freqmin)/freqdel+1)
         ebc  = freqmin*10**(freqdel*(imin-0.5)) 
         IF(EAJ(NO,J) .GT. ebc) IMIN = IMIN+1
         IF(IMIN .GT. nib) GOTO 1001
C
         IF(ISO .LE. 10) GOTO 500
         A2 = A2P(ISO-10,N-10)
         B2 = B2P(ISO-10,N-10)
         DO I = IMIN,IMAX
            HNU=  freqmin*10**(freqdel*(i-0.5)) 
            X  = EAJ(NO,J)/HNU
C           IF(X.GT.1.0) GOTO 200
            SG = A2*X*X + B2*X*X*X
            PHOT1 = SG * FACT1 *  ABIN(I) / HNU
            PA1 = PA1 + PHOT1
            HE = HE + PHOT1 * HNU * FACT
  200       CONTINUE
         ENDDO
  500    CONTINUE
C
         IMIN = MAX(IK+1,1)
         IF(IMIN .GT. nib) GOTO 1001 
         SSS = 1.66 * EKS ** (0.071)
         DO I=IMIN, nib
            HNU= freqmin*10**(freqdel*(i-0.5)) 
            X  = EKS / HNU
C           IF(X.GT.1.0) GOTO 600
            SG = (295. / EKS) * X ** SSS
            PHOT1 = SG * FACT1 * ABIN(I) / HNU
            PA2 = PA2 + PHOT1
            HE = HE + PHOT1 * HNU * FACT
  600       CONTINUE
         ENDDO
 1001    CONTINUE
C
         HEATR(NO,J) = HEATR(NO,J) + HE
         AP1(NO,J) = PA1
         AP2(NO,J) = PA2
      ENDDO
C
      RETURN
      END
C
C
C
C     NAME=PHTR
      SUBROUTINE PHTR(NO,N,PHOTR,HEATR)
      INCLUDE "MTLPARAM.h"
C
      REAL PHOTR(12,30),HEATR(12,30)
C
      DIMENSION INNER(30),AINNER(12,9),BINNER(12,9),EINNER(12,9)
      DATA INNER/4*0,1,2,3,5*0,4,5,6,3*0,7,8,9,9*0/
      DATA EINNER/0.,30.9,55.8,87.6,172.,283.,422.,589.,
     .  784.,1006.,1842.,2178.,
     .  0.,15.6,36.7,63.8,139.,241.,371.,528.,713.,925.,1731.,2056.,
     .  2*0.,20.3,42.6,108.,201.,321.,469.,644.,847.,1622.,1938.,
     .  6*0.,22.9,57.6,105.2,165.,421.,531.,
     .  6*0.,13.5,43.8,87.6,144.,388.,493.,
     .  7*0.,30.7,70.4,123.,356.,458.,
     .  9*0.,28.,213.,296.,
     .  9*0.,38.,190.,271.,
     .  10*0.,169.,246./
      DATA AINNER/0.,2.357,2.745,1.443,.6024,.3391,.2058,.1385,.09804,
     .  .0764,.04443,.0376,
     .  0.,18.26,6.950,2.817,.9583,.5001,.2672,.1685,.1184,.0913,.04837,
     .  .0407,
     . 2*0.,27.80,8.139,1.630,.7415,.3621,.2394,.1478,.112,.05442,.0455,
     .  6*0.,4.961,2.695,1.401,.8416,.2929,.2254,
     .  6*0.,10.489,5.507,1.575,1.409,.4378,.3263,
     .  7*0.,12.06,4.726,2.351,.6321,.4799,
     .  9*0.,29.6,2.452,1.788,
     .  9*0.,53.8,3.265,2.070,   10*0.,4.078,2.351/
      DATA BINNER/0.,-.7352,.0201,.4321,.4094,.2605,.1935,.1423,.111,
     .  .0865,.04050,.0343,
     .  0.,-6.256,-.2887,1.463,1.029,.5726,.4080,.2885,.2132,.164,
     .  .07817,.0654,
     .  2*0.,-14.93,1.194,2.298,1.106,.7338,.4742,.3510,.267,.1244,.104,
     .  6*0.,-4.607,-1.863,-.6435,-.2594,-.0225,-.0111,
     .  6*0.,-10.412,-4.613,-1.118,-.4648,-.0221,.0105,
     .  7*0.,-11.09,-3.357,-1.014,-.0255,-.0078,
     .  9*0.,-28.3,-.386,-.360,
     .  9*0.,-45.8,-.610,-.178,  10*0.,-.834,.0037/
C
C
C
C
C     INITIALIZE PHOTOIONIZATION RATE AND HEATING TO ZERO
C     FIRST PART, DETERMINE THE PHOTOIONIZATION RATES
C
C     FACT1= 1.E-21/1.6E-12*1.E-18
C     FIRST FACTOR: 1.E-21 DUE TO THE ADOPTED UNITS OF THE PHOTOIONIZATION
C                   FIELD: 1.E-21 ERG/CM^2/SEC
C     SECOND FACTOR: CONVERTING eV to ERG
C     THIRD FACTOR: UNITS OF THE PHOTOIONIZATION CROSS SECTION
C
      FACT1= 1.E-9/1.6*1.E-18
C
C
C
C     FACT = 1.6E-12*1.E23 
C     FIRST FACTOR: CONVERTING BACK ERG to eV      
C     SECOND FACTOR: THE UNITS OF THE PHOTOHEATING IS OUTPUTED IN 1.E-23 ERG/SEC
      FACT = 1.6E11
C
C
C
      DO J=1,N
C        S2,S3,S4,S5: COEFFICIENTS FOR PHOTOIONZATION CROSS SECTION
         SS2 = S2J(NO,J)
         SS3 = S3J(NO,J)
         SS4 = S4J(NO,J)
         SS5 = S5J(NO,J)
C
c
         PHOT = 0.
         HE   = 0.
C
         IMIN = log10(EJ(NO,J)/freqmin+1.e-6)/freqdel + 1.0
         IMIN = max(1,IMIN)
         ebc  = freqmin*10**(freqdel*(imin-0.5)) 
         IF(EJ(NO,J) .GT. ebc) IMIN = IMIN + 1
         IF(IMIN .GT. nib) GOTO 1000 
         ISO  = N-J+1
         IA   = log10(EAJ(NO,J)/freqmin+1.e-6)/freqdel
         IA   = max(1,IA)
         IMAX = nib 
         IF(ISO .GE. 3) IMAX = MIN(nib,IA)
         IISO = INNER(ISO)
         IF(IISO .EQ. 0) GOTO 991
         EIN = EINNER(NO,IISO)
         IIN = log10(EIN/freqmin+1.e-6)/freqdel+1
         IIN = max(1,IIN)
         ebd = freqmin*10**(freqdel*(iin-0.5)) 
         IF(EIN .GT. ebd) IIN = IIN + 1
         IMAX = MIN(IIN-1,Nib)
         IMAX = MAX(1,IMAX) 
  991    CONTINUE
C
         DO I = IMIN,IMAX
            HNU   = freqmin*10.0**(freqdel*(i-0.5)) 
            X     = EJ(NO,J)/HNU
C           IF(X.GT.1.0) GOTO 20
            X2    = X*X
            SG    = SS2*X2 + SS3*X2*X + SS4*X2*X2 + SS5*X2*X2*X
            PHOT1 = SG*FACT1*ABIN(I)/HNU
            PHOT  = PHOT + PHOT1
            HE    = HE + PHOT1 * HNU * FACT
   20       CONTINUE
         ENDDO
C
C        INNER SUBSHELL
         IF(IISO .EQ. 0) GOTO 980
         IF(IIN .GT. nib) GOTO 980
         IMAX = MIN(IA,nib)
         SS2 = AINNER(NO,IISO)
         SS3 = BINNER(NO,IISO)
C
C
         DO I=IIN,IMAX
            HNU=  freqmin*10**(freqdel*(i-0.5)) 
            X = EIN / HNU
C           IF(X.GT.1.0) GOTO 979 
            SG = SS2 * X*X + SS3 * X*X*X
            PHOT1 = SG*FACT1*ABIN(I)/HNU
            PHOT = PHOT + PHOT1
            HE = HE + PHOT1 * HNU * FACT
  979       CONTINUE
         ENDDO
  980    CONTINUE
C
 1000    CONTINUE
         PHOTR(NO,J)=PHOT
         HEATR(NO,J)=HE
      ENDDO
      RETURN
      END
C
C
C
C     NAME=FELINE
      SUBROUTINE FELINE(T,PHOTR,AP1,AP2,CNC,PCOOL)
C
C     T  : TEMPERATURE
C
      INCLUDE "MTLPARAM.h"
C
      REAL E(30),EA(30),S2(30),S3(30),S4(30),S5(30)
     .    ,WAVE(MAXLN),E3(MAXLN),F(MAXLN)
      INTEGER LLL(30)
C
      REAL PHOTR(12,30),AP1(12,30),AP2(12,30)
      REAL CNC(12,30,NID)
      REAL PCOOL(NID)
      REAL PCL(12,NID) 
      REAL C2J(12,30)
      REAL REDUC(4)
      REAL TU(NID)
C
C
C
      DO NO=1,NUME
         DO ID=1,NID 
            PCL(NO,ID)= 0.
         ENDDO
      ENDDO
C
C
C     SOMEHOW THIS STATEMENT PREVENTS THE CODE FROM HAVING SIGMENTATION FAULT
c     TOLD = T
      ST   = SQRT(T)
C
C
C     ELEMENT LOOP
      DO NO=2,NUME
         DO ID=1,NID
            TU(ID) = 0.
         ENDDO
C
         ABUND = 10**(ABUNJ(NO)-12.)
C
         DO L = 1,MAXLN
C           WAVELENGTH
            WAVE(L) = WJ(NO,L)
C           OSCILLATOR STRENGTH
            F(L)    = FJ(NO,L)
C           LINE EXCITATION ENERGY 
            E3(L)   = E3J(NO,L)
         ENDDO
C
         N = NJ(NO)
         NN = N+1
         DO J=1,N
            E(J) = EJ(NO,J)
            EA(J) = EAJ(NO,J)
            S2(J) = S2J(NO,J)
            LLL(J)= LLJ(NO,J)
            S3(J) = S3J(NO,J)
            S4(J) = S4J(NO,J)
            S5(J) = S5J(NO,J)
         ENDDO
C
         E(N+1) = 0.
         LLL(N+1) = 0.
C
C        CALCULATE THE IONIZATION FRACTION.
C        NO: RANK [1...12]
C        N:  ATOMIC NUMBER (e.g., HE:2, C:6, etc)
         CALL NQUIL(T,NO,N,PHOTR,AP1,AP2,CNC,PCOOL,C2J,REDUC,PM)
C
C        CALCULATE THE PERMITTED LINE SPECTRUM AND COOLING 
         IF(N.EQ.2) GOTO 140 
         DENO = AN(NID2)
         IX=0
         DO IL=1,N
            IY=3*LLL(IL)
            IF(IY .LE. 0.) GOTO 270
            PW = 0.
C
            DO L = 1,IY
               JBB=IX+L
               IF(WAVE(JBB) .LE. 0.) GOTO 260
C
               G = GAUNT(T,E3(JBB),N,IL,L,DENO,REDUC,PM)
C
               ALPHAD=ALPHADI(NO,N,IL,L,JBB,T)
               PW = ABUND*1.0E23*(ALPHAD*E3(JBB)*1.602E-12
     .         +2.72E-15*F(JBB)*EXP(-MIN(50.0,11590*E3(JBB)/T))*G/ST)
     .         *12398./(WAVE(JBB)*E3(JBB))
C
               DO ID = 1, NID
                  CN = CNC(NO,IL,ID)
                  PCOOL(ID) = PCOOL(ID) + PW*CN
               ENDDO
  260          CONTINUE
            ENDDO
C
C
            IF(IL .EQ. N) THEN 
              EN = N
              Y  = E3(IX+1) * 11590. / T
              CC = LOG((Y+1)/Y) - 0.4/(Y+1)**2
              GBAR = 0.047
C             GBAR = 0.20 * Y * CC + 0.276 * CC
              IF(N .EQ. 1) GBAR = .017 + .1*CC
C             1.2 FOR CASCADES : TU IS DIRECT EXCITATION TWO PHOTON
              OM   =1.2 * 14.5 * .416 * GBAR * 13.6 / E3(IX+1)
              P2PPH=(8.63E-6*OM/SQRT(T))*EXP(-MIN(50.,Y))
     .             *E3(IX+1)*1.6E11*ABUND
              DO ID=1,NID
                 CN=CNC(NO,IL,ID)
                 TU(ID) = TU(ID) + P2PPH*CN
              ENDDO
            ENDIF
            IF(IL .EQ. N-1) THEN 
              EN = N
C             EXCITATION OF 1S2S 1S LEVEL PRADHAN ; 
C             BALUJA,CALLAWAY,HENRY: 1.2 FOR CASCADES
              Y  = 11590. * E3(IX+1) / T
              CC = LOG((Y+1.)/Y) - 0.4/(Y+1)**2
              GBAR = (.06168+.0002522*N)+(-.02404-.001948*N)*Y*CC
     .             + (-.007604+.002416*N)*(Y-Y*Y*CC)
C             GBAR = (.03014+.0009855*N)+(.2392-.007019*N)*Y*CC
C    .             +(-.1292+.003543*N)*(Y-Y*Y*CC)
              IF(N .EQ. 2) GBAR = .048*(T/10000.)**(.201)
              OM = 1.2 * 14.5 * F(IX+1) * GBAR * 13.6 / E3(IX+1)
              P2PPHE = (8.63E-6*OM/SQRT(T) )*EXP(-MIN(50.,Y))*ABUND
     .                *E3(IX+1)*1.6E11
C
              DO ID = 1, NID
                 CN = CNC(NO,IL,ID)
                 TU(ID) = TU(ID) + P2PPHE*CN
              ENDDO
            ENDIF
C
C
C
  270       CONTINUE 
            IX=IX+IY
         ENDDO
C
         DO ID = 1, NID
            PCL(NO,ID) = PCOOL(ID) + TU(ID)
         ENDDO
C
  140    CONTINUE
      ENDDO
C
C 
      DO ID = 1, NID
         PCOOL(ID) = 0.
         DO NO=2,NUME
            NN=NJ(NO)
            RECENER=0.
            DO JK=2,NN+1
               RECENER=RECENER+1.6E11*10**(ABUNJ(NO)-12.)
     .                 *CNC(NO,JK,ID)*C2J(NO,JK)*EJ(NO,JK-1)
C              EJ=E: READIN CONSTANT          
C              C2J: RECOMBINATION COEFFICIENT FROM nquil.f
            ENDDO
            PCL(NO,ID)=PCL(NO,ID)-RECENER
            IF(PCL(NO,ID).LE.0) PCL(NO,ID)=0.
         ENDDO
C
         DO NO=2,NUME
C           PCL() IS LOCAL TO THIS SUBROUTINE
            PCOOL(ID) = PCOOL(ID)+PCL(NO,ID)
         ENDDO
      ENDDO
C
      RETURN
      END
C
C
C     NAME NQUIL
      SUBROUTINE NQUIL(T,NO,N,PHOTR,AP1,AP2,CNC,PCOOL,C2J,REDUC,PM)
C     THIS ROUTINE COMPUTES THE CONCENTRATIONS OF ELEMENTS AT ALL LEVELS,
C     CNC(), AND PERMITTED LINE COOLING (PCOOL)
C
C     T:  TEMPERATURE
C     NO: RANK OF THE METAL ([1...12])
C     N:  ATOMIC NUMBER OF THE METAL
C
C     AS INPUT THIS ROUTINE TAKES THE ELECTRON NUMBER DENSITY IN THE FORM
C     AN(), IN UNITS OF THE HYDROGEN ATOM NUMBER DENSITY
C
      INCLUDE "MTLPARAM.h"
C
      REAL PHOTR(12,30),AP1(12,30),AP2(12,30)
      REAL CNC(12,30,NID)
      REAL PCOOL(NID)
      REAL C2J(12,30)
      REAL CONCE(30)
      REAL C1(30),C2(30),C1D(30,NID),CA(30,NID)
      REAL REDUC(4)
C
C
      DIMENSION EMETJ(30,4),CD1(28),CD2(28)
      DIMENSION RAT(30),PROD(30),C5(30),C6(30),OVER(30),IMETA(30)
      DIMENSION SEC(30),C0(30)
      DATA EMETJ/5*0.,6.5,8.35,10.2,0.,14.,0.,17.6,0.,21.6,0.,25.3,0.,
     .  29.0,0.,32.7,5*0.,47.0,0.,51.8,0.,0.,
     .  5*0.,5.34,7.10,8.87,0.,12.5,0.,16.0,0.,19.8,0.,24.0,0.,28.6,0.,
     .  34.2,5*0.,57.1,0.,64.8,0.,0.,
     .  5*0.,4.18,5.78,7.48,0.,11.0,0.,14.7,0.,18.6,0.,23.0,0.,28.0,0.,
     .  33.9,5*0.,60.4,0.,69.2,0.,0.,  11*0.,
     .  2.72,0.,6.55,0.,10.4,0.,14.2,0.,18.1,5*0.,29.7,0.,38.9,2*0./
C
      DATA IMETA/0,0,0,1,2,3,5*0,4,18*0/
      DATA CD1/.0024,-.0005,.00,.061,.027,.011,.005,.005,0.,.0107,
     .  .09,.13,.11,.081,.075,.066,.051,11*0./
      DATA CD2/0.,.01485,.30,.108,.107,.024,.075,.051,.054,.0167,.36,
     .  17*0./
      DATA SEC/4*1.,.9,.9,6*1.,.6,.5,1.,.7,.9,13*1./
C
C
C
      DENO  = AN(NID2)
      IX    = 0
      ABHE  = 10.0**(ABUNJ(1) - 12.)
      ABUND = 10.0**(ABUNJ(NO)-12.)

C
      PM    = 0.
C
      ORANGE  = 1.
      C5(1)   = 0.
      C5(N+1) = 0.
      JLOW    = 1
      JMAX    = N+1
      JTOP    = N+1
      C6(1)   = 0
      T6      = T*1.0E-6
      ST      = SQRT(T)
C
C     LOOP THROUGH ALL THE IONS
      DO J=1,N
C        COLLISIONAL IONIZATION RATE FROM YOUNGER, CALCULATED IN CION
C        E(J): IONIZATION POTENTIAL
         EION  = EJ(NO,J)
         C1(J) = CION(N,J,EION,T)
C
         IF(N-J .LE. 1) GOTO 104
C        IONIZATION FROM METASTABLE STATES
         IMET  = IMETA(N-J+1)
         IF(IMET .EQ. 0) GOTO 105
         EMET  = EMETJ(N,IMET)
C        POPMET: POPULATION OF METASTABLE LEVEL
         FMET  = POPMET(N,IMET,EMET,DENO,T,REDUC)
         IF(IMET.EQ.1) PM = FMET
         EM    = EJ(NO,J) - EMET
         C1(J) = C1(J)*(1.-FMET) + FMET * CION(N,J,EM,T)
  105    CONTINUE
C        INNERSHELL IONIZATION
C        EAJ(N,J): INNERSHELL IONIZATION POTENTIAL
         EION  = EAJ(NO,J)
         JP    = N-1
         IF(N-J .GE. 10) JP = N-9
         C1(J) = C1(J) + CION(N,JP,EION,T)
  104    CONTINUE
C        IONIZATION FOLLOWING INNERSHELL EXCITATION: COWAN AND MANN
         C1(J) = C1(J) + AUTOIN(N,J,T)
C  
C
C
C        NOW GO ON TO CALCULATING THE RECOMBINTION RATES.
C        RADIATIVE RECOMBINATION
C
         ZEFF  = J
         BETA  = ZEFF*ZEFF/(6.34*T6)
C
C        C2 IS THE RECOMBINATION RATE TO FORM HYDROGENIC IONS,
C        SEE COX AND TUCKER, APJ 157:1159. 
         C2(J+1)= 2.07E-11*ZEFF*ZEFF*(0.4288 + 0.5*LOG(BETA)
     .           +0.469*(BETA**(-1./3.)))/SQRT(T)
         RECOM1 = C2(J+1)
         PLUM   = 0
         RECN   = 0.
         RGND   = 0.
C    
         IF(N .EQ. J) GOTO 50
         NVAL   = GND(N-J) + 0.1
         STARN  = SQRT(13.6/EJ(NO,J))*ZEFF
C  
         DO NQM=1,NVAL
            QM      = NQM
            BE      = BETA/(QM*QM)
            CHECK   = 5.197E-14*J*SQRT(BE)*SEATON(BE,QM)
            C2(J+1) = C2(J+1) - 5.197E-14*J*SQRT(BE)*SEATON(BE,QM)
         ENDDO
C
         IF(C2(J+1) .LT. 0)  C2(J+1) =0
         RECOM2 = C2(J+1)
         I      = N-J
C
C        RECOMBINATION TO GROUND STATE
C        RGND= RECOMBINATION TO GROUND
         IF(I .LE. 17) RGND = GREC(NO,N,J,EJ(NO,J),T)
         IF(I .LE. 1) GOTO 50
         IF(I .GE. 4 .AND. I .LE. 9) GOTO 50
         EI     = EJ(NO,J)
         IF(I .GE. 17) GOTO 40
C
C        RECOMBINATION TO OTHER STATES IN SAME SHELL
         LION   = 1
         IF(I .EQ. 12) LION = 3
         EI     = EI - E3J(NO,IX+LION)
   40    CONTINUE
         STARN  = SQRT(13.595/EI)*ZEFF
         BE     = BETA / (STARN*STARN)
         RECN   = EFFNEW(N-J,T,ZEFF,N)*SEATON(BE,STARN)*2.599E-14*J
     .           *SQRT(EI*11590./T)/(STARN*STARN)
   50    CONTINUE
         C2(J+1)= C2(J+1) + RECN + RGND
C  
         IF(BETA .GE. 0.5) THEN 
           CHI  = 0.5 * (0.735 + LOG(BETA) + 1./(3.*BETA))
         ELSE
           CHI  = 1.2*BETA+(0.55+1.04*LOG(BETA))*BETA*BETA
     .           +(-0.43+1.01*LOG(BETA))*BETA**3
         ENDIF
C 
         DO NQM=1,NVAL
            QM  = NQM
            CHI = CHI - P(BETA,QM)
         ENDDO
C
         CHI    = CHI + EFFN(N-J,ZEFF,T) * P(BETA,STARN)
     .                 /(2.*STARN*STARN)
         C6(J+1)= EJ(NO,J)*1.6027E-12*(C2(J+1)
     .          + 4.22E-22*SQRT(T/1.037E-11)*CHI)
C  
C        DIELECTRONIC RECOMBINATION ***** DENSITY DEPENDENCE 
C        AND SECONDARY  ***** AUTOIONIZATION
         IY     = LLJ(NO,J) * 3
         IF(J .EQ. 1) GOTO 370
         ZEFF   = J-1
         RHO    = (DENO/ZEFF**7.0)**0.2
         DD     = 1.
         IF(RHO.GE.3.0) DD=1./(1.+(CD1(N-J+1)+CD2(N-J+1)/ZEFF)*RHO)
C
         DO L=1,IY
            LN   = IX+L
            PLUM = PLUM + ALPHADI(NO,N,J,L,LN,T)*MIN(DD,SEC(N-J+1))
         ENDDO
C
C        DIELECTRONIC RECOMBINATION FOR METASTABLE LEVELS OF BE,B,C,MG ISO
         IMET   = IMETA(N-J+1)
         IF(IMET .EQ. 0) GOTO 321
         ZF     = J-1
         B      = SQRT(ZF)*(ZF+1.)**2.5/SQRT(ZF*ZF+13.4)
         PLUM   = PLUM * (1.-FMET)
         PLUM   = PLUM + FMET*DIMET(N,J,T,B,DD)
  321    CONTINUE
C
C        ADD DIELECTRONIC RECOMB AND STOREYS LOW T DIELEC REC
         C2(J)  = C2(J) + PLUM +ALFLO(N,J,T)
  370    CONTINUE
         OVER(J)= PLUM/ORANGE
         IF(J .NE. 1) C5(J) = PLUM*EJ(NO,J-1)*ABUND*1.6027E11
         IX     = IX + IY
         ORANGE = RECOM1
      ENDDO
C
      C1(N+1)= 0.
      C2(1)  = 0.
      NN     = N + 1
C
C     CHARGE TRANSFER CONTRIBUTION
      FPLUS  = 1.0
      HENEUT = CNC(1,1,NID2)  
      HEPLUS = CNC(1,2,NID2)  
      DENH   = 1./(.003 +FPLUS+ABHE*HEPLUS+ABHE*2.*(1.-HENEUT-HEPLUS))
C  
C
C     S. BUTLER CHARGE TRANSFER RATES
      DO J=1,N
         C1(J)  = C1(J)+CTHI(N,J,T)*DENH*FPLUS
     .          + CTHEI(N,J,T,DENO)*HEPLUS*ABHE
         C2(J+1)= C2(J+1)+CTHR(N,J+1,T,PM)*DENH*(1.-FPLUS)  
     .          + CTHER(N,J+1,T)*HENEUT*ABHE
         C0(J)  = C1(J)
C
C        ADDED 11/21/96, R.CEN
         DO ID=1,NID
            C1D(J,ID) = C1(J)
         ENDDO
      ENDDO
C
C
C
C     PHOTOIONIZATION
      IF(IPHOT.GT.0) THEN
        DO J=1,N
C          DENSITY DEPENDENCE
           DO ID=1,NID
C             PHOTR(NO,J): PHOTOIONIZATION RATE  FROM phtr.f
C             CHANGED 11/21/96, R.CEN
              C1D(J,ID) = C1D(J,ID) + PHOTR(NO,J)/(AN(ID))
           ENDDO
        ENDDO
C
        DO J=1,30
           DO ID=1,NID
              CA(J,ID) = 0.
           ENDDO
        ENDDO
        RN=N
        FLUOR = 0.5 * (RN/26.)**4
        DO J=1,N-2
           ISO = N-J+1
C          AP1(NO,J), AP2(NO,J): PHOTOIONIZATION RATE  FROM aphtr.f
           PA1 = AP1(NO,J) 
           PA2 = AP2(NO,J)
           DO ID=1,NID
              CA(J,ID)  = (PA1+PA2)*(1.-FLUOR)/(AN(ID))
              C1D(J,ID) = C1D(J,ID)+(PA1+PA2)*FLUOR
     .                             /(AN(ID))
              IF(ISO .EQ. 3) THEN 
                CA(J,ID)  = 0.
                C1D(J,ID) = C1D(J,ID) + PA2/(AN(ID))
              ENDIF
              IF(ISO .EQ. 11) THEN 
                CA(J,ID)  = (1.-FLUOR)*PA2/(AN(ID))
                C1D(J,ID) = C1D(J,ID)+(PA2*FLUOR+PA1)
     .                               /(AN(ID))
              ENDIF
           ENDDO
        ENDDO
C
        DO ID=1,NID
           Q = 0.
           X = CA(1,ID)
           DO J=1,N-1
              C1D(J,ID) = C1D(J,ID) + X
              IF(C1D(J,ID) .GT. 0.) THEN 
                Q = C2(J+1) * CA(J,ID) / C1D(J,ID)
                X = CA(J+1,ID) + Q
              ENDIF
           ENDDO
        ENDDO
      ENDIF 
C     END OF PHOTOIONIZATION SECTION
C
C
C
C
      DO ID=1,NID
         DO J=1,N
            C2(J+1)=ABS(C2(J+1))
            IF(C1D(J,ID) / C2(J+1) .LE. 0.) THEN
              RAT(J)=1.0E+6
            ELSE
              RAT(J)=C2(J+1)/C1D(J,ID)
            ENDIF
            IF(J.EQ.1 .AND. T.GE.1.0E8) RAT(J)=MIN(RAT(J),1.0E4)
            IF(RAT(J) .GE. 1.0E-5) GOTO 34
            JLOW=J+1
   34       CONTINUE
            IF(RAT(J) .LT. 1.0E+5) GOTO 400
            JMAX=J
            GOTO 41
  400       CONTINUE
         ENDDO
C
         IF(JLOW .EQ. JMAX) GOTO 525
   41    CONTINUE
         JUMP=JMAX - 1
         PROD(JUMP)=RAT(JUMP)
         KTOP=JUMP - JLOW
         SUM = 1.00000 + PROD(JUMP)
         IF(KTOP .EQ. 0)  GOTO 550
C
         DO K=1,KTOP
            PROD(JUMP-K)=RAT(JUMP-K)*PROD(JMAX-K)
            PROD(JUMP-K) = MIN(PROD(JUMP-K),1.0E30)
            SUM = SUM + PROD(JUMP-K)
         ENDDO
C
         GOTO 550
  525    CONTINUE 
         SUM = 1.0000
  550    CONTINUE
C
         DO J=1,JTOP
            CONCE(J)=0
            C5(J) = 0.
         ENDDO
C
         CONCE(JMAX) = 1.000/SUM
         KMAX=JMAX-JLOW
         IF(KMAX .EQ. 0)  GOTO 725
C
         DO K=1,KMAX
            CONCE(JMAX-K)=PROD(JMAX-K)*CONCE(JMAX)
         ENDDO
C
  725    CONTINUE
C
         PCOOL(ID) = 0.
         DO J=1,NN
            CNC(NO,J,ID) = CONCE(J)
            PCOOL(ID)=PCOOL(ID) + C6(J)*CNC(NO,J,ID)*ABUND*1.0E23
         ENDDO
C
      ENDDO
C
      DO J=1,N +1
         C2J(NO,J) = C2(J)
      ENDDO
C
      RETURN
  100 FORMAT(2I5,E112.4)
      END
C
C
C     NAME=P
      FUNCTION P(Y,Q)
      A = Y/(Q*Q)
      IF(A .LE. 75) GO TO 20
      P = 1./Q
      RETURN
   20 CONTINUE
      P = A * (1.-EXINT1(A,3)) / Q
      RETURN
      END
C
C
C
C     NAME=ALFLO
      FUNCTION ALFLO(N,J,T)
C  LOW TEMPERATURE DIELECTRONIC RECOMBINATION FROM NUSSBAUMER AND
C  STOREY: C - SI FOR 1000 - 60,000 K
C
      DIMENSION IION(8,16),A(25),B(25),C(25),D(25),F(25)
C
      DATA IION/40*0, 0,1,2,3,0,0,0,0,  0,4,5,6,7,0,0,0,
     .  0,8,9,10,11,12,0,0,  8*0,  0,13,14,15,16,17,18,19,
     . 8*0, 0,20,6*0, 0,21,22,5*0,  0,23,24,25,4*0,  8*0, 8*0/
C
      DATA A/.0108,1.8267,2.3196,  0.0,0.0320,-.8806,.4134,
     .  0.0,-.0036,0.0,.0061,-2.8425,
     .  0.0,.0129,3.6781,-.0254,-.0141,19.928,5.4751,
     .  1.2044,.0219,.7086,  -.0219,3.2163,.1203/
C
      DATA B/-.1075,4.1012,10.7328,  .6310,-.6624,11.2406,-4.6319,
     .  .0238,.7519,21.879,.2269,.2283,
     .  0.0,-.1779,14.1481,5.5365,33.8479,235.0536,203.9751,
     .  -4.6836,-.4528,-3.1083,.4364,-12.0571,-2.690/
C
      DATA C/.2810,4.8443,6.8830,  .1990,4.3191,30.7066,25.9172,
     .  .0659,1.5252,16.2730,32.1419,40.4072,
     .  0.,.9353,17.1175,17.0727,43.1608,152.5096,86.9016,
     .  7.6620,2.5427,7.0422,.0684,16.2118,19.1943/
C
      DATA D/-.0193,.2261,-.1824,  -.0197,.0003,-1.1721,-2.2290,
     .  .0349,-.0838,-.7020,1.9939,-3.4956,
     .  0.,-.0682,-.5017,-.7225,-1.6072,9.1413,-7.4568,
     .  -.5930,-.1678,.5998,-.0032,-.5886,-.1479/
C
      DATA F/-.1127,.5960,.4101,.4398,.5946,.6127,.2360,
     .  .5334,.2769,1.1899,-.0646,1.7558,
     .  0.,.4516,.2313,.1702,.1942,.1282,2.5145,
     .  1.6260,.2276,.4194,.1342,.5613,.1118/
C
      ALFLO = 0.
      IF(J .EQ. 1 .OR. J .GT. 8) RETURN
      IF(N .GT. 14) RETURN
C
      T4 = T*0.0001
C
      IF(T4 .LT. 0.1 .OR. T4 .GT. 6.) RETURN
C
      IF(N .EQ. 6 .AND. J .EQ. 2 .AND. T4 .LT. 0.2) RETURN
      IF(N .EQ. 8 .AND. J .EQ. 4 .AND. T4 .LT. 0.2) RETURN
      IF(N .EQ. 8 .AND. J .EQ. 6 .AND. T4 .LT. 0.4) RETURN
      IF(N .EQ. 10 .AND. J .EQ. 8 .AND. T4 .LT. 0.2) RETURN
      IF(N .EQ. 12 .AND. J .EQ. 2 .AND. T4 .LT. 0.15) RETURN
      IF(N .EQ. 14 .AND. J .EQ. 3 .AND. T4 .LT. 0.15) RETURN
C
      IJ = IION(J,N)
      IF(IJ .EQ. 0) RETURN
C
      ALFLO = 1.0E-12*(A(IJ)/T4 + B(IJ) + C(IJ)*T4 + D(IJ)*T4*T4)
     .               *EXP(-MIN(50.0,F(IJ)/T4))/T4**1.5
C
      RETURN
      END
C
C
C     NAME=ALPHADI
      FUNCTION ALPHADI(NO,N,J,L,LN,T)
C     DIELECTRONIC RECOMBINATION : BURGESS AND TWORKOWSKI FOR H-LIKE
C     N: ATOMIC NUMBER OF THE METAL
C     J: [1...N]
C     L: [1...3*LLJ(NO,J)]
C     LN: IX+IL (IL=[1...L]), (IX=0, IX=IX+3*LLJ(NO,J))
C     T: TEMPERATURE
C
      INCLUDE "MTLPARAM.h"
C
C     BURGESS AND TWORKOWSKI
      Z12 = J-12.
      TWORK = .84 + .5/J**2 + .03*Z12/(1.+4.5E-5*Z12**3)
      IF(N .NE. J) TWORK = 1.
      ALPHADI = 0.
      DL = DELTT(N,J,L)
      IF(DL .LE. 0.) RETURN
      ZF = J-1.
C     B = IS FOUND IN EQ 4 APJ 141:1590
      B = SQRT(ZF)*(ZF+1.)**2.5/SQRT(ZF*ZF+13.4)
C     X = DOES NOT LOOK THE SAME AS IN THE DEFINITION OF X IN APJ 141:1590
      X = E3J(NO,LN) / ((ZF+1.)*13.6)
C     A = IS FOUND IN EQ 3 APJ 141:1590
      A = SQRT(X) / (1.+.105*X + .015*X*X)
C     EBAR = IS ACTUALLY DEFINED AS E3/A WHERE A IS IN EQ 6 APJ 141:1590
      EBAR = E3J(NO,LN) / (1.+.015*ZF**3/(ZF+1.)**2)
      IF(N-J .NE. 1) GOTO 25
      B = (0.5 * ZF/SQRT(ZF)) * B
      X = 0.75 * J
      A = SQRT(X) * 0.5 / (1.+0.21*X+0.03*X*X)
C     YOUNGER JQSRT 29, 67 FOR HE-LIKE IONS
   25 CONTINUE
      ALPHADI = 0.0030*T**(-1.5)*A*B*FJ(NO,LN)*TWORK*DL
     .                *EXP(-MIN(50.0,EBAR*11590./T))
      RETURN
      END
C
C
C
C
C     NAME=ATREAD
      SUBROUTINE ATREAD
      INCLUDE "MTLPARAM.h"   

      REAL E(30),EA(30),S2(30),S3(30),S4(30),S5(30)
     .    ,WAVE(MAXLN),E3(MAXLN),F(MAXLN)
      INTEGER LLL(30)
      INTEGER INPUT_DONE

      integer io_err

      DATA INPUT_DONE / 0 /

      save input_done

      if ( input_done < 1 ) then

      write(57,'("Read atomic data")')

      open(unit=18, file='ATOMIC.DAT', status='old', form='formatted',
     &     iostat=io_err)

      if ( io_err .ne. 0 ) then
        write(0,'("Error opening ATOMIC.DAT")')
        ERROR_MESSAGE
      end if

      rewind(18)

!     OPEN(UNIT=18,FILE='ATOMIC.DAT')
!     REWIND(18)

      DO NO=1,NUME
!        E(J) : IONIZATION POTENTIAL
!        EA(J): INNERSHELL IONIZATION POTENTIAL
         IX = 0
   25    READ(18,26) N,J,E(J),EA(J),S2(J),S3(J),S4(J),S5(J),LLL(J)
         IY = LLL(J)
         IF(IY .LE. 0) GOTO 30

         DO L=1,IY
            IAA = IX + 3*L-3
            READ (18,29) (WAVE(IAA+K),E3(IAA+K),F(IAA+K),K=1,3)
!           WAVELENGTH, LINE EXCITATION ENERGY OSCILLATOR STRENGTH
         ENDDO
 
         IX = IX + 3 * IY
   30    IF(N .GT. J) GOTO 25
         NJ(NO) = N

         DO J=1,N
            EJ (NO,J) = E(J)
            EAJ(NO,J) = EA(J)
            S2J(NO,J) = S2(J)
            S3J(NO,J) = S3(J)
            S4J(NO,J) = S4(J)
            S5J(NO,J) = S5(J)
            LLJ(NO,J) = LLL(J)
         ENDDO
         DO L=1,MAXLN
            WJ (NO,L) = WAVE(L)
            E3J(NO,L) = E3(L)
            FJ (NO,L) = F(L)
         ENDDO
!        write(*,*)'NO, = ',NO,IX,IY,N,MAXLN
      ENDDO

      CLOSE(18)

      input_done = 1

      end if

      RETURN

   26 FORMAT(2I5,F5.0,5F7.0,I5)
   29 FORMAT(9F6.0)

      END
C
C
C     NAME=AUTOIN
      FUNCTION AUTOIN(N,J,T)
C     INNERSHELL EXCITATION FOLLOWED BY AUTOIONIZATION FOR NA,MG,AL,SI
C     SEQUENCES. CRANDALL ET AL PRA 25,143 FOR MG,SI.  MANN FOR FE* .75.
C     COMPROMISE BETWEEN COWAN AND MANN AND DOUBLE THAT.
C     USE 2S AT 2P ENERGY, ETC FOR B,C,... SEQUENCES
      DIMENSION ET(20),OM(20)
      DATA ET/0.,55.,0.,115.,0.,190.,0.,280.,0.,375.,5*0.,802.,
     .  0.,1002.,0.,0./
      DATA OM/0.,.34,0.,.16,0.,.18,0.,.20,0.,.22,5*0.,.29,0.,.3,0.,0./
      AUTOIN = 0.
      IF(N .LE.10) RETURN
      I = N-J+1
      IF (I .LE. 10 .OR. I .GE. 15) RETURN
      AUTOIN = 8.63E-6*OM(N-10)*EXP(-11590.*ET(N-10)/T)/SQRT(T)
C     ASSUMES THAT NUMBER OF 32,3P ELECTRONS DOES NOT MATTER
      RETURN
      END
C
C
C     NAME=CION
      FUNCTION CION(N,J,E,T)
C
C     CALCULATES THE COLLISIONAL IONIZATION SOMETHING OR ANOTHER
C     SOME OF THE FORMULAE HERE ARE STRAIGHT OUT OF YOUNGER
C     SM YOUNGER JQSRT 26, 329; 27, 541; 29, 61   WITH MOORES FOR UNDONE
C     A0 FOR B-LIKE ION HAS TWICE 2S PLUS ONE 2P  AS IN SUMMERS ET AL
C     CHI = KT / I
C
      DIMENSION A0(30),A1(30),A2(30),A3(30),B0(30),B1(30),B2(30),
     .  B3(30),C0(30),C1(30),C2(30),C3(30),D0(30),D1(30),D2(30),D3(30)
C
      DATA A0/13.5,27.0,9.07,11.80,20.2,28.6,37.0,45.4,53.8,62.2,
     .11.7,38.8,37.27,46.7,57.4,67.0,77.8,90.1, 106.,120.8,135.6,150.4,
     .165.2,180.0,194.8,209.6,224.4,239.2,154.0,268.8/
      DATA A1/-14.2,-60.1,4.30,27*0./
      DATA A2/40.6,140.,7.69,27*0./
      DATA A3/-17.1,-89.8,-7.53, 27*0./
C
      DATA B0/-4.81,-9.62,-2.47,-3.28,-5.96,-8.64,-11.32,-14.00,
     .  -16.68,-19.36,  -4.29,-16.7,  -14.58,-16.95,-19.93,-23.05,
     .  -26.00,-29.45, -34.25,-38.92,-43.59,-48.26,-52.93,-57.60,
     .  -62.27,-66.94,-71.62,-76.29,-80.96,-85.63/
      DATA B1/9.77,33.1,-3.78, 27*0./
      DATA B2/-28.3,-82.5,-3.59, 27*0./
      DATA B3/11.4,54.6,3.34, 27*0./
C
      DATA C0/1.85,3.69,1.34,1.64,2.31,2.984,3.656,4.328,5.00,5.672,
     .1.061,1.87, 3.26,5.07,6.67,8.10,9.92,11.79,  7.953,8.408,
     .8.863,9.318,9.773,10.228,10.683,11.138,11.593,12.048,12.505,12.96/
      DATA C1/0.,4.32,.343, 27*0./
      DATA C2/0.,-2.527,-2.46, 27*0./
      DATA C3/0.,.262, 1.38, 27*0./
C
      DATA D0/-10.9,-21.7,-5.37,-7.58,-12.66,-17.74,-22.82,-27.9,-32.98,
     .-38.06, -7.34,-28.8, -24.87,-30.5,-37.9,-45.3,-53.8,-64.6,
     .-54.54,-61.70,-68.86,-76.02,-83.18,-90.34,-97.50,-104.66,-111.82,
     .-118.98,-126.14,-133.32/
      DATA D1/8.90,42.5,-12.4, 27*0./
      DATA D2/-35.7,-131.,-8.09, 27*0./
      DATA D3/16.5,87.4,1.23, 27*0./
C
      CION = 0.
      CHIR = T / (11590. * E)
      IF( CHIR .LE. .0115) RETURN
      CHI = MAX(CHIR,0.1)
      CH2 = CHI * CHI
      CH3 = CH2 * CHI
      ALPHA = (.001193 + .9764*CHI + .6604*CH2 + .02590*CH3)
     .  /(1.0 + 1.488*CHI + .2972*CH2 + .004925*CH3)
      BETA = (-.0005725 + .01345*CHI + .8691*CH2 + .03404*CH3)
     .  /(1.0 + 2.197*CHI + .2457*CH2 + .002503*CH3)
      J2 = J*J
      J3 = J2*J
      ISO = N-J+1
C
      A = A0(ISO) +  A1(ISO)/J + A2(ISO)/J2 + A3(ISO)/J3
      B = B0(ISO) +  B1(ISO)/J + B2(ISO)/J2 + B3(ISO)/J3
      C = C0(ISO) +  C1(ISO)/J + C2(ISO)/J2 + C3(ISO)/J3
      D = D0(ISO) +  D1(ISO)/J + D2(ISO)/J2 + D3(ISO)/J3
C
C     FE II EXPERIMENTAL IONIZATION MONTAGUE ET AL: D. NEUFELD FIT
      IF(N .NE. 26 .OR. J .NE. 2) GOTO 10
      A = -13.825
      B = -11.0395
      C = 21.07262
      D = 0.
   10 CONTINUE
      CH = 1./CHI
      FCHI = 0.3*CH*(A+B*(1.+CH) + (C-(A+B*(2.+CH))*CH)*ALPHA 
     .     + D * BETA * CH)
      IF(ISO .GE. 4 .AND. ISO .LE. 10) FCHI = FCHI * 1.59
C     CORRECT YOUNGER JQSRT 27, 541 FROM TABLE TO GRAPHS
      CION = 2.2E-6*SQRT(CHIR)*FCHI*EXP(-MIN(50.0,1./CHIR))
     .      /(E*SQRT(E))
      RETURN
      END

      FUNCTION CTHEI(N,J,T,DENO)
C  SCOTTS RATES FOR HE+  + J  TO HE0  + J+1
C  THROUGH ARGON  0.1 < T4 < 10.
      CTHEI = 0.
      T4 = 0.0001 * T
      IF(J .EQ. 1 .OR. J .GE. 4) RETURN
      IF(J .EQ. 3) GOTO 10
      IF(N .EQ. 6) CTHEI = 5.3E-10 * EXP(-MIN(50.0,13.2/T4))
C  ALBERTS RATES FOR N+ FROM 1D LEVEL
      IF(N .NE. 7) GOTO 70
      QDW = 5.14E-06 * DENO / SQRT(T)
      QUP = 2.86E-6 * EXP(-MIN(50.0,21826./T))
      R21 = QUP / (QDW + .0041)
      D1 = R21 / (1.+R21)
      P3 = 1. - D1
      CTHEI = P3*4.1E-10*EXP(-MIN(50.0,10.4/T4))
     .      + D1*1.1E-10*EXP(-MIN(50.0,3.59/T4))
   70 CONTINUE
      IF(N .EQ. 14) CTHEI = 2.9E-10*EXP(-MIN(50.0,8.7/T4))
      IF(N .EQ. 18) CTHEI = 1.1E-10*EXP(-MIN(50.0,3.57/T4))
      RETURN
   10   CONTINUE
      IF(N .EQ. 14) CTHEI = 3.4E-9*EXP(-MIN(50.0,10.5/T4))
      IF(N .EQ. 16) CTHEI = 4.8E-10*EXP(-MIN(50.0,15.7/T4))
      IF(N .EQ. 26) CTHEI = 1.2E-9*SQRT(T4)*EXP(-MIN(50.0,12./T4))
      RETURN
      END
C
C
C     NAME CTHER
      FUNCTION CTHER(N,J,T)
C     FITS TO SCOTTS RATES FOR  HE I + J   TO   HE II + J-1
C     1000  TO 100,000 K THROUGH ARGON
      DIMENSION A(2,9)
      DATA A/4*0.,.059,.00001,1.1,0.7,.00001,1.7,.075,2.2,.096,1.2,
     .  1.1,.0008,.00001,1.0/
C
      CTHER = 0.
      IF(J .EQ. 1) RETURN
      T4 = MAX(.0001 * T,0.1)
      T4 = MIN(T4,10.)
      CR = 5.4E-10 * T4
      IF(J .EQ. 2 .OR. N .GT. 18) RETURN
      CTHER = MAX(CR,(-.1+.4*J)*1.0E-9)
      IF(J .GE. 6) RETURN
      IF(J .NE. 3) GO TO 10
      CTHER = 0.
C     INCLUDE ALBERTS RATE TO 1D OF N+
      IF(N .EQ. 7) CTHER = 3.0E-10 * (1.+.26*T4) + 6.2E-11
      IF(N .EQ. 8) CTHER = 3.3E-10 * T4 ** .7
      IF(N .EQ. 18) CTHER = 1.3E-10
      RETURN
   10 IF(N .NE. 7) GO TO 11
      CTHER = 1.1E-10
      IF(J .EQ. 5) CTHER = 2.0E-9
      RETURN
   11 CONTINUE
      IN = N/2
      CTHER = A(J-3,IN) * 1.0E-9
      IF(J .EQ. 5) GO TO 12
      IF(N .EQ. 6) CTHER = 5.1E-11 * T4 ** 1.46
      IF(N .EQ. 14) CTHER = 9.6E-10 * T4 ** .55
      IF(N .EQ. 16) CTHER = 1.1E-9 * T4 ** .63
      RETURN
   12 IF(N .EQ. 18) CTHER = 1.0E-9 * T4 ** .91
      RETURN
      END
C
C
C     NAME=CTHI
      FUNCTION CTHI(N,J,T)
C     SCOTTS RATES FOR H II + J TO HI + J+1 THROUGH ARGON
      T4 = MAX(.0001 * T,0.1)
      T4 = MIN(T4,10.)
      CTHI = 0.
      IF(J .GT. 2) RETURN
      IF(J .EQ. 2) GOTO 2
      IF(N .EQ. 6) CTHI = 2.5E-15
      IF(N .EQ. 7) CTHI = 3.3E-13 * T4 ** .12
      IF(N .EQ. 8) CTHI = 3.3E-10
      IF(N .EQ. 16) CTHI = 6.7E-10
      RETURN
    2 IF(N .EQ. 12) CTHI = 2.6E-11 * EXP(-MIN(50.0,6.0/T4))
      IF(N .EQ. 14) CTHI = 9.6E-10 * EXP(-MIN(50.0,2.74/T4))
C     NEUFELDS FE II - FE III RATES
      IF(N .EQ. 26) CTHI = 1.0E-9*1.666*(1.25+.25*LOG10(T4))
     .                           *EXP(-MIN(50.0,3.005/T4))
      RETURN
      END
C
C
C
C     NAME=CTHR
      FUNCTION CTHR(N,J,T,PM)
C     SCOTTS RATES FOR H I + J TO H II + (J-1)
C     THROUGH ARGON : NEUFELTS FE II - FE III
C
      INCLUDE "MTLPARAM.h"
C
      DIMENSION A(30),B(30)
      DATA A/0.,0.,.00001,0.,5.9,6.2,7.8,9.4,11.,13.,14.,16.1,17.6,
     .  19.,20.4,21.8,23.3,13*0./
      DATA B/4*0.,.25,.18,.13,.1,.08,.05,.04,19*0./
C
      CTHR = 0.
      IF(J .EQ. 1) RETURN
      T4 = MAX(T * .0001,0.1)
      T4 = MIN(T4,10.)
C     NEUFELDS FE II - FE III, FE III-FE IV, AND NI II - NI III RATES
      IF(N .EQ. 26 .AND. J .EQ. 3) CTHR=1.0E-9*(1.25+.25*LOG10(T4))
      IF(N .EQ. 26 .AND. J .EQ. 4) CTHR = 3.4E-9 * SQRT(T4)
      IF(N .EQ. 28 .AND. J .EQ. 3) CTHR = 1.0E-9*(.34+1.86*T4)
      IF(N .GT. 18) RETURN
      CTHR = A(J) * (1.+B(J)*T4) * 1.0E-9
      IF(J .GE. 6) RETURN
      IF(N .NE. 2) GO TO 6
      CTHR = 1.9E-15
      IF(J .EQ. 3) CTHR = 1.7E-13
      RETURN
    6 IF(N .NE. 6) GO TO 7
C     ALBERTS C III TRIPLET P
      IF(J .EQ. 3) CTHR=2.5E-9*PM*EXP(-MIN(50.0,15000./T))
      IF(J .EQ. 4) CTHR = 2.9E-9
      IF(J .EQ. 5) CTHR = 7.6E-10 * T4 ** 1.48
      RETURN
    7 IF(N .NE. 7) GO TO 8
      CTHR = 1.1E-12 / (1.+.1*T4)
      IF(J .EQ. 3) CTHR = 5.2E-10
      IF(J .EQ. 4) CTHR = 2.7E-9 * T4 ** .926
      IF(J .EQ. 5) CTHR = 1.7E-10 * T4 ** 1.40
C     ALBERTS N VI
      IF(J .EQ. 6) CTHR = 6.6E-10
      RETURN
    8 IF(N .NE. 8) GO TO 10
      CTHR = 4.0E-10
      IF(J .EQ. 3) CTHR = 7.7E-10 * SQRT(T4)
      IF(J .EQ. 4) CTHR = 2.1E-9
      IF(J .EQ. 5) CTHR = 1.4E-10 * (1.+T4)
C     ALBERTS O VII : NEED O VI
      IF(J .EQ. 7) CTHR = 5.4E-8
      RETURN
   10 IF(N .NE. 10) GO TO 12
      IF(J .EQ. 4) CTHR = 3.8E-9 * SQRT(T4)
      RETURN
   12 IF(N .NE. 12) GO TO 14
      IF(J .EQ. 4) CTHR = 4.4E-9 * (1. + .37*T4)
      RETURN
   14 IF(N .NE. 14) GO TO 16
      IF(J .EQ. 3) CTHR = 4.0E-9 * T4**.23
      IF(J .EQ. 4) CTHR = 4.1E-10
      IF(J .EQ. 5) CTHR = 2.2E-9 * (1.+.1*T4)
      RETURN
   16 IF(N .NE. 16) GO TO 18
      IF(J .EQ. 4) CTHR = 2.5E-9
      IF(J .EQ. 5) CTHR = 7.0E-9
      RETURN
   18 IF(N .NE. 18) GO TO 20
      IF(J .EQ. 4) CTHR = 4.4E-8 * T4 ** .27
      RETURN
   20 RETURN
      END
C
C
C     NAME=DELTT
      FUNCTION DELTT(N,J,L)
C     OUR VERY OWN MODIFICATIONS TO JACOBS
C     DENSITY DEPENDENCE EXTERNAL
C     TAKE 0.5 FOR MERTS ET AL   DELTTA N .NE. 0 EXCEPT S
C     ASSUME 3D OF N,O,F INTERP BETWEEN C AND NE ISOSEQUENCES
C     THIS USES YOUNGERS CLAIM THAT DELTT=1 FOR HE-LIKE RESONANCE LINE;
C     TRY BELY-DUBAU??????
      DELTT = 0.
      IF(J .EQ. 1) RETURN
      I = N-J+1
      IF(L .EQ. 12 .AND. I .NE. 10) GOTO 20
      IF(L .NE. 1 .OR. I .LE. 2) GOTO 50
      IF(I .EQ. 13) GOTO 50
   20 DELTT = 1.
      GOTO 200
   50 CONTINUE
      IF(I .GE. 19) GOTO 19
      IF(I .EQ. 18) GOTO 18
      IF(I .GE. 15) GOTO 14
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14) , I
    1 CONTINUE
C     BURGESS & TWORKOWSKI ARE IN ALPHADI
      IF(L .EQ. 1) DELTT = 1.
      IF(L .GT. 1) DELTT = 0.1
      IF(L .EQ. 5) DELTT = 0.
      GOTO 200
    2 IF (L .EQ. 4 .OR. L .EQ. 5) DELTT = 0.1
      IF (L .EQ. 9) DELTT = 1.0
      GOTO 200
    3 IF (L .GE. 2 .AND. L .LE. 4) DELTT = .5 * (-.35+1.26*N/(N+11.))
      GOTO 200
    4 IF (L .EQ. 7) DELTT = MAX(0.5 * (-.55+1.50*N/(N+11.)), 0.05)
        GOTO 200
    5   IF (L .LE. 6) DELTT = 1.
        IF (L .EQ. 4 .OR. L .EQ. 5) DELTT = 0.5 * (-.54 + 2.2*N/(N+11.))
        GOTO 200
    6   IF (L .LE. 5) DELTT = 1.
        IF (L .EQ. 4) DELTT = MAX(0.5 * (-1.28+3.2*N/(N+11)), 0.05)
        GOTO 200
    7   IF (L .EQ. 4) DELTT = 1.
        IF (L .EQ. 5) DELTT = MAX(0.5*(-1.44+3.4*N/(N+11.)),0.05)
        IF (L .EQ. 7) DELTT = MAX(0.25*(-1.44+3.4*N/(N+11.)), 0.025)
        GOTO 200
    8   IF (L .LE. 3) DELTT = 0.5*(-1.60 + 3.6*N/(N+11.))
        IF (L .EQ. 4 .OR. L .EQ. 5) DELTT = 1.
        IF (L .EQ. 7) DELTT = .25 * (-1.6+3.6*N/(N+11.))
        IF (L .EQ. 15) DELTT = 1.
        GOTO 200
    9   IF (L .LE. 3) DELTT = 0.5 * (-1.75 + 3.8*N/(N+11.))
        IF (L .EQ. 5) DELTT = 1.
        IF (L .EQ. 6) DELTT = 0.25* (-1.75 + 3.8*N/(N+11.))
        GOTO 200
   10   IF (L .LE. 2) DELTT = 1.
        IF (L .GE. 3 .AND. L .LE. 5) DELTT = 0.5 * (-1.9+4.0*N/(N+11.))
        IF (L .GE. 6 .AND. L .LT. 8) DELTT = .25*(-1.9+4.*N/(N+11.))
        IF (L .EQ. 14) DELTT = 1.
        GOTO 200
   11   CONTINUE
        IF (L .EQ. 2 .OR. L .EQ. 3 .OR. L .EQ. 6) DELTT = .25
        GOTO 200
   12   IF (L .EQ. 2) DELTT = .25
        GOTO 200
   13   IF (L .EQ. 2 .OR. L .GE. 9) DELTT = 1.
        IF (L .EQ. 3 .OR. L .EQ. 5) DELTT = .5
        IF (L .EQ. 7) DELTT = 0.25
        GOTO 200
   14   IF (L .LE. 2) DELTT = 1.
        IF (L .EQ. 3 .OR. L .EQ. 4) DELTT = .25
        GOTO 200
   18   DELTT = 0.25
        IF (L .LE. 3) DELTT = 1.
        IF (N .NE. 20) RETURN
        DELTT = 0.
        IF (L .EQ. 6) DELTT = 0.05
        GOTO 200
   19   DELTT = 0.25
  200 CONTINUE
        RETURN
        END
C
C
C     NAME=DIMET
      FUNCTION DIMET(N,J,T,B,DD)
      DIMENSION NOJ(30),E(12,3),F(12,3)
      DATA NOJ/0,1,3*0,2,3,4,0,5,0,6,0,7,0,8,0,9,0,10,5*0,11,0,12,2*0/
      DATA E/0,10.55,13.4,16.3,22.1,27.9,40.0,40.5,46.7,52.8,76.3,86.2,
     . 0.,12.28,16.06,19.8,27.4,35.0,42.7,51.0,60.9,73.5,122.,145.,
     .  5*0.,4.46,9.55,14.5,19.5,24.5,40.7,51./
      DATA F/0.,.26,.21,.18,.16,.12,.10,.091,.081,.075,.003,.002,
     .  0.,.16,.17,.14,.084,.076,.071,.067,.063,.059,.049,.046,
     . 5*0.,.61,.564,0.458,.408,.366,.30,.25/
C
      IJ = N-J-2
      IF(IJ .EQ. 9) IJ = 3
      Z  = J-1
      NO = NOJ(N)
      X  = E(NO,IJ) / (J*13.6)
      A  = SQRT(X) / (1.+.105*X+.015*X*X)
      EBAR  = E(NO,IJ) / (1.+.015*Z*Z*Z/(J*J))
      DIMET = 0.003*T**(-1.5)*A*B*F(NO,IJ)
     .       *EXP(-MIN(50.0,11590.*EBAR/T))*DD
      RETURN
      END
C
C
C     NAME=EFFN
      FUNCTION EFFN(I,ZEFF,T)
      T3 = T/(ZEFF*ZEFF*1000.0)
      XX = 0.4342 * LOG(T3)
      IF(T3-1.0)101,101,102
  101 F1 = 0.266
      F2 = 0.13
      F3 = 0.13
      GOTO 105
  102 IF(T3-10.**5)103,104,104
  103 F1 = 0.266 + 0.1068 * XX - 0.074 * SIN(1.2566*XX)
      F2 = 0.130 + 0.1160 * XX - 0.074 * SIN(1.2566*XX)
      F3=0.130-0.0120*XX+0.050*EXP(-MIN(50.0,(XX-2.)*(XX-2.)))
      GOTO 105
  104 F1 = 0.80
      F2 = 0.71
      F3 = 0.07
  105 CONTINUE
      IF(I .EQ. 0) GOTO 1
      IF(I-18) 110,110,18
  110 GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),I
    1 EYE = I
      EFFN = 2. - EYE
      GOTO 1000
    2 EFFN = 8.0
      GOTO 1000
    3 EFFN = 8. - (4.*F1)
      GOTO 1000
    4 EFFN = 8.*(1.-F1)
      GOTO 1000
    5 EFFN = 6.6667*(1.-F1)
      GOTO 1000
    6 EFFN = 5.33333*(1.-F1)
      GOTO 1000
    7 EFFN = 4. * (1.-F1)
      GOTO 1000
    8 EFFN =  2.6667*(1.-F1)
      GOTO 1000
    9 EFFN =  1.33333*(1.-F1)
      GOTO 1000
   10 EFFN = 18.
      GOTO 1000
   11 EFFN =  18.-(9.*F2)
      GOTO 1000
   12 EFFN = 18.0 * (1.-F2)
      GOTO 1000
   13 EFFN = 18.*(1.-F2) -1.*(9.*F3)
      GOTO 1000
   14 EFFN = 18.*(1.-F2) -2.*(9.*F3)
      GOTO 1000
   15 EFFN = 18.*(1.-F2) -3.*(9.*F3)
      GOTO 1000
   16 EFFN = 18.*(1.-F2) - 4.*(9.*F3)
      GOTO 1000
   17 EFFN = 18.*(1.-F2) -45.0*F3
      GOTO 1000
   18 NO = 28 -I
      EFFN = NO * 1.8 * (1.-3.*F3-F2)
      IF(EFFN)106,1000,1000
  106 EYE = I
      EFFN = 60. -EYE
C     GUARD PACKAGE FOR I = 17
C     PROBABLY UNNECESSARY
      IF(EFFN .LE. 0) EFFN = 1.0
 1000 CONTINUE
      RETURN
      END
C
C
C     NAME=EFFNEW
      FUNCTION EFFNEW(I,T,Z,N)
C
      T3 = .001 * T / Z**2
      XX = .4342 * LOG(T3)
      IF(T3-1.) 101,101,102
  101 F1 = .266
      F2 = .13
      F3 = .13
      GOTO 105
  102 IF( (T3 - 1.E5) .GE. 0.) GOTO 104
      F1 = .266 + .1068*XX - .074*SIN(1.2566*XX)
      F2 = .130 + .1160*XX - .074*SIN(1.2566*XX)
      F3 = .130 - .012 * XX + .05*EXP(-MIN(50.0,(XX-2.)*(XX-2.)))
      GOTO 105
  104 F1 = .80
      F2 = .71
      F3 = .07
  105 CONTINUE
      IF(I .EQ. 2 .OR. I .EQ. 3) EFFNEW = 8. * (1.-F1)
      IF(I .EQ. 10 .OR. I .EQ. 11) EFFNEW = 18. * (1.-F2)
      IF(I .GE. 12 .AND. I .LE. 17) EFFNEW = 18. * (1.-3.*F3 - F2)
      IF(I .GE. 18) EFFNEW = (28. - N + Z) * 1.8 * (1.-3.*F3 - F2)
      RETURN
      END
C
C
C     NAME EXINT1
      FUNCTION EXINT1(X,JUMP)
C
C     R. MOORE OCTOBER 1976
C     JUMP=1    EXINT1=E1(X)
C     JUMP=2    EXINT1=EXP(X)*E1(X)
C     JUMP=3    EXINT1=X*EXP(X)*E1(X)
      IF(X.GE.1.) THEN 
        X2=X*X
        X3=X2*X
        X4=X3*X
        EXINT1=(X4+8.5733287401D0*X3+18.059016973D0*X2
     .         +8.6347608925D0*X+0.2677737343D0
     .         )/(X4+9.5733223454D0*X3+25.6329561486D0*X2
     .           +21.0996530827D0*X+3.9584969228D0)
        IF(JUMP.EQ.1) THEN
          EXINT1=EXINT1*EXP(-MIN(50.0,X))/X
        ELSEIF(JUMP.EQ.2) THEN 
          EXINT1=EXINT1/X
        ENDIF
      ELSE
        EXINT1 =((((((((7.122452D-7*X-1.766345D-6)*X
     .                +2.928433D-5)*X-0.0002335379D0)*X
     .              +0.001664156D0)*X-.01041576D0)*X
     .            +0.05555682D0)*X-0.2500001D0)*X
     .          +0.9999999D0)*X-LOG(X)-0.57721566490153D0
        IF(JUMP.EQ.2) THEN
          EXINT1=EXP(X)*EXINT1
        ELSEIF(JUMP.EQ.3) THEN 
          EXINT1=X*EXP(X)*EXINT1
        ENDIF
      ENDIF
      RETURN
      END
C
C
C
C     NAME=GAUNT
      FUNCTION GAUNT(T,E33,N,J,L,DENO,REDUC,PM)
      INCLUDE "MTLPARAM.h"
C
C     MAY 82 VERSION RELYING HEAVILY ON MANN AND ROBB CALCULATIONS,
C     BHATIA FOR DELTA N = 0 AND BHATIA AND MASON DELTA N = 1
C     DEC 17, 1991 : MODIFY TO IMPROVE TRANS 1 OF LI-LIKE IONS
C
      DIMENSION AH(5),BH(5),GHE(24),GBE(8,5)
      DIMENSION GBE1(4),GBE2(4),GBE3(11),GBE4(11),GBE5(11),BESPEC(4,3)
      DIMENSION GB(8,6),GC(8,9),GN(8,4),GOX(8,8),GF(8,9),GNE(5,12)
      DIMENSION GNA(8,2),GNA4(3,3),GMGII(12)
      DIMENSION GMG(10),GMGP(8,2)
      DIMENSION GAL(10),GAL2(10),A(10),B(10),C(10)
      DIMENSION GMSHELL(6,4),GK(10),GCA(10)
      DIMENSION A1(12),B1(12),C1N(12)
      REAL REDUC(4)
C
C
      DATA AH/0.,.08,.08,.08,0./
      DATA BH/.22,.16,.19,.20,0./
      DATA GHE/.03785,.0002038,.04214,-.001901,-.03263,.004462,.28,0.,
     .  .02327,-.0007424,-.06369,.001865,.07711,.000058,0.,0.,
     .  -.03168,.0008439,.182,-.007289,-.05873,.009605,0.,0./
      DATA GBE/.7114,.00558,-.9277,-.00153,.2711,.00523,0.0,0.0,
     .   0.0,0.0,-.044,.00735,.482,-.01855,0.0,0.0,
     .   .44,0.0,-.010,5*0.0,
     .   .077,.00666,6*0.,
     .   -.2943,-.00084,.2574,.01167,.04852,-.00777,.2671,.00555/
      DATA GBE1/-.0264,0.03,0.0,0.20/
      DATA GBE2/.0881,0.10,.117,0.10/
      DATA GBE3/0.0,.064,.069,.048,.035,.0303,.025,.02,.02,.032,.032/
      DATA GBE4/0.0,.342,.39,.192,.103,.0836,.063,.043,.043,0.0,0.0/
      DATA GBE5/0.0,3*0.0,55000.,64000.,86000.,97000.,108000.,
     .  156000.,168000./
        DATA BESPEC/.00069,.122,.0059,.0063,  .0016,.13,.0073,.0058,
     .  .0017,.127,.0087,.0040/
      DATA GB/.2392,.00543,.3723,.0335,.0321,-.0253,.4591,-.00232,
     .  .2317,.00654,-.1479,.03186,.3038,-.02130,.3904,-.00302,
     .  .0811,.01318,-.2523,.04591,.2122,-.02192,.2304,.00639,
     .  -.112,.0063,.050,.011,.0808,-.005,.239,.0023,
     .  .126,0.0,.287,0.0,0.0,0.0,.28,0.0,
     .  .128,-.025,-.46,.133,.138,-.023,-.70,.128/
      DATA GC/.3539,.00039,.2314,.02314,.0358,-.01534,.5449,-.00858,
     .        .2005,.00395,-.3012,.04332,-.09599,.01606,.3727,-.001517,
     .        .3229,-.002883,.0212,.02481,.0783,-.01326,.4671,-.008106,
     .        -.112,.0063,.050,.011,.0808,-.005,.239,.0023,
     .        .128,-.025,-.46,.133,.138,-.023,-.70,.128,
     .        -.269,-.00117,.0318,.0233,-.0102,-.0075,.235,.0172,
     .        .22,.0056,6*0.,
     .        .198,.0069,6*0.,
     .        8*0./
      DATA GN/-.01876,.0011,.0083,.00182,.0135,-.00083,.040,.00040,
     .  .51,-.1,-1.84,.53,.55,-.091,-2.80,.51,
     .  .128,-.025,-.46,.133,.138,-.023,-.70,.128,
     .  -.07504,.0042,.033,.00726,.054,-.0033,.159,.0015/
C
      DATA GOX/.69,0.0,4*0.0,.46,0.0,
     .  -.112,.0063,.050,.011,.0808,-.005,.2391,.0023,
     .  -.112,.0063,.050,.011,.0808,-.005,.2391,.0023,
     .  .128,-.025,-.46,.133,.138,-.023,-.70,.128,
     .  .51,-.1,-1.84,.53,.55,-.091,-2.80,.51,
     .  .22,.0056,6*0.0,
     .  8*0.,
     .  .1261,0.0,.2869,0.0,0.0,0.0,.28,0.0/
      DATA GF/.2246,.0075,.3055,.0062,-.0575,-.0017,.4248,-.0049,
     .  -.063,0.,.407,0.,-.0238,0.,.478,0.,
     .  -.0157,0.,.188,0.,.0195,0.,.283,0.,
     .  -.001,0.,.134,0.,.123,0.,.158,0.,
     .  .51,-.1,-1.84,.53,.55,-.091,-2.8,.51,
     .  .039,.0154,6*0.,
     .  8*0.,
     .  .22,.056,6*0.,
     .  8*0./
        DATA GNE/-.72117,1.2406,11.746,8.2169,-7.7772,
     .  1.0227,.70828,4.5400,4.1450,-3.8657,
     .  -1.04290,.84411,6.7031,3.1927,-3.1693,
     .  -.039582,.26837,.25803,.086929,-.086641,
     .  -.017323,.29514,.29301,.13223,-.13071,
     .  -.071593,.15504,.12598,-.000449,.000523,
     .  .040202,.25113,.14858,.030780,-.030745,
     .  .45761,.38477,.52142,.92153,-.91649,
     .  -.17862,.32249,.28172,.040677,-.040639,
     .  -.062670,.14921,1.5354,1.0586,-1.0219,
     .  -.057871,.030701,.36471,.14784,-.14293,
     .  .093106,-.001108,-.067652,-.021663,.021657/
C
      DATA GNA/.1586,.007375,.1866,.04156,.02403,-.02416,.3100,-.00098,
     .  -.3245,0.0,.5548,0.0,-.1562,0.0,.266,0./
      DATA GNA4/1.153,-.0333,0.0,  .001,.0053,.058,  .67,-.0133,0.0/
      DATA GMGII/.24,37.3,68.3,96.,96.,2.39,5.99,5.99,0.,0.,0.,.142/
      DATA GMG/.9,.2,.25,.23,.23,.2,.2,.2,.14,.0094/
      DATA GMGP/-.1333,.0155,-.6758,.0577,.5057,-.0323,.314,0.0,
     .          -.294,0.0,.5043,0.0,-.1619,0.0,.2606,0.0/
      DATA GAL/0.0139,.01,.06,0.,0.2,3*0.,1.64,.11/
      DATA GAL2/0.0,.28,.28,0.,0.,3*0.,.00,.28/
      DATA A/.022,.625,.660,.240,.061,.256,.0368,.271,.833,.696/
      DATA B/.0035,.360,.432,.019,.364,.354,.343,.794,.404,.387/
      DATA C/.1261,.1261,.1261,.477,.162,.0108,.138,.0719,.1261,.1261/
      DATA GMSHELL/.453,.91,.13,.93,.2,.2,  .7,.98,.13,.93,0.,0.,
     .  .59,.95,.13,.93,0.,0.,  .51,.95,.20,.29,0.,0./
      DATA GK/.35,.35,1.1,.092,.91,.97,1.1,3*0./
      DATA GCA/.35,.21,43.,.14,.12,.43,37.,.20,.11,4./
C
      DATA A1/0.,.4508,.3957,.5584,.5055,.5055,.5537,.5537,.5108,.5108,
     .       1.235,1.235/
      DATA B1/0.,-.02787,.4532,.1363,.3461,.3461,.3468,.3468,.3304,
     .        .3304,-.6625,-.6625/
      DATA C1N/0.,.1570,-.1219,.1158,-.01872,-.01872,.00516,.00516,
     . .1007,.1007,.8310,.8310/
C
C     SHOULD PUT IN SENSIBLE 'F'S FOR POTASSIUM ISOSEQUENCE
C
      T4    = T / 10000.
      ST    = SQRT(T)
      GAUNT = 0.2
      Y     = E33 * 11590. / T
      CC    = EXINT1(Y,2)
      IF(J .EQ. 1 .AND. N .GT. 2) GOTO 150
C
      II = N - J + 1
      IF(II .GE. 18) GOTO 50
      GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14,14) , II
C
    1 CONTINUE
C     HYDROGENIC IONS: MEWE; HAYES AND SEATON: 3S,3D FROM THOMAS; OTHERS MEWE
      GAUNT = AH(L) + BH(L)*Y*CC + .276*CC
      IF(L .EQ. 5) GAUNT = .212 -.203*(Y-Y*Y*CC) + .0256*CC
      IF(N .GT. 2) RETURN
      IF(L .EQ. 5) GAUNT = .1198*Y*CC - .1194*(Y-Y*Y*CC) + .0328*CC
      GOTO 200
    2 CONTINUE
C     HELIUM-LIKE IONS:  PRADHAN+CASCADES FOR N=2, MEWE N=3,4
C     NEUTRAL HE FROM BERRINGTON ET AL J PHYS B 18,4135; AGGARWAL AP J 278,874
C     L = 9 HAS GAUNT = 0.0,  SATELLITE LINES : L = 6 DONE BY BRANCH FROM L=3
      IF(N .EQ. 2) GOTO 25
      IF(L .GT. 3) GOTO 21
      GAUNT = GNCRCH(GHE,L,Y,CC,N)
      RETURN
   21 IF(L .EQ. 4) GAUNT=(1.+7./N)*(.053+.022*(Y*Y*Y*CC-Y*Y+Y)+.276*CC)
      IF(L .EQ. 5) GAUNT=(1+1.5/N)*(.053+.022*(Y*Y*Y*CC-Y*Y+Y)+.276*CC)
      IF(L .EQ. 6) GAUNT = 0.
      IF(L .EQ. 7) GAUNT = .04 * (Y-Y*Y*CC)
      IF(L .EQ. 8) GAUNT = .02
      IF(L .EQ. 9) GAUNT = 0.
      RETURN
   25 IF(L .EQ. 1) GAUNT = MIN(.00656*T4**(.89),.30*CC)
      IF(L .EQ. 4 .OR. L .EQ. 5) GAUNT = .30 * CC
      IF(L .EQ. 2) GAUNT = MIN(.0283*T4**(-.10),.359*T4**(-.569))
      IF(L .EQ. 3) GAUNT = MIN(.0105*T4**(.30),.095*T4**(-.479))
      IF(L .EQ. 6) GAUNT = 0.0
      IF(L .EQ. 7) GAUNT = 4.08E-7 * T ** .941
      IF(L .EQ. 8) GAUNT = 2.46E-8 * T ** 1.22
      IF(L .EQ. 9) GAUNT = 0.0
      GOTO 200
    3 CONTINUE
C     LITHIUM-LIKE IONS: MEWE MODIFIED TO FIT CLOSE-COUPLING 
C     CALCULATIONS OF ROBB, HENRY L=1 MODIFIED DEC 17, 1991
      I = N-4
      IF(N .GE. 10) I = N/2
      IF(N .GE. 26) I = N/2-2
      IF(L .EQ. 1)GAUNT=A1(I)+(B1(I)*Y-C1N(I)*Y*Y+0.279)*CC+C1N(I)*Y
      X = 1. / (N-3.)
C     IF(L .EQ. 1) GAUNT = (0.68 + .02*N) * ((.7+.35*X) +
C    .  ((1.-.8*X)*Y+(.5-.5*X)*Y*Y +.28)*CC - (.5-.5*X)*Y)
      IF(L .EQ. 2) GAUNT = .053 + .16 * CC
      IF(L .EQ. 3 .OR. L .EQ. 4) GAUNT = -(.16+.32*X)+
     .  ((.8-.56*X)*Y+.2*Y*Y+.28)*CC-.2*Y
      IF(L .EQ. 5) GAUNT = (.19+.25*X) + .079*CC
      IF(L .EQ. 6) GAUNT = .31 - .1*Y*CC
      IF(L .EQ. 7) GAUNT = .096 + .32 * X
      IF (L .EQ. 8) GAUNT = .13
      GOTO 200
    4 CONTINUE
C     BERYLLIUM-LIKE IONS:  QUB FOR N=2 UP THROUGH NE; MANN, ROBB & 
C     SAMPSON PLUS RESONANCES FOR OTHERS MANN AND MALINOVSKY FOR N = 3
C     AND QUB O V N=3 FROM WIDING N = 4 FROM JOHNSTON & KUNZE, 
C     MASON & STOREY (FE XXII) AND LI-LIKE
      I = N-5
      IF(N .GE. 10) I = N / 2 - 1
      IF(N .GE. 26) I = N / 2 - 3
      IF(L .NE. 1) GOTO 41
      TL = LOG10(T)
      GAUNT = .54 + .0125 * N + .135 * CC
      IF(N .LE. 10) GAUNT = GBE1(I) + GBE2(I) * TL
      GOTO 200
   41 IF(L .NE. 2) GOTO 42
      GNT = GBE3(I) / (1. + GBE4(I)/Y)  +  GBE5(I) / T
      IF(N .EQ. 6) GNT = .046 / (1. + T*T/6.25E10)
      GAUNT = GNT * REDUC(1)
      GOTO 200
   42 CONTINUE
      IF(L .LE. 7) GAUNT = GNCRCH(GBE,L-2,Y,CC,N)
      IF(L .EQ. 8) GAUNT = .1261 + .2869*Y*CC + .276*CC
      IF(N .GT. 8) GOTO 200
      EM = -4.67 + 1.86 * N
      EX = EXP(-MIN(50.0,EM*11590./T))
      M = N-5
      IF(L .EQ. 10) GAUNT=BESPEC(1,M)*(1.-PM)*EX+BESPEC(2,M)*PM
      IF(L .EQ. 11) GAUNT=BESPEC(3,M)*(1.-PM)*EX+BESPEC(4,M)*PM
      GOTO 200
    5 CONTINUE
C     BORON ISOSEQUENCE: N=2 INTERPOLATED FROM O IV AND FE XXII OF ROBB; 
C     GOOD TO 10% WRT FLOWER & NUSSBAUMER NA-S, DERE ET AL AR, MANN C II; 
C     OVEREST BY 36% NEAR THRESHOLD WRT ROBB CC FOR C II
C     N = 3 INTERPOLATED C II AND FE XXII FROM MANN; 2S-3L SCALED 
C     FROM BE-, F- & NE-LIKE IONS, AGREES W MASON & STOREY
C     TO 8% FOR 2S-3P.  N = 4 FROM MASON & STOREY
C     L = 9 INTERCOMBINATION LINE
      IF(L .LE. 6) GAUNT = GNCRCH(GB,L,Y,CC,N)
      IF(L .EQ. 9) GAUNT = REDUC(2)
      IF(L .EQ. 12) GAUNT = 0.
      GOTO 200
    6 CONTINUE
C     CARBON ISOSEQUENCE: N=2 INTERP FROM MANN AND ROBB;  
C     AGREES TO 5-10% WITH BHATIA AR, MG
C     N=3,4 GENERAL B-F; AGREES 3-10% WITH MANN 2P-3D; L=6 INTERCOMB
      IF(L .LE. 5) GAUNT = GNCRCH(GC,L,Y,CC,N)
      IF(L .EQ. 7) GAUNT = .198 + .0069*N
      IF(L .EQ. 8) GAUNT = .22 + .0056 * N
      IF(L .EQ. 9) GAUNT = .1261 + (.2869*Y + .28) * CC
      IF(L .EQ. 6) GAUNT = REDUC(3)
      IF(L .EQ. 6 .AND. N .EQ. 6) GAUNT=SQRT(.0001*T)*REDUC(3)
      IF(L .EQ. 12) GAUNT = 0.
      GOTO 200
    7 CONTINUE
C     NITROGEN SEQUENCE:   N=2 FROM MANN FE XX AND MASON &
C     BHATIA MG,SI,S,AR,CA N=3 AND 4 GENERAL B-F
      IF(L .EQ. 1) GAUNT = (1.086-1.71/J)*(.3642+.9358*Y*CC-.3758*
     .  (Y-Y*Y*CC) + .3586 * CC)
      IF(L .GE. 2 .AND. L .LE. 5) GAUNT = GNCRCH(GN,L-1,Y,CC,N)
      IF(L .EQ. 12) GAUNT = 0.
      IF(L .EQ. 8) GAUNT = .22 + .0056*N
      IF(L .EQ. 9) GAUNT = .1261 + .2869*Y*CC + .276*CC
      GOTO 200
    8 CONTINUE
C     OXYGEN ISOSEQUENCE:  N=2 FROM MANN, ROBB FE XIX AND FROM 
C     BHATIA,F&D SI, S, AR OTHERS GENERIC B-F
      IF(L .LE. 8) GAUNT = GNCRCH(GOX,L,Y,CC,N)
      IF(L .EQ. 10) GAUNT = .16 + .0015*N
      IF(L .EQ. 15) GAUNT = 0.
      GOTO 200
    9 CONTINUE
C     FLUORINE SEQUENCE:   N=2 FROM MANN AL V AND ROBB FE XVIII;
C     OTHERS GENERIC B-F ISO; N=3 FROM MANN
      IF(L .LE. 8) GAUNT = GNCRCH(GF,L,Y,CC,N)
      IF(L .EQ. 9) GAUNT = .1261 + (.2869*Y + .28) * CC
      GOTO 200
   10 CONTINUE
C     NEON SEQUENCE: SMITH ET AL INCLUDING CASCADES AND RESONANCES
C     ASSUME THEY ARE ALL LIKE FE XVII
      IF(L .LE. 12) GAUNT = GNE(1,L)+(GNE(2,L)+GNE(3,L)*Y+
     .                      GNE(4,L)*Y*Y)*EXINT1(Y,2)+GNE(5,L)*Y
      GOTO 200
   11 CONTINUE
C     SODIUM SEQUENCE:  MANN, FLOWER & NUSSBAUMER, AND BLAHA
      IF(N .NE. 12) GOTO 111
      GAUNT = GMGII(L)
      IF(L .EQ. 1) GAUNT=.112+(.0269*Y-.0998*Y*Y+.318)*CC+.0998*Y
      IF(L .EQ. 2) GAUNT=141+(59.3*Y-671*Y*Y+.858)*CC+671*Y
      RETURN
  111 CONTINUE
      IF(L .GT. 2) GOTO 112
      GAUNT = GNCRCH(GNA,L,Y,CC,N-10)
      IF(N .EQ. 14 .AND. L .EQ. 2) GAUNT = -0.0172 + (.832*Y
     .                         +.029*Y*Y+.3513)*CC - .029*Y
      RETURN
  112 CONTINUE
      IF(L .GT. 5) GOTO 113
      LLR = L-2
      GAUNT = GNA4(1,LLR) + GNA4(2,LLR)*N + GNA4(3,LLR)*CC
      RETURN
  113 IF(L .EQ. 6) GAUNT = -.16 + .8*Y*CC - .2*(Y-Y*Y*CC)+.276*CC
      IF(L .EQ. 7) GAUNT = .44 - 0.1*Y*CC
      IF(L .EQ. 8) GAUNT = .15 - .05*Y*CC
      IF(L .EQ. 9 .OR. L .EQ. 10) GAUNT = 0.03
      IF(L .EQ. 11) GAUNT = 0.07
      IF(L .EQ. 12) GAUNT = 0.15
      RETURN
   12 CONTINUE
C     MAGNESIUM SEQUENCE:   INTERPOLATE SI TO FE FOR 3P, 4P EXCITATIONS;
C     USE MANN FE GAUNT FACTORS FOR 3D, 4S, 4F, AND INTERCOMB.  ASSUME
C     G FOR 4D = G FOR 4F;  L=10 IS INTERCOMBINATION LINE
      GAUNT = GMG(L)
      IF(L .LE. 2) GAUNT = GNCRCH(GMGP,L,Y,CC,N)
      IF(L .EQ. 10) GAUNT=REDUC(4)
      IF(L .EQ. 10 .AND. N .EQ. 14) GAUNT=REDUC(4)*(20000./T)**(.406)
      RETURN
   13 CONTINUE
C     ALUMINUM ISO-SEQUENCE:   SI II FROM ROBERTS, S IV 3P FROM MANN,
C     3D FROM BHATIA USED FOR AR & CA
C     FE XIV FROM BLAHA USED FOR NI AND FOR ALL N=4
C     RESONANCES INCLUDED FOR METASTABLE USING BHATIA COLLISION STRENGTHS
C     SI II 1814 LINES FROM BROWN, FERRAZ AND JORDAN
C     FE XIV RESONANCES COMPUTED AS IN SMITH ET AL FROM BLAHA OMEGAS
      IF(L .GE. 4 .AND. L .LE. 8) GOTO 132
      IF(N .GT. 14) GOTO 131
      GAUNT = GAL(L) + GAL2(L) * CC
      IF(L .EQ. 1) GAUNT=GAL(1)*1.9E7*ST /(DENO+ST*1.9E7)
      RETURN
  131 IF(N .GT. 20) GOTO 132
      CN = ST * 2.7 * 10. ** (7.+J/2.)
      IF(L .EQ. 1) GAUNT =( (.128 + (.3449*Y+.3544*Y*Y)*CC-.3544*Y) *
     .   1.8 / (1.+.25/Y) )  * .069 *  CN/(DENO+CN)
      IF(L .EQ. 2) GAUNT = .0405+(.2052*Y+.0328*Y*Y+.2311)*CC-.0328*Y
      IF(L .EQ. 3) GAUNT = .142 + .147*CC
      IF(L .EQ. 9) GAUNT = .1330+(.3833*Y+.0934*Y*Y+.1611)*CC-.0934*Y
      IF(L .EQ. 10) GAUNT =.1684+(.2432*Y+.0484*Y*Y+.2638)*CC-.0484*Y
      RETURN
  132 CONTINUE
      GAUNT = A(L) + B(L) * (1.-1. / (1. + C(L)/Y))
      IF(L .EQ. 1) GAUNT=(.0121+.0036/(1.+.1261/Y)) 
     .                  *6.0E17/(DENO+6.0E17)
      RETURN
C
C     SILICON THROUGH CHLORINE SEQUENCES
   14 GAUNT = GMSHELL(L,II-13)
      IF(L .EQ. 3 .AND. II .NE. 17)  GAUNT = GMSHELL(L,II-13) 
     .                                     - 0.114/(1.+.138/Y)
      RETURN
   50 CONTINUE
      IF(II .GT. 18) GOTO 19
      RETURN
   19 IF(II .GT. 19) GOTO 20
      GAUNT = GK(L)
      IF(N .EQ. 20) GAUNT = GCA(L)
      RETURN
   20 IF(L .LE. 3) GAUNT = MAX(.2,-.606*LOG10(Y)-.052)
      RETURN
  150 CONTINUE
      GAUNT = .01
      IF(Y .LE. 10.) GAUNT = MIN(EXP(-MIN(50.0,.7*LOG(Y)-2.3)),
     .                          .28*CC)
  200 CONTINUE
      RETURN
      END
C
C
C     NAME=GNCRCH
      FUNCTION GNCRCH(R,L,Y,CC,N)
      DIMENSION R(1)
      IL = 8 * (L-1)
      A = R(IL+1) + R(IL+2) * N
      B = R(IL+3) + R(IL+4) * N
      C = R(IL+5) + R(IL+6) * N
      D = R(IL+7) + R(IL+8) * N
      GNCRCH = A + B*Y*CC + C*(Y-Y*Y*CC) + D*CC
      RETURN
      END
C
C
C     NAME=GND
      FUNCTION GND(NUM)
      GND = 4.
      IF(NUM .LT. 28) GND = 3.
      IF(NUM .LT. 10) GND = 2.
      IF(NUM .LT. 2 ) GND = 1.
      RETURN
      END
C
C
C     NAME=GREC
      FUNCTION GREC(NO,N,J,EE,T)
C
      INCLUDE "MTLPARAM.h"
C
      DIMENSION DRAT(30)
      DATA DRAT/2.,.5,2.,.5,6.,2.5,2.22,2.25,.67,.167,2.,.5,6.,2.5,
     .  2.22,2.25,.67,.167,12*0./
      DEG = DRAT(N-J + 1)
      X1  = EE*11590./T
      X12 = X1*X1
      X14 = X12*X12
      G = 5.238E-23*T**1.5*DEG*(S2J(NO,J)*X12+S4J(NO,J)*X12*X1+
     .  X12*X1*S5J(NO,J)*(0.5-X1/2.)+EXINT1(X1,2)
     . *(S3J(NO,J)*X12*X1-X14*S4J(NO,J)+X1*X14*S5J(NO,J)/2.))
      GREC = G
      IF(X1 .GE. 300.) GREC = 5.238E-23*T**1.5*DEG*(S2J(NO,J)
     .  + S3J(NO,J) + S4J(NO,J) + S5J(NO,J)) * X12
      RETURN
      END
C
C
C
C     NAME=NOSON
      FUNCTION NOSON (A,X,L,LMAX)
C     RUDOLF LOESER, 28 JULY 1965
C     NOSON IS A MATRIX INVERSION ROUTINE, USING THE METHOD OUTLINED ON
C     P 434 OF HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,
C     (NEW YORK,1956).
C     A IS A MATRIX OF ORDER L, WITH COLUMN DIMENSION LMAX, ITS ELEMENTS
C     ARE ASSUMED TO BE STORED COLUMNWISE, IN THE USUAL FORTRAN MANNER,
C     X IS WORKING STORAGE OF LENGTH L.
C     THE INVERSE OF A WILL REPLACE A.
C     UPON RETURN, NOSON=1 IF INVERSION WENT OK, =0 IF A DIVISOR IS
C     ZERO (IN WHICH CASE A MAY CONTAIN GARBAGE UPON RETURN).
      DIMENSION A(1),X(1)
      N = L - 1
      MAX = N * LMAX + L
      DO I = 1,L
         X(I) = 1.
      ENDDO
      K1 = -LMAX
      DO K=1,L
         K1 = K1 + LMAX
         K2 = K1 + K
         IF (A(K2) )101,115,101
  101    CONTINUE
         DO I = 1,L
            J1 = K1 + I
            IF ( A(J1) )102,104,102
  102       F = 1. / A(J1)
            X(I) = X(I) * F
            DO J1 = I,MAX,LMAX
               A(J1) = A(J1) * F
            ENDDO
  104       CONTINUE
         ENDDO
         A(K2) = X(K)
         X(K) = 1.
         DO I=1,L
            KI = K-I
            IF( KI )105,108,105
  105       CONTINUE
            J1 = K1 + I
            IF( A(J1) )106,108,106
  106       CONTINUE
            A(J1) = 0.
            DO J2 = I, MAX, LMAX
               J1 = J2 + KI
               A(J2) = A(J2) - A(J1)
            ENDDO
  108       CONTINUE
         ENDDO
      ENDDO
      DO I = 1,N
         IF ( X(I) )111,115,111
  111    CONTINUE
         F = 1. / X(I)
         DO J1 = I, MAX, LMAX
            A(J1) = A(J1) * F
         ENDDO
      ENDDO
      NOSON = 1
      RETURN
  115 NOSON = 0
      RETURN
      END
C
C
C
C     NAME=POPMET
      FUNCTION POPMET(N,IJ,EMET,DENO,T,REDUC)
C     RETURNS POPULATION OF METASTABLE LEVEL OF BE,B,C, OR MG-LIKE ION;  
C     CARBON-NICKEL RESONANCE CONTRIBUTIONS CALC FOR HIGH Z WITH 
C     SMITH ET AL METHOD PASSES REDUC BACK TO GAUNT FOR INTERCOMBINATION LINES
      INCLUDE "MTLPARAM.h"
C
      REAL MGMM,MGM1,MGM5,MG21,MG5,MG1
      DIMENSION RM5(12),EBE(12),POPS(4)
      DIMENSION MGMM(3,7),MGM1(7),MG5(7),MG21(2,7),EMG(7)
      DIMENSION NK(30),A21(2,12),B21(3,12),C21(3,12),RM1(12)
      DIMENSION RMM(3,12),BM1(12),BMM(3,12),CM1(3,12),CMM(2,12)
      DIMENSION ECM(2,12),R(3,3),X(3),C(3),CPR(3)
      DIMENSION RB5(10),RC5(10),RRES(4,10),ERES(4,10)
      DIMENSION RBE(10),BEE(10)
      REAL REDUC(4)
C
      DATA RM5/0.,18.6,11.7,6.09,2.79,1.54,.93,.61,.41,.31,.13,.10/
      DATA EBE/0.,6.20,7.87,9.51,12.81,16.1,19.5,23.11,26.9,28.9,55.5,
     .  63.1/
      DATA MG21/430.,5.8E-4,22000.,.008,1.63E+5,.000,9.E+5,.12,3.3E+6,
     .  .26,4.5E+7,1.00,8.E+7,2./
      DATA MG5/0.2,0.5,.26,.035,.026,.0058,.0039/
C
C     GILES QUOTED BY MENDOZA FOR S V, BALUJA FOR SI III
C
      DATA MGM1/.045,.35,.072,.026,.016,.0069,.0047/
      DATA MGMM/.5,.38,1.5,1.4,1.06,4.1,.25,.39,1.32,
     .  .29,.21,.82,.17,.13,.51,.033,.25,.43,.02,.13,.30/
      DATA EMG/1.63,3.72,5.50,7.13,8.72,13.9,17./
      DATA NK/0,1,0,0,0,2,3,4,0,5,0,6,0,7,0,8,0,9,0,10,5*0,11,0,12,0,0/
      DATA A21/2*0.,77.2,.00282,471.,.00629,1880.,.0118,16200.,.0306,
     .  84500.,.0629,320000.,.112,972600.,.1823,2.546E+6,.2781,
     .  5.914E+6,.4029,7.3E+7,1.21,1.4E+8,1.8/
      DATA B21/3*0.,332.,130.,82.,  1183.,130.,347.,4210.,499.,1470.,
     . 28400.,3170.,8910.,152000.,19800.,56200.,705000.,102000.,288000.,
     . 43790.,11613.,27300.,4.28E+06,7.08E+05,2.42E+06,
     .  1.02E+07,1.55E+06,5.62E+06,8.E+7,1.E+7,5.E+7,2.E+8,2.E+7,1.E+8/
      DATA C21/3*0.,.00031,.0026,32.,.0041,.0342,190.,.022,.16,385.,
     .  .599,4.64,6910.,4.67,37.1,45200.,26.1,200.,197000.,115.,856.,
     . 680000.,424.,3000.,2.0E+06,1370.,9090.,5.32E+06,
     .  29400.,128000.,6.7E+7,95000.,380000.,1.7E+8/
C
C     REDUCED COLLISION RATES * 10. ** 7
C
      DATA RM1/0.,9.384,4.944,2.814,1.340,.7265,.4718,.3164,.2243,
     .  .1732,.08,.062/
C
C     INCLUDES PROTON RATES FROM MUNRO THROUGH SI, EXT F OR S, NONE ABOVE A
C
      DATA RMM/3*0.,29.9,13.4,50.8,17.7,8.41,32.2,11.2,6.00,21.8,
     .  6.48,4.19,13.9,3.96,3.04,9.68,2.91,2.67,7.96,2.19,2.42,6.83,
     .  1.25,0.40,1.79,.533,.317,1.40,.046,.16,.67,.020,.12,.52/
      DATA BM1/0.,14.,8.1,5.3,4.1,3.01,2.0,1.5,1.1,.95,.67,.58/
      DATA BMM/3*0.,23.,7.9,31.,17.1,5.73,20.8,12.3,4.17,14.1,
     .  6.52,2.39,7.91,3.84,1.47,3.01,2.65,1.01,3.16,18.6,13.9,38.,
     .  1.68,.57,1.9,1.34,.43,1.45,.69,.18,.69,.55,.14,.54/
      DATA CM1/3*0.,12.1,12.9,10.,51.,29.5,23.,12.0,8.6,4.18,
     . 9.1,6.4,2.92,6.75,4.72,2.1,5.39,3.86,1.55,4.17,3.04,1.19,
     .  2.84,2.02,1.00,2.43,1.57,0.90,1.1,.75,.60,1.0,.56,.54/
      DATA CMM/2*0.,12.9,.0005,32.5,.001,29.,.0015,17.5,.00205,11.2,
     .  .00173,7.42,.00173,5.44,.00518,4.1,.0069,3.52,.0138,1.7,.03,1.5,
     .  .05/
      DATA ECM/2*0.,1.26,2.68,1.90,4.05,2.42,5.34,3.76,7.93,5.09,10.6,
     . 6.58,13.4,8.34,16.5,10.6,20.1,13.5,24.5,30.1,45.,51.,62./
      DATA RB5/0.,25.,21.6,16.9,8.80,5.23,3.64,15.5,2*0./
      DATA RC5/10*0./
      DATA BEE/0.,10.5,13.4,16.3,22.1,28.,34.,40.5,2*0./
      DATA RBE/0.,115.,95.4,78.4,47.7,32.1,24.6,18.2,2*0./
      DATA ERES/4*0.,12.7,9.25,7.94,0.,16.3,12.5,11.5,0.,
     .  19.7,15.8,14.9,0.,27.0,22.3,21.9,0.,33.8,28.8,28.9,4.35,
     .  40.9,35.8,36.3,10.3,48.4,11.6,44.1,15.8,  8*0./
      DATA RRES/4*0.,327.,95.,5.9,0.,306.,111.,74.9,0.,
     .  229.,97.8,139.,0.,141.,60.7,111.,0.,96.,41.3,70.3,68.0,
     .  71.7,30.3,37.1,956.,52.5,22.0,26.5,64.3,8*0./
C
      DO I = 1,3
         DO J=1,3
            R(I,J) = 0.
         ENDDO
      ENDDO
C
C     AVOID OVERFLOW FOR LOW T
      REDUC(IJ) = 0.
      POPMET = 0.
      IF(EMET*11590./T .GE. 25.) RETURN
      ST = SQRT(T)
      NO = NK(N)
      V  = MAX(DENO * 1.E-7 / ST, 1.E-7/ST)
      GOTO (1,2,3,4) , IJ
    1 CONTINUE
      R(1,1) = -(RM1(NO)+RMM(1,NO)*3.+RMM(2,NO)*5.) * V
      R(2,2) = -A21(1,NO) - V*(RM1(NO)+RMM(1,NO)+RMM(3,NO)*1.6666667)
      R(3,3) = -A21(2,NO) - V * (RM1(NO)+RMM(2,NO)+RMM(3,NO))
      R5 = RM5(NO)*V*EXP(-MIN(50.0,11590.*EBE(NO)/T))
      R(1,1) = R(1,1) - R5
      R(2,2) = R(2,2) - R5
      R(3,3) = R(3,3) - R5
      R(2,1) = V * RMM(1,NO) * 3.
      R(3,1) = V * RMM(2,NO) * 5.
      R(1,2) = V * RMM(1,NO)
      R(1,3) = V * RMM(2,NO)
      R(3,2) = V * 1.666667 * RMM(3,NO)
      R(2,3) = V * RMM(3,NO)
      C(1) =-RM1(NO) * EXP(-MIN(50.0,EMET*11590./T)) * V
      C(2) = C(1) * 3.
      C(3) = C(1) * 5.
      GOTO 10
    2 CONTINUE
C     B-LIKE IONS
      R5 = RB5(NO)*V*EXP(-MIN(50.0,11590.*(ERES(2,NO)-EMET)/T))
      R(1,1) = -(BM1(NO)+BMM(1,NO)*2.+BMM(2,NO)*3.)*V - B21(1,NO)
      R(2,2) = -(BM1(NO)+BMM(1,NO)+BMM(3,NO)*1.5)*V - B21(2,NO)
      R(3,3) = -(BM1(NO)+BMM(2,NO)+BMM(3,NO))*V - B21(3,NO)
      R(2,1) = V * BMM(1,NO) * 2.
      R(3,1) = V * BMM(2,NO) * 3.
      R(1,2) = V * BMM(1,NO)
      R(1,3) = V * BMM(2,NO)
      R(3,2) = V * 1.5 * BMM(3,NO)
      R(2,3) = V * BMM(3,NO)
      C(1) = -BM1(NO) * EXP(-MIN(50.0,EMET*11590./T))*V / 3.
      C(2) = C(1) * 2.
      C(3) = C(1) * 3.
      GOTO 10
    3 CONTINUE
C     C-LIKE IONS
      R5 = RC5(NO) * V * EXP(-MIN(50.0,11590.*(ERES(3,NO)-EMET)/T))
      E1 = EXP(-MIN(50.0,ECM(1,NO)*11590./T))
      E2 = EXP(-MIN(50.0,ECM(2,NO)*11590./T))
      E3R = EXP(-MIN(50.0,EMET*11590./T))
      R(1,1)=-(CM1(1,NO)+CMM(1,NO)*.2*E2/E1
     .      +CMM(2,NO)*E3R/E1)*V-C21(1,NO)
      R(2,2) = -(CM1(2,NO)+CMM(1,NO))*V- C21(2,NO)
      R(3,3) = -(CM1(3,NO) + CMM(2,NO)) * V - C21(3,NO)
      R(2,1) = V * CMM(1,NO) * .2 * E2/E1
      R(1,2) = V * CMM(1,NO)
      R(3,1) = V * CMM(2,NO) * E3R/E1
      R(1,3) = V * CMM(2,NO)
      R(3,2) = 0.
      R(2,3) = 0.
      C(1) = -CM1(1,NO) * V * 5. * E1 / 9.
      C(2) = -CM1(2,NO) * V * E2 / 9.
      C(3) = -CM1(3,NO) * V * E3R * 5. / 9.
      GOTO 10
    4 CONTINUE
      V = 8.63E-6 * DENO / ST
      NP = NO - 5
      MG1 = V * MGM1(NP)
      MGM5 = V * MG5(NP) * EXP(-MIN(50.0,11590.*EMG(NP)/T))
      IF(N .EQ. 14) MGM5 = MGM5 / (1.+8.6E-6*T)
      IF(N .EQ. 14) MG1 = MG1 / (1.+T*3.5E-6)
      R(1,1) = -MG1 - MGM5 - MGMM(1,NP) * V - MGMM(2,NP) * V
      R(2,2) = -MG1-MGM5-MGMM(1,NP)*V/3-MGMM(3,NP)*V/3-MG21(1,NP)
      R(3,3) = -MG1-MGM5-MG21(2,NP)-MGMM(2,NP)*V/5.-MGMM(3,NP)*V/5.
      R(2,1) = V*MGMM(1,NP)
      R(3,1) = V*MGMM(2,NP)
      R(1,2) = V*MGMM(1,NP) / 3.
      R(1,3) = V*MGMM(2,NP) / 5.
      R(3,2) = V*MGMM(3,NP) / 3.
      R(2,3) = V * MGMM(3,NP) / 5.
      C(1) = -MG1 * EXP(-MIN(50.0,11590.*EMET/T))
      C(2) = C(1) * 3.
      C(3) = C(1) * 5.
   10 CONTINUE
      NN = NOSON(R,X,3,3)
      IF(NN .EQ. 0) PRINT 99
   99 FORMAT (' BAD INVERSION')
c     IF(NN .EQ. 0) THEN
c       WRITE(*,*)'BAD INVERSION T,DENO= ',T,DENO
c     ENDIF
C
      DO I = 1,3
         CPR(I) = 0.
         DO J = 1,3
            CPR(I) = CPR(I) + R(I,J) * C(J)
         ENDDO
      ENDDO
      SUM = 1. + CPR(1) + CPR(2) + CPR(3)
      POPS(1) = 1./SUM
      POPS(2) = CPR(1) / SUM
      POPS(3) = CPR(2) / SUM
      POPS(4) = CPR(3) / SUM
      POPMET = (CPR(1) + CPR(2) + CPR(3)) / SUM
      IF(IJ .EQ. 3) POPMET = CPR(3) / SUM
      GOTO (51,52,53,54) , IJ
C     ERGS * 10-23 / S / ATOM
   51 SUMC = -C(1) - C(2) - C(3)
      REDUC(1) = (POPS(3)*A21(1,NO)+POPS(4)*A21(2,NO))/SUMC
      RETURN
   52 SUMC = -C(1) - C(2) - C(3)
      REDUC(2)=(POPS(2)*B21(1,NO)+POPS(3)*B21(2,NO)+POPS(4)*B21(3,NO))
     .  /SUMC
      RETURN
   53 REDUC(3) = -POPS(4)*C21(3,NO) / C(3)
      RETURN
   54 SUMC = -C(1) - C(2) - C(3)
      REDUC(4) = (POPS(3)*MG21(1,NP)+POPS(4)*MG21(2,NP))/SUMC
      RETURN
      END
C
C
C     NAME=SEATON
      FUNCTION SEATON(X,Q)
C     R. MOORE OCTOBER 1976
      IMPLICIT REAL*8(A,B)
      DOUBLE PRECISION XD,XS1,XS2,X3
      DATA A10/8.7469604697013D-4/    ,B10/-1.2007397521051D-4/
      DATA A11/.20406771231267D0/       ,B11/-.055223640825293D0/
      DATA A12/-2.1524482972354D0/      ,B12/.029171138841798D0/
      DATA A13/12.663578339302D0/       ,B13/.25091093604147D0/
      DATA A14/-43.153566859883D0/      ,B14/-.94344532109356D0/
      DATA A15/61.113098339262D0/
      DATA A20/.011834608639468D0/      ,B20/-1.2467753947278D-3/
      DATA A21/-2.9889195903436D-3/   ,B21/-.043871158058636D0/
      DATA A22/-.10946027271945D0/      ,B22/.013617064094285D0/
      DATA A23/.097410292577482D0/      ,B23/-5.2824884665512D-3/
      DATA A24/-.039676727608179D0/     ,B24/8.9490487211065D-4/
      DATA A25/6.2318768197420D-3/
      DATA A30/.018985613015081D0/      ,B30/-.010963734653233D0/
      DATA A31/-.064707516794785D0/     ,B31/-.028119928428050D0/
      DATA A32/6.2172804659938D-3/    ,B32/1.1971209805431D-3/
      DATA A33/-4.1777961107942D-4/   ,B33/-4.8739343472085D-5/
      DATA A34/1.5264422582645D-5/    ,B34/8.4274360230135D-7/
      DATA A35/-2.2708774951499D-7/
C
      IF(X.GE.20.)  GOTO 40
      IF(X.GT.2.)   GOTO 30
      IF(X.GT.0.2)  GOTO 20
      IF(X.GT.0.02) GOTO 10
      XS1=0.4629*X*(1.+4.*X)-1.0368*X**1.3333333333333
     .          *(1.+1.875*X)
      XS2=-0.0672*X*(1.+3.*X)+.1488*X**1.6666666666667
     .           *(1.+1.8*X)
      GOTO 50
   10 XS1=A10+(A11+(A12+(A13+(A14+A15*X)*X)*X)*X)*X
      XS2=B10+(B11+(B12+(B13+B14*X)*X)*X)*X
      GOTO 50
   20 XS1=A20+(A21+(A22+(A23+(A24+A25*X)*X)*X)*X)*X
      XS2=B20+(B21+(B22+(B23+B24*X)*X)*X)*X
      GOTO 50
   30 XS1=A30+(A31+(A32+(A33+(A34+A35*X)*X)*X)*X)*X
      XS2=B30+(B31+(B32+(B33+B34*X)*X)*X)*X
      GOTO 50
   40 X3=3.*X
      XS1=-0.1728*X**.33333333333333*(1.+(-8.
     .         +(70.+(-800.+11440./X3)/X3)/X3)/X3)
      XS2=-0.0496*X**.66666666666667*(1.+(-3.
     .                      +(32.-448./X3)/X3)/X3)
   50 X3=X**0.33333333333333*Q**0.66666666666667
      SEATON=EXINT1(X,3)+(XS1+XS2/X3)/X3
#endif /* CEN_METALS */
      RETURN
      END
