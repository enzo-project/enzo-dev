/***********************************************************************
/
/  Modified Bessel Functions of various order and kind
/
/  adopted by: Tom Abel 
/  date:       May, 2010
/  modified1:
/
/
************************************************************************/

#include <stdio.h>
#include <math.h>


/* --------------------------------------------------------------------
! Program to calculate the first kind modified Bessel function
! of integer order N, for any REAL X, using the function BESSI(N,X).
! ---------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                         Visual C++ Release 1.0 By J-P Moreau, Paris
! --------------------------------------------------------------------*/

  double BESSI0(double X);

  double BESSI1(double X);

// ---------------------------------------------------------------------
  double BESSI(int N, double X) {
/*
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
*/

      int IACC = 40; 
	  double BIGNO = 1e10, BIGNI = 1e-10;
      double TOX, BIM, BI, BIP, BSI;
      int J, M;

      if (N==0)  return (BESSI0(X));
      if (N==1)  return (BESSI1(X));
      if (X==0.0) return 0.0;

      TOX = 2.0/X;
      BIP = 0.0;
      BI  = 1.0;
      BSI = 0.0;
      M = (int) (2*((N+floor(sqrt(IACC*N)))));
      for (J = M; J>0; J--) {
        BIM = BIP+J*TOX*BI;
        BIP = BI;
        BI  = BIM;
        if (fabs(BI) > BIGNO) {
          BI  = BI*BIGNI;
          BIP = BIP*BIGNI;
          BSI = BSI*BIGNI;
        }
        if (J==N)  BSI = BIP;
      }
      return (BSI*BESSI0(X)/BI);
  }

// ----------------------------------------------------------------------
//  Auxiliary Bessel functions for N=0, N=1
  double BESSI0(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067429;
      P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
      Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
      Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
      Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

// ---------------------------------------------------------------------
  double BESSI1(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=0.5; P2=0.87890594; P3=0.51498869; P4=0.15084934;
      P5=0.2658733e-1; P6=0.301532e-2; P7=0.32411e-3;
      Q1=0.39894228; Q2=-0.3988024e-1; Q3=-0.362018e-2;
      Q4=0.163801e-2; Q5=-0.1031555e-1; Q6=0.2282967e-1;
      Q7=-0.2895312e-1; Q8=0.1787654e-1; Q9=-0.420059e-2;
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return(X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

// ---------------------------------------------------------------------

/**********************************************************
* Reference:   From Numath Library By Tuan Dang Trong     *
* in Fortran 77 [BIBLI 18].                               *
*                                                         *
*                      C++ Version By J-P Moreau, Paris.  *
**********************************************************/

#define BIG  1e30

    double X, Y;
    int i, N;

    double BESSK0(double X);
    double BESSK1(double X);
    double BESSI0(double X);
    double BESSI1(double X);


    double BESSK(int N, double X) {
/* -----------------------------------------------------------------------
!     THIS ROUTINE CALCULATES THE MODIFIED BESSEL FUNCTION OF THE THIRD
!     KIND OF INTEGER ORDER, N FOR ANY POSITIVE REAL ARGUMENT, X. THE
!     CLASSICAL RECURSION FORMULA IS USED, STARTING FROM BESSK0 AND BESSK1.
! ------------------------------------------------------------------------ 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
! ----------------------------------------------------------------------*/
      double TOX,BK,BKM,BKP;
      int J;
      if (N == 0) return BESSK0(X);
      if (N == 1) return BESSK1(X);
      if (X == 0.0) return BIG;  //arbitrary big value
      TOX = 2.0/X;
      BK  = BESSK1(X);
      BKM = BESSK0(X);
      for (J=1; J<N; J++) {
        BKP = BKM + J*TOX*BK;
        BKM = BK;
        BK  = BKP;
      }
      return BK;
    }

    double BESSK0(double X)  {
/* ---------------------------------------------------------------------
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
! --------------------------------------------------------------------*/
      double Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,tmp;
      P1=-0.57721566; P2= 0.42278420; P3=0.23069756; P4=0.3488590e-1;
      P5= 0.262698e-2; P6=0.10750e-3; P7=0.74e-5;
      Q1= 1.25331414; Q2=-0.7832358e-1; Q3=0.2189568e-1; Q4=-0.1062446e-1;
      Q5= 0.587872e-2; Q6=-0.251540e-2; Q7=0.53208e-3;
      if (X == 0.0) return BIG;  //arbitrary big value
      if (X <= 2.0) {
        Y=X*X/4.0;
        AX=-log(X/2.0)*BESSI0(X);
        tmp = AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
	return tmp;
      }
      else {
        Y=2.0/X;
        AX=exp(-X)/sqrt(X);
        tmp = AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))));
	return tmp;
      }
    }

    double BESSK1(double X)  {
/* ------------------------------------------------------------------------
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
! -----------------------------------------------------------------------*/
      double Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7;
      P1=1.0; P2=0.15443144; P3=-0.67278579; P4=-0.18156897;
      P5=-0.1919402e-1; P6=-0.110404e-2; P7=-0.4686e-4;
      Q1=1.25331414; Q2=0.23498619; Q3=-0.3655620e-1; Q4=0.1504268e-1;
      Q5=-0.780353e-2; Q6=0.325614e-2; Q7=-0.68245e-3;
      if (X == 0.0) return BIG;
      if (X <= 2.0) {
        Y=X*X/4.0;
        AX=log(X/2.0)*BESSI1(X);
        return (AX+(1.0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
      }
      else {
        Y=2.0/X;
        AX=exp(-X)/sqrt(X);
        return (AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7)))))));
      }
    }


