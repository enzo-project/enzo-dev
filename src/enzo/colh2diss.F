#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c/////////////////////  SUBROUTINE COLH2DISS  \\\\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine colh2diss(tgas, f1, f2, f3, f4, f5, f6, f7)
c
c  COMPUTE DENSITY DEPENDENT COLLISIONAL H2 DISSOCIATION RATE
c
c  written by: Tom Abel
c  date:       
c  modified1: Feb, 2000 by Greg Bryan; adapted to AMR
c
c  PURPOSE:
c    Computes the 7 temperature dependant functions required to
c      generate the density-dependant k13 rate.
c
c     compute density dependent collisional H2 dissociation by HI
c     data from Martin, Schwartz, Mandy, 1996, ApJ 461, 265
c     returns log (base 10) of the rate coefficient in cm^3/s
c     of the reaction:
c     H2   +   H   ->  3 H
c       Tom Abel 2000
c
c  UNITS:
c    log10(cgs)
c
c  PARAMETERS:
c
c  INPUTS:
C     T is the gas temperature in Kelvin
c
c  OUTPUTS:
c     f1-7: rates as given below
c
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c  Arguments
c
      R_PREC f1,f2,f3,f4,f5,f6,f7
      real*8 t, tgas
c
c  Parameters
c
c
c  Locals
c
      real*8 tl, a1, a, b, b1, c, c1, d, ftd
      real*8 y(21)
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
c
      t  = tgas
      f1 = tiny
      f2 = tiny
      f3 = tiny
      f4 = tiny
      f5 = 1._RKIND
      f6 = 1._RKIND
      f7 = 0._RKIND
c
c     do not use tiny values for temperatures below 500 K also do not
c     return values for temperatures above 1 million Kelvin. Note that     
c     data and fits are only accurate for t< 5 10^5. However,
c     collisional dissociation by electrons will be dominant above this
c     temperature, anyhow. 
c
      if (tgas .le. 500._RKIND) then
         t=500._RKIND
c         CID = -60._RKIND
      endif
      if (t .ge. 1.e6_RKIND) then
         t = 1.e6_RKIND
c         CID = -60._RKIND
      endif
c
c     fitting parameters
c
      y(1)   =   -1.784239e2_RKIND
      y(2)   =   -6.842243e1_RKIND
      y(3)   =    4.320243e1_RKIND
      y(4)   =   -4.633167e0_RKIND
      y(5)   =    6.970086e1_RKIND
      y(6)   =    4.087038e4_RKIND
      y(7)   =   -2.370570e4_RKIND
      y(8)   =    1.288953e2_RKIND
      y(9)   =   -5.391334e1_RKIND
      y(10)  =    5.315517e0_RKIND
      y(11)  =   -1.973427e1_RKIND
      y(12)  =    1.678095e4_RKIND
      y(13)  =   -2.578611e4_RKIND
      y(14)  =    1.482123e1_RKIND
      y(15)  =   -4.890915e0_RKIND
      y(16)  =    4.749030e-1_RKIND
      y(17)  =   -1.338283e2_RKIND
      y(18)  =   -1.164408e0_RKIND
      y(19)  =    8.227443e-1_RKIND
      y(20)  =    5.864073e-1_RKIND
      y(21)  =   -2.056313e0_RKIND
c
      tl=log10(t)
c high density limit
      a =   y(1)+y(2)*tl+y(3)*tl*tl+y(4)*tl*tl*tl
     $     +y(5)*dlog10(1._RKIND+y(6)/t)
      a1=   y(7)/t
c low density limit
      b =   y(8)+y(9)*tl+y(10)*tl*tl+y(11)*dlog10(1._RKIND+y(12)/t)       
      b1=   y(13)/t
c critical density
      c =   y(14)+y(15)*tl+y(16)*tl*tl+y(17)/t
      c1 =  y(18)+c
      d =   y(19)+y(20)*exp(-t/1850._RKIND)+y(21)*exp(-t/440._RKIND)    
c     tabulate the following temperature dependent coefficients:
      f1 = a
      f2 = (a-b)
      f3 = a1   
      f4 = (a1-b1)
      f5 = 10._RKIND**c
      f6 = 10._RKIND**c1
      f7 = d
c
c and then get the dissociation rate in cm^3/s with this
c    (note: this is in log10)
c
c      CID = f1-f2/(1._RKIND+(nh/f5)**f7)
c     $     +f3-f4/(1._RKIND+(nh/f6)**f7)
c
c
      return
      end
