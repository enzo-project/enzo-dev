#include "../enzo/fortran.def"
!=======================================================================
! Set R_PREC to match what's in fortran_types.def (since it cannot be 
! included before the first function declaration)
#ifdef CONFIG_BFLOAT_4
#define R_PREC real*4
#endif

#ifdef CONFIG_BFLOAT_8
#define R_PREC real*8
#endif
!=======================================================================

      R_PREC function ran1(idum)

      implicit none
#include "../enzo/fortran_types.def"

!     Parameters

      INTG_PREC, parameter :: m1 = 259200, ia1 = 7141, ic1 = 54773
      INTG_PREC, parameter :: m2 = 134456, ia2 = 8121, ic2 = 28411
      INTG_PREC, parameter :: m3 = 243000, ia3 = 4561, ic3 = 51349

      R_PREC, parameter :: rm1 = 3.8580247e-6_RKIND
      R_PREC, parameter :: rm2 = 7.4373773e-6_RKIND

!     Argument

      INTG_PREC :: idum

!     Locals

      INTG_PREC :: iff, ix1, ix2, ix3, j
      R_PREC :: r(97)

      data iff /0/
      save iff, r, ix1, ix2, ix3

      if ((idum < 0) .or. (iff == 0)) then
        iff = 1
        ix1 = mod(ic1-idum, m1)
        ix1 = mod(ia1*ix1+ic1, m1)
        ix2 = mod(ix1, m2)
        ix1 = mod(ia1*ix1+ic1, m1)
        ix3 = mod(ix1, m3)

        do j = 1, 97
          ix1 = mod(ia1*ix1+ic1, m1)
          ix2 = mod(ia2*ix2+ic2, m2)
          r(j) = (REAL(ix1,RKIND)+REAL(ix2,RKIND)*rm2)*rm1
        end do

        idum = 1
      endif

      ix1 = mod(ia1*ix1+ic1, m1)
      ix2 = mod(ia2*ix2+ic2, m2)
      ix3 = mod(ia3*ix3+ic3, m3)
      j = 1+(97*ix3)/m3

      if((j > 97) .or. (j < 1)) stop 'j>97 or j<1'

      ran1 = r(j)
      r(j) = (REAL(ix1,RKIND)+REAL(ix2,RKIND)*rm2)*rm1

      return
      end
