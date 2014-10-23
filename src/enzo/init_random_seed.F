c     A routine useful for seeding fortran's built-in RNG.
c     This will first try to seed from /dev/urandom but falls back
c     to xoring the current time with the PID of the process.
c     Only call this function ONCE.

      subroutine init_random_seed()
#ifdef __INTEL_COMPILER
      use ifport, only: GETPID
#endif
      use iso_fortran_env, only: int64
      implicit none
#include "fortran_types.def"
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid
      integer :: t

      call random_seed(size = n)
      allocate(seed(n))
c     First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream",
     &     form="unformatted", action="read", status="old",
     &     iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
c     Fallback to XOR:ing the current time and pid. The PID is
c     useful in case one launches multiple instances of the same
c     program in parallel.
         call system_clock(t)
         if (t == 0) then
            call date_and_time(values=dt)
            t = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000
     &      + dt(2) * 31 * 24 * 60 * 60 * 1000
     &      + dt(3) * 24 * 60 * 60 * 1000
     &      + dt(5) * 60 * 60 * 1000
     &      + dt(6) * 60 * 1000 + dt(7) * 1000
     &      + dt(8)
         end if
         pid = GETPID()
         t = ieor(t, int(pid, kind(t)))
         do i = 1, n
            seed(i) = lcg(t)
         end do
      end if
      call random_seed(put=seed)
      contains
c     This simple PRNG might not be good enough for real work, but is
c     sufficient for seeding a better PRNG.
      function lcg(s)
      integer :: lcg
      integer :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_IKIND)
      end if
      s = mod(s * 279470273_IKIND, 4294967291_IKIND)
      lcg = int(mod(s, int(HUGE(0), IKIND)), kind(0))
      end function lcg
      end subroutine init_random_seed
