#include "fortran.def"
#include "phys_const.def"
#include "error.def"
c=======================================================================
c/////////////////////////  SUBROUTINE TURBOINIT \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine turboinit(rank, nbox, u, v, w, in, jn, kn, ig, jg, kg      
     +     ,random_forcing)
c
c  PROVIDES QUASI-ISOTROPIC ZERO-MEAN SOLENOIDAL VELOCITIES FOR LARGE-SCALE  
c  TURBULENCE DRIVING (AND/OR INITIAL CONDITIONS) IN A PERIODIC BOX
c
c  written by: Alexei Kritsuk
c
c  date:       August 2007
c  Modified:   dec. 2007, dcollins: added mach number as a parameter.
c
c  PURPOSE: 
c
c For definitions of Multidimensional DFT see
c http://en.wikipedia.org/wiki/Discrete_Fourier_transform
c
c  EXTERNALS:
c
c  INPUTS:
c     in,jn,kn  - dimensions of field arrays
c     ig,jg,kg  - global zone number (for the left-most zone) for each dimension
c     rank      - dimension of problem (unused)
c     u         - x-velocity field
c     v         - y-velocity field
c     w         - z-velocity field
c     random_forcing - random forcing mach number.
c
c  LOCALS:
c
c-----------------------------------------------------------------------
      implicit NONE
#include "fortran_types.def"
c-----------------------------------------------------------------------
c
c  Arguments
c
      INTG_PREC rank, nbox, in, jn, kn, ig, jg, kg
      R_PREC u(in,jn,kn), v(in,jn,kn), w(in,jn,kn)
      P_PREC random_forcing
c
c  Locals
c
      INTG_PREC i, j, k, imo, nmode
      parameter (nmode=16)
      INTG_PREC mode(3,nmode)
      data    mode/
     &             1,1,1,
     &             -1,1,1,
     &             1,-1,1,
     &             1,1,-1,
     &             0,0,1,
     &             0,1,0,
     &             1,0,0,
     &             0,1,1,
     &             1,0,1,
     &             1,1,0,
     &             0,-1,1,
     &             -1,0,1,
     &             -1,1,0,
     &             0,0,2,
     &             0,2,0,
     &             2,0,0/
c
c A set of randomly selected phases for seed=12398L that provide good isotropy
c Phases are uniformly sampled from [0, 2pi)
c Phases for x, y, and z velocities for each mode
c
      R_PREC phax(nmode)
      data phax/       
     &       4.8827171_RKIND , 4.5501628_RKIND ,  3.6897256_RKIND ,
     &       5.760673_RKIND  , 2.0264773_RKIND ,  0.83200777_RKIND, 
     &       1.9374901_RKIND , 0.014175551_RKIND, 5.1355696_RKIND , 
     &       2.7778759_RKIND , 2.0290945_RKIND ,  0.66376913_RKIND, 
     &       1.805125_RKIND  , 3.3130596_RKIND ,  1.0506331_RKIND , 
     &       1.7523085_RKIND/
c
      R_PREC phay(nmode)
      data phay/
     &       1.4011313_RKIND, 5.7180996_RKIND , 3.8207288_RKIND , 
     &       1.0026506_RKIND, 2.2681668_RKIND , 2.8144622_RKIND , 
     &       0.99058449_RKIND, 2.9458065_RKIND, 3.9271564_RKIND, 
     &       0.89623797_RKIND, 1.853578_RKIND , 2.846061_RKIND, 
     &       1.6346333_RKIND, 3.4661922_RKIND , 5.5859957_RKIND , 
     &       1.5948143_RKIND/
c
      R_PREC phaz(nmode)
      data phaz/
     &       5.6059551_RKIND, 4.1390905_RKIND, 6.2273364_RKIND, 
     &       5.9263325_RKIND, 3.5187488_RKIND, 5.4222918_RKIND, 
     &       5.7706189_RKIND, 4.9518018_RKIND, 4.4614434_RKIND, 
     &       5.2936754_RKIND, 5.5074186_RKIND, 2.394968_RKIND, 
     &       4.5948687_RKIND, 2.2385154_RKIND, 3.1959155_RKIND, 
     &       4.470665_RKIND/
c
c Random Gaussian amplitudes for each mode for seed=12398L, solenoidalized
c
      R_PREC amp(nmode,3)
      data amp/ 
     &      0.075595722_RKIND, -1.3572438_RKIND,   0.37845582_RKIND,      ! X
     &     -0.383104_RKIND,     0.11698084_RKIND, -1.1607968_RKIND,       ! X
     &      0._RKIND,          -0.028096508_RKIND, 0._RKIND,              ! X
     &      0._RKIND,          -0.23279878_RKIND,  0._RKIND,              ! X
     &      0._RKIND,          -0.87953436_RKIND, -0.60458595_RKIND,      ! X
     &      0._RKIND,                                                     ! X
     &      1.0322379_RKIND,    0.53098691_RKIND, -0.24294342_RKIND,      ! Y
     &     -0.83271527_RKIND,  -0.60710335_RKIND,  0._RKIND,              ! Y
     &     -0.27813554_RKIND,   0._RKIND,          -1.1801908_RKIND,      ! Y
     &      0._RKIND,          0._RKIND,          0.97667843_RKIND,       ! Y
     &      0._RKIND,          -0.69450939_RKIND,  0._RKIND,              ! Y
     &     -0.60800761_RKIND,                                             ! Y
     &      1.018258_RKIND,   -0.96607661_RKIND,  0.21195602_RKIND,       ! Z
     &     -0.60592365_RKIND,  0._RKIND,          0.31490606_RKIND,       ! Z
     &      0.10941788_RKIND,  0._RKIND,          0._RKIND,               ! Z
     &     -1.5361234_RKIND,   0._RKIND,          0._RKIND,               ! Z
     &      0.81321216_RKIND,   0._RKIND,         -0.36861938_RKIND,      ! Z
     &     -0.37148938_RKIND/                                             ! Z
c
c signs of choice in eqs. (10.6) and (10.7) in Crockett (2005), p.96
c
      R_PREC sign1(4)
      data sign1/1._RKIND,-1._RKIND,-1._RKIND, 1._RKIND/
c
      R_PREC sign2(4)
      data sign2/-1._RKIND,-1._RKIND, 1._RKIND, 1._RKIND/
c
      R_PREC  aa, pi, phayy, phazz, k1
      parameter (pi=pi_val) 
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////////////
c=======================================================================
c
c  Error checks
c
      if (rank .ne. 3) then
         write(6,*) 'TURBOINIT: Only 3D grids are supported.'
         write(0,*) 'stop_all_cpus in turboinit'
         ERROR_MESSAGE
      endif
c
      if (max(in,jn,kn) .gt. MAX_ANY_SINGLE_DIRECTION) then
         write(6,*) 'TURBOINIT: A grid dimension is too long.'
         write(6,*) '   (increase max_any_single_direction.)'
         write(0,*) 'stop_all_cpus in turboinit'
         ERROR_MESSAGE
      endif
c

c      write(6,*) 'TURBOINIT: ig =',ig,  ' jg =',jg,  ' kg =',kg

c   this is for large-scale force 1<k<2
      aa    = 2._RKIND*pi/nbox

c   this is for the force  8<k<16
c      aa    = 8.0*2.0*pi/nbox
c
c fill-in the velocity arrays
c
      do k=1, kn
         do j=1, jn
            do i=1, in
c
c   fill in 0s first
c
               u(i,j,k) = 0._RKIND
               v(i,j,k) = 0._RKIND
               w(i,j,k) = 0._RKIND
c
c   start with first four modes
c
               do imo=1,4
                  k1 = mode(1,imo)*(i+ig) + 
     &                 mode(2,imo)*(j+jg) + 
     &                 mode(3,imo)*(k+kg) 
                  u(i,j,k) = u(i,j,k) + 
     &                            amp(imo,1)*cos(aa*k1 + phax(imo))
c
c get solenoidal corrections for y- and z-phases of modes with
c k=(1,1,1), (-1,1,1), (1,-1,1), and (1,1,-1)
c 
                  phayy = phax(imo) + sign1(imo)
     &                 *acos((amp(imo,3)**2-amp(imo,1)**2-amp(imo,2)**2)
     &                 /2._RKIND/amp(imo,1)/mode(1,imo)
     &                 /mode(2,imo)/amp(imo,2))
                  v(i,j,k) = v(i,j,k) + 
     &                 amp(imo,2)*cos(aa*k1 + phayy)
                  phazz = phax(imo) + sign2(imo)
     &                 *acos((amp(imo,2)**2-amp(imo,1)**2-amp(imo,3)**2)
     &                 /2.0/amp(imo,1)/mode(1,imo)
     &                 /mode(3,imo)/amp(imo,3))
                  w(i,j,k) = w(i,j,k) + 
     &                 amp(imo,3)*cos(aa*k1 + phazz)
               enddo
c    
c continue with other modes
c
               do imo=5,nmode
                  k1 = mode(1,imo)*(i+ig) + 
     &                 mode(2,imo)*(j+jg) + 
     &                 mode(3,imo)*(k+kg) 
                  u(i,j,k) = u(i,j,k) + 
     &                               amp(imo,1)*cos(aa*k1 + phax(imo))
                  v(i,j,k) = v(i,j,k) + 
     &                               amp(imo,2)*cos(aa*k1 + phay(imo))
                  w(i,j,k) = w(i,j,k) + 
     &                               amp(imo,3)*cos(aa*k1 + phaz(imo))
               enddo
c
c normalize to get rms 3D Mach = 3.0
c
	       u(i,j,k) = u(i,j,k)/2.84832_RKIND*random_forcing
               v(i,j,k) = v(i,j,k)/2.84832_RKIND*random_forcing
               w(i,j,k) = w(i,j,k)/2.84832_RKIND*random_forcing
            enddo
         enddo
      enddo
c
      return
      end
c-------------------------------------------------------------------------
