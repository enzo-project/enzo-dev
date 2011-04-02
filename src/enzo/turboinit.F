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
c-----------------------------------------------------------------------
c
c  Arguments
c
      integer rank, nbox, in, jn, kn, ig, jg, kg
      real u(in,jn,kn), v(in,jn,kn), w(in,jn,kn)
      REALSUB random_forcing
c
c  Locals
c
      integer i, j, k, imo, nmode
      parameter (nmode=16)
      integer mode(3,nmode)
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
      real phax(nmode)
      data phax/       
     &       4.8827171d0 , 4.5501628d0 , 3.6897256d0  , 5.7606730d0, 
     &       2.0264773d0 , 0.83200777d0, 1.9374901d0 , 0.014175551d0, 
     &       5.1355696d0 , 2.7778759d0 , 2.0290945d0 , 0.66376913d0, 
     &       1.8051250d0 , 3.3130596d0 , 1.0506331d0 , 1.7523085d0/
c
      real phay(nmode)
      data phay/
     &       1.4011313d0, 5.7180996d0 , 3.8207288d0 , 1.0026506d0, 
     &       2.2681668d0, 2.8144622d0 , 0.99058449d0, 2.9458065d0, 
     &       3.9271564d0, 0.89623797d0, 1.8535780d0 , 2.8460610d0, 
     &       1.6346333d0, 3.4661922d0 , 5.5859957d0 , 1.5948143d0/
c
      real phaz(nmode)
      data phaz/
     &       5.6059551d0, 4.1390905d0, 6.2273364d0, 5.9263325d0, 
     &       3.5187488d0, 5.4222918d0, 5.7706189d0, 4.9518018d0, 
     &       4.4614434d0, 5.2936754d0, 5.5074186d0, 2.3949680d0, 
     &       4.5948687d0, 2.2385154d0, 3.1959155d0, 4.4706650d0/
c
c Random Gaussian amplitudes for each mode for seed=12398L, solenoidalized
c
      real amp(nmode,3)
      data amp/ 
     &      0.075595722d0, -1.3572438d0,   0.37845582d0, -0.38310400d0,   ! X
     &      0.11698084d0,  -1.1607968d0,   0.d0,         -0.028096508d0,  ! X
     &      0.d0,           0.d0,         -0.23279878d0,  0.d0,           ! X
     &      0.d0,          -0.87953436d0, -0.60458595d0,  0.d0,           ! X
     &      1.0322379d0,    0.53098691d0, -0.24294342d0, -0.83271527d0,   ! Y
     &     -0.60710335d0,   0.d0,         -0.27813554d0,  0.d0,           ! Y
     &     -1.1801908d0,    0.d0,          0.d0,          0.97667843d0,   ! Y
     &      0.d0,          -0.69450939d0,  0.d0,         -0.60800761d0,   ! Y
     &      1.0182580d0,   -0.96607661d0,  0.21195602d0, -0.60592365d0,   ! Z
     &      0.d0,           0.31490606d0,  0.10941788d0,  0.d0,           ! Z
     &      0.d0,          -1.5361234d0,   0.d0,          0.d0,           ! Z
     &      0.81321216d0,   0.d0,         -0.36861938d0, -0.37148938d0/   ! Z
c
c signs of choice in eqs. (10.6) and (10.7) in Crockett (2005), p.96
c
      real sign1(4)
      data sign1/1.d0,-1.d0,-1.d0, 1.d0/
c
      real sign2(4)
      data sign2/-1.d0,-1.d0, 1.d0, 1.d0/
c
      real  aa, pi, phayy, phazz, k1
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
      aa    = 2.d0*pi/nbox

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
               u(i,j,k) = 0.d0
               v(i,j,k) = 0.d0
               w(i,j,k) = 0.d0
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
     &                 /2.d0/amp(imo,1)/mode(1,imo)
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
	       u(i,j,k) = u(i,j,k)/2.84832d0*random_forcing
               v(i,j,k) = v(i,j,k)/2.84832d0*random_forcing
               w(i,j,k) = w(i,j,k)/2.84832d0*random_forcing
            enddo
         enddo
      enddo
c
      return
      end
c-------------------------------------------------------------------------
