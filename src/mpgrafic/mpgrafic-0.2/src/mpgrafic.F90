!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#include <../config.h>
program grafic1
  !  Generate initial conditions for cosmological N-body integration
  !  as a Gaussian random field.
  !  This version does not do constraints.  Instead, it produces output
  !  compatible with grafic2 for multiscale initial conditions.
  !  Disk use: lun=10 is output file, lun=11 is temp files.
  !
  use grafic_types
  use transform
  use grafic_io
  use normalize
  implicit none
  include 'grafic1.inc'
#ifdef ADD1US
#define  rfftw3d_f77_create_plan  rfftw3d_f77_create_plan_
#define  rfftwnd_f77_destroy_plan rfftwnd_f77_destroy_plan_
#define  rfftwnd_f77_one_real_to_complex rfftwnd_f77_one_real_to_complex_
#define  rfftw3d_f77_mpi_create_plan  rfftw3d_f77_mpi_create_plan_
#define  rfftwnd_f77_mpi_destroy_plan rfftwnd_f77_mpi_destroy_plan_
#define  rfftwnd_f77_mpi rfftwnd_f77_mpi_
#define  rfftwnd_f77_mpi_local_sizes rfftwnd_f77_mpi_local_sizes_
#endif

  !
  !	real f(np1,np2,np3),slice(np1,np2)
#ifdef DOUB
  real(dp), dimension(:), allocatable :: f
#else
  real(sp), dimension(:), allocatable :: f
#endif
  real dx,x1o,x2o,x3o,xoff,fm
  integer irand,iseed,itide,m1s,m2s,m3s,m1off,m2off,m3off
  integer m1t,m2t,m3t,m1offt,m2offt,m3offt
  integer i1,i2,i3,icomp,idim,nrefine,l1,l2,l3,m1,m2,m3
  integer ipad, nxs, nys, nzs, iseeds, irefine
  double precision twopi,sigma,lsigma,dsigma,rombin,sigmadispl
  double precision dsigmadispl
  parameter (twopi=6.283185307179586d0)
  logical hanning
  character(len=128) :: filename, in_file, small_file_name, small_kfile_name
  real astart,omegam,omegav,h0,dladt,fomega,vfact,pbar,pcdm
  real omegab,omegahdm,degen_hdm
  real asig,dxr,dpls,dplus,adp,vunit,boxsize,disp_factor
  integer ierr, myid, status
  integer(i8b) :: plan, iplan
  integer :: local_z_start,local_nz,local_y_start,local_ny
  integer :: nx, ny, nz
  integer :: total_local_size
  integer :: local_z_starts, local_nzs, local_y_starts, local_nys, total_local_sizes
  type(taille) :: headt, shtaille
  type(cosmo) :: headc, shcosmo
  character(len=128) :: outputfile
  common /cosmoparms/ omegam,omegav,h0
  common /cosmoparms2/ omegab,omegahdm,degen_hdm
  common /dsig/ asig,dxr
  external pbar,pcdm,dsigma,dsigmadispl,rombin,dplus,adp
  !
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call init_grafic_io
  if (myid==0) then
     print*
     print*,'Will generate initial conditions on grid of size ', &
          &     np1,np2,np3
     print*
  endif

  ! We follow grafic conventions: slowest varying index is nz
  ! The domain decomposition will be done in this direction
  nx=np1
  ny=np2
  nz=np3


  !
  !  Initialize power spectrum.
  call pini
  !  Initialize grid spacing.
  if (myid==0) then
     print*,'Enter dx (initial particle spacing in Mpc, not Mpc/h)'
     print*,'   or enter -boxlength in Mpc/h'
     print*,'   i.e. dx in Mpc if > 0, or -boxlength in Mpc/h if < 0'
     read(*,*) dx
     print*, dx
  endif
  call mpi_bcast(dx,1,mpi_real,0,mpi_comm_world,ierr)
  if (dx.lt.0.0) then
     boxsize=-dx
     dx=-dx*100./(h0*nx)
     if (myid==0) print*,'  dx=',dx,' Mpc'
  else
     boxsize = dx * nx * (h0/100.0)
  endif

  !  Set astart based on smallest grid spacing.
  if (myid==0) then
     print*,'Enter 1 if these initial conditions will not be ', &
          &    'refined, otherwise'
     print*,'  enter the ultimate refinement factor for the ', &
          &    'smallest grid spacing'
     print*,'  (This is used only to set astart.)'
     read(*,*) nrefine
     print*, nrefine
  endif
  call mpi_bcast(nrefine,1,mpi_integer,0,mpi_comm_world,ierr)
  if (nrefine.lt.1) then
     if (myid==0) print*,'Error! refinement factor must be >= 1'
     call mpi_finalize(ierr)
     stop
  end if
  if (nrefine.gt.1) then
     if (myid==0) print*,'nrefine>1 not implemented yet'
     call mpi_finalize(ierr)
     stop
  end if
  dxr=dx/nrefine
  asig=1.0d0
  sigma=2.0*twopi*rombin(dsigma,0.0d0,0.5d0*twopi/dxr,1.0d-7)
  !  This is sigma at a=1.
  sigma=sqrt(sigma)

!!!!!!!! ALL following !!$ are to remove R. Teyssier's astart and
!!!!!!!! sigma computations....

!!$  sigmadispl=2.0*twopi* &
!!$       &     rombin(dsigmadispl,twopi/dxr/np1,0.5d0*twopi/dxr,1.0d-7)
  !  This is sigma at a=1.
!!$  sigmadispl=sqrt(sigmadispl)
!!$  sigma=MAX(sigma,sigmadispl/dxr)
  !  Normalize so that rms density flutuation=sigstart at starting
  !  redshift scaling back the fluctuations as if they grew like cdm.
  dpls=sigstart/sigma*dplus(1.0,omegam,omegav)
  astart=adp(dpls,omegam,omegav)
  asig=astart
  sigma=2.0*twopi*rombin(dsigma,0.0d0,0.5d0*twopi/dxr,1.0d-7)
!!$  sigma=2.0*twopi*rombin(dsigma,twopi/dxr/np1,0.5d0*twopi/ &
!!$       &                         dxr,1.0d-7)
  sigma=sqrt(sigma)
  if (myid==0) then
     print*,'Scaling initial conditions to starting a=',astart
     print*,'  when sigma(delta) at ultimate refinement scale=' &
          &          ,real(sigma)
  endif
  !sigmadispl=2.0*twopi*rombin(dsigmadispl,0.0d0,0.5d0*twopi/dxr,1.0d-7)
!!$  sigmadispl=2.0*twopi* &
!!$       &   rombin(dsigmadispl,twopi/dxr/np1,0.5d0*twopi/dxr,1.0d-7)
!!$  sigmadispl=sqrt(sigmadispl) 
!!$  if (myid==0) then
!!$     print*,'  when sigma(displ) at ultimate refinement scale=', &
!!$          &    real(sigmadispl),' Mpc '
!!$     print*,'                                                 ', &
!!$          &    real(sigmadispl)/dxr,' dx'
!!$     print*
!!$  endif

!!!!!!!! END OF astart/sigma modifications. Code should be back to
!!!!!!!! original grafic1 starting conditions....


  !  velocity (proper km/s) =  Displacement (comoving Mpc at astart) * vfact.
  !  vfact = dln(D+)/dtau where tau=conformal time.
  vfact=fomega(astart,omegam,omegav)*h0*dladt(astart,omegam,omegav)/astart
  !
  !  Now set output parameters.  There are two cases:
  !  hanning=T: no further refinement.
  !  hanning=F: prepare for further refinement.
  !
  if (myid==0) then
     print*,'Enter 0 for final output or 1 for further refinement'
     read(*,*) irefine
     if (irefine.ne.0) irefine = 1
     print*, irefine
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(irefine,1,mpi_integer,0,mpi_comm_world,ierr)

#ifndef ENZO  
  if (irefine.ne.1) then
     if (myid==0) print*,'Further refinement (input to grafic2) not implemented yet'
     call mpi_finalize(ierr)
     stop
  endif
#endif

!!$  if (irand.eq.1) then
!!$     hanning=.false.
!!$     if (myid==0) print*,'Enter size (n1c,n2c,n3c) of the level-1 subgrid ', &
!!$          &      'to be extracted.'
!!$     if (myid==0) print*,'  Sizes must be even numbers and be no larger than ', &
!!$          &      int(0.5*np1),int(0.5*np2),int(0.5*np3)
!!$     read(*,*) m1s,m2s,m3s
!!$     if (mod(m1s,2).ne.0.or.mod(m2s,2).ne.0.or.mod(m3s,2).ne.0) then
!!$        if (myid==0) print*,'Error! Sizes must be even numbers!'
!!$        call mpi_finalize(ierr)
!!$        stop
!!$     end if
!!$     if (2*m1s.gt.np1.or.2*m2s.gt.np2.or.2*m3s.gt.np3) then
!!$        if (myid==0) print*,'Error!  Subgrid is too large'
!!$        call mpi_finalize(ierr)
!!$        stop
!!$     end if
!!$     if (myid==0) then
!!$        print*,'Enter offset of level-1 subgrid (m1off,m2off,m3off).'
!!$        print*,'  Offsets are relative to level-0 grid corner.'
!!$        print*,'  Offsets may be positive or negative, with absolute ', &
!!$             & ' values no larger than',int(0.5*np1),int(0.5*np2),int(0.5*np3)
!!$     endif
!!$     read(*,*) m1off,m2off,m3off
!!$     !  Now get coordinates for tidal volume.
!!$     if (myid==0) print*,'Enter size (m1t,m2t,m3t) of the final subvolume, ', &
!!$          &      'in units of top grid spacing'
!!$     read(*,*) m1t,m2t,m3t
!!$     if (m1t.gt.m1s.or.m2t.gt.m2s.or.m3t.gt.m3s) then
!!$        if (myid==0) print*,'Error! Final subvolume cannot be larger than ', &
!!$             &        'level-1 subgrid'
!!$        call mpi_finalize(ierr)
!!$        stop
!!$     end if
!!$     if (myid==0) then
!!$        print*,'Enter offset of final subvolume (m1offt,m2off2,m3offt).'
!!$        print*,'  Offsets are relative to level-0 grid corner.'
!!$        print*,'  Final subvolume must lie within level-1 subgrid'
!!$     endif
!!$     read(*,*) m1offt,m2offt,m3offt
!!$     if (m1offt.lt.m1off.or.m1offt+m1t.gt.m1off+m1s.or. &
!!$          &        m2offt.lt.m2off.or.m2offt+m2t.gt.m2off+m2s.or. &
!!$          &        m3offt.lt.m3off.or.m3offt+m3t.gt.m3off+m3s) then
!!$        if (myid==0) print*,'Error! Final subvolume isn''t contained within ', &
!!$             &             'level-1 subgrid'
!!$        call mpi_finalize(ierr)
!!$        stop
!!$     end if
!!$     !  Coordinates of the subgrid corner.
!!$     x1o=m1off*dx
!!$     x2o=m2off*dx
!!$     x3o=m3off*dx
!!$     if (myid==0) then
!!$        print*,'Grafic1 will output level-1 subgrid as follows:'
!!$        print*,'  xmin, xmax = ',x1o,x1o+m1s*dx
!!$        print*,'  ymin, ymax = ',x2o,x2o+m2s*dx
!!$        print*,'  zmin, zmax = ',x3o,x3o+m3s*dx
!!$        print*,'Grafic1 will compute tides assuming final subvolume:'
!!$        print*,'  xmin, xmax = ',m1offt*dx,(m1offt+m1t)*dx
!!$        print*,'  ymin, ymax = ',m2offt*dx,(m2offt+m2t)*dx
!!$        print*,'  zmin, zmax = ',m3offt*dx,(m3offt+m3t)*dx
!!$     endif
!!$     !  Open output file and write header.
!!$     open(10,file='grafic2.top',form='unformatted',status='unknown')
!!$     rewind 10
!!$     write(10)  2*m1s,2*m2s,2*m3s,dx,x1o,x2o,x3o,m1t,m2t,m3t, &
!!$          &      m1offt-m1off,m2offt-m2off,m3offt-m3off,hanning,astart, &
!!$          &      omegam,omegav,h0
!!$     if (myid==0) print*,'Output file is grafic2.top; use as input to subsequent', &
!!$          &           ' run of grafic2'
!!$  else
  hanning=.false.
  m1s=np1 !Unused here...
  m2s=np2 !Unused here...
  m3s=np3 !Unused here...
  x1o=0.0
  x2o=0.0
  x3o=0.0
  if (myid==0) then
     print*,'Setting output grid to ',m1s,m2s,m3s
     print*,'Enter <RETURN> to skip output grid size'
     read(*,*)
     print*,'Enter <RETURN> to skip output grid offset'
     read(*,*)
!     read(*,*) x1o, x2o, x3o
     print*,'Enter <RETURN> to skip final grid size'
     read(*,*)
     print*,'Enter <RETURN> to skip final grid offset'
     read(*,*)
#ifdef ENZO
     print*,'Will produce output files GridDensity, GridVelocities, ', &
          'ParticlePositions, ParticleVelocities'
#else
     print*,'Will produce output files ic_deltab, ic_velb[xyz],', &
          &           ' ic_velc[xyz]'
#endif
  endif
!!$end if
!  call mpi_bcast(x1o, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
!  call mpi_bcast(x2o, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
!  call mpi_bcast(x3o, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
!  ! input is in box length fraction units.  Convert to physical units
!  x1o = x1o * boxsize / (h0/100.0)
!  x2o = x2o * boxsize / (h0/100.0)
!  x3o = x3o * boxsize / (h0/100.0)

  !  Set parameters for subgrid noise.
  if (myid==0) then
     print*
     print*,'Subgrid white noise:'
     print*,'  Choose irand (1 or 2) from the following list:'
     !	print*,'    irand=0 to generate new noise and don''t save it'
     print*,'    irand=1 to generate new noise and save to file'
     print*,'    irand=2 to read noise from existing file'
     print*,'Enter irand'
     read(*,*) irand
     print*, irand
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(irand,1,mpi_integer,0,mpi_comm_world,ierr)
  if (irand.lt.0.or.irand.gt.2) then
     if (myid==0) print*,'Illegal value of irand'
     call mpi_finalize(ierr)
     stop
  end if

  if (myid==0) then 
     print*,'  Enter random number seed (9-digit integer, ignored ', &
          &    'if irand=2)'
     read(*,*) iseed
     print*, iseed
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(iseed,1,mpi_integer,0,mpi_comm_world,ierr)
  filename=''
  if (myid==0) then 
     print*,'  Enter filename of white noise file (or <RETURN> ', &
          &    'if irand=0)'
     read(*,'(a)') filename
     print*,filename
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  ! filename='white.dat'
  ! filename is needed by other cpus (to read or write)
  call mpi_bcast(filename,128,mpi_character,0,mpi_comm_world,ierr)
  ! print*,'in mpgrafic, filename is ','***'//trim(filename)//'***'

  ! Ask if padding of small k modes is required
  if (myid==0) then
     print*,' Enter 1 if padding of large scales is required, 0 otherwise'
     read(*,*) ipad
     print*, ipad
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(ipad,1,mpi_integer,0,mpi_comm_world,ierr)
  small_file_name=''
  if (myid==0) then
     print*,' Enter padding file name (white noise format, ignored if ipad=0)'
     read(*,'(a)') small_file_name
     print*, small_file_name
  endif
  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_bcast(small_file_name,128,mpi_character,0,mpi_comm_world,ierr)

  ! Compute small k modes from small file and store them, if needed
  small_kfile_name='ksmall'
  if (ipad==1) then
     if (myid==0) then
        print*,' Will compute the large scales from file ', trim(small_file_name)
        print*,' Beware that the large scales are assumed to be prewhitened, '
        print*,' and that the constraints are assumed to be imposed with the same '
        print*,' power spectrum: cosmologies of small and large k modes need to be '
        print*,' consistent...'
     endif
!!$     call grafic_read_header(small_file_name,shtaille,shcosmo)
!!$     nxs=shtaille%nx
!!$     nys=shtaille%ny
!!$     nzs=shtaille%nz
     call grafic_read_header_white(small_file_name,nxs,nys,nzs,iseeds)
     call rfftw3d_mpi_create_plan(plan,MPI_COMM_WORLD,nxs,nys,nzs, &
          & FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE)
     call rfftwnd_mpi_local_sizes(plan,local_nzs,local_z_starts, &
          & local_nys, local_y_starts, total_local_sizes)
     call init_fftw(total_local_sizes)
     allocate(f(total_local_sizes))
     call grafic_read(f,local_nzs,local_z_starts,nzs,nys,nxs,small_file_name,&
          white_in=.true.)
     ! Make sure avg is 0 and rms is 1
     call mpnorm(local_nzs,nzs,nys,nxs,total_local_sizes,f)
     ! Normalise a la Bertschinger
     f = f/sqrt(real(nxs*nys*nzs,dp))
     call fft_mpi(plan,f)
     ! Scale small grid power to large grid
     ! f = f*real(np1*np2,dp)*np3/real(nxs*nys*nzs,dp)
     call rfftwnd_mpi_destroy_plan(plan)
     call cleanup_fftw
     shtaille%nx=nxs
     shtaille%ny=nys
     shtaille%nz=nzs
     if (myid==0) call grafic_write_header(small_kfile_name,shtaille, &
          & shcosmo)
     call mpi_barrier(MPI_COMM_WORLD,ierr)
     call grafic_write(f,local_nzs,local_z_starts,nzs,nys,nxs,small_kfile_name, &
          & padding_in=.true.)
     deallocate(f)
  endif

  ! Compute local buffer size from total grid size, and prepare fftw
  ! plans. Beware that domain decomposition is done along nx...
  ! (Slowest varying dimension in the array)

  call rfftw3d_mpi_create_plan(plan,MPI_COMM_WORLD,nx,ny,nz, &
       FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE)
  call rfftw3d_mpi_create_plan(iplan,MPI_COMM_WORLD,nx,ny,nz, &
       FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE)
  call rfftwnd_mpi_local_sizes(plan,local_nz,local_z_start,local_ny, &
       local_y_start,total_local_size)
  call init_fftw(total_local_size)

  if (total_local_size <= 0) then
	if (myid==0) then
	  print*,'Beware: size of local buffers does not fit in 4 bits !'
          print*,'Use with a bigger number of processes ...'
	  print*,'Aborting'
        endif
     call mpi_finalize(ierr)
     stop

  endif  
!
  !  Outer loop over components to be refined.
  !
  do icomp=0,12
     !  0: baryon density.
     !  1,2: inner,outer baryon x-velocity.
     !  3,4: inner,outer baryon y-velocity.
     !  5,6: inner,outer baryon z-velocity.
     !  7,8: inner,outer CDM x-velocity.
     !  9,10: inner,outer CDM y-velocity.
     !  11,12: inner,outer CDM z-velocity.
     !
     !  If not generating grafic2.dat for further refinement, and
     !    if icomp=odd, then skip.
     if (mod(icomp,2).eq.1) cycle
     !
     !  Optional half-grid offset to velocity fields.
     if (icomp.gt.0.and.icomp.le.6) then
        xoff=0.5*dx*offvelb
     else if (icomp.gt.6) then
        xoff=0.5*dx*offvelc
     else
        xoff=0.0
     end if
     !  Don't need to recompute random numbers.
     if (icomp.gt.0.and.irand.eq.1) irand=2
     !
     !  Generate random fields.  idim=0,1,2,3 for density, displacement. 
     itide=0
     if (mod(icomp,2).eq.0) then
        idim=icomp/2
        !  Outer case.
        ! if (.not.hanning.and.idim.gt.0) itide=1
     else
        !  Inner case.
        idim=(icomp+1)/2
        ! if (.not.hanning.and.idim.gt.0) itide=-1
     end if
     allocate(f(total_local_size),stat=status)
     if (status /= 0) then
        print*,'Could not allocate array, aborting'
        call mpi_finalize(ierr)
        stop
     endif
     ! First Feed the headers...
     headt%nx=np1
     headt%ny=np2
     headt%nz=np3
     headt%dx=dx
     headt%lx=x1o+xoff
     headt%ly=x2o+xoff
     headt%lz=x3o+xoff
     headc%astart=astart
     headc%omegam=omegam
     headc%omegav=omegav
     headc%h0=h0
     headc%vfact=vfact

     if (idim.le.3) then
        call mpic4(idim,irand,iseed,itide,m1t,m2t,m3t,m1offt,m2offt, &
             &     m3offt,hanning,filename,astart,pbar,dx,xoff,f,fm, &
             &     plan,iplan,local_nz,local_z_start,total_local_size, &
             &     headt,headc,small_kfile_name,ipad)
     else
        idim=idim-3
        call mpic4(idim,irand,iseed,itide,m1t,m2t,m3t,m1offt,m2offt, &
             &     m3offt,hanning,filename,astart,pcdm,dx,xoff,f,fm, &
             &     plan,iplan,local_nz,local_z_start,total_local_size, &
             &     headt,headc,small_kfile_name,ipad)
     end if
     !
     !  Prepare data for output.
     !
     if (.true.) then
        ! COMMENT OUT THIS LINE TO USE WITH RAMSES
        !if (idim > 0) f = f*vfact
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef ENZO
        select case(idim)
        case(0)   ! baryon density
           f = (1+f) * (omegab/omegam)
!           where (f <= 0)
!              f = 1e-2 * (omegab/omegam)
!           endwhere
        case(1:3) ! baryon/DM velocity
           vunit = vfact/(1.225e2 * boxsize * sqrt(omegam/astart))
           f = f*vunit
        end select
#endif
        !  Output files.
        if (icomp.eq.0) then
           outputfile='ic_deltab'
        else if (icomp.eq.2) then
           outputfile='ic_velbx'
        else if (icomp.eq.4) then
           outputfile='ic_velby'
        else if (icomp.eq.6) then
           outputfile='ic_velbz'
        else if (icomp.eq.8) then
           outputfile='ic_velcx'
        else if (icomp.eq.10) then
           outputfile='ic_velcy'
        else if (icomp.eq.12) then
           outputfile='ic_velcz'
        end if
        if (myid==0) then
           headt%nx=np1
           headt%ny=np2
           headt%nz=np3
           headt%dx=dx
           headt%lx=x1o+xoff
           headt%ly=x2o+xoff
           headt%lz=x3o+xoff
           headc%astart=astart
           headc%omegam=omegam
           headc%omegab=omegab
           headc%omegav=omegav
           headc%h0=h0
           headc%vfact=vfact
           call grafic_write_header(outputfile,headt,headc)
        endif
        call mpi_barrier(mpi_comm_world,ierr)

        ! If using enzo and there is no more refinement (irefine ==
        ! 0), write particles as a 1D array in HDF and compute and
        ! output positions.

!        if (irefine.eq.0 .and. icomp.ge.8) then
!
!           ! (vfact/vunit) converts back to km/s.
!           ! (1/vfact/(boxsize/(h0/100))) converts km/s to a particle
!           ! displacement
!           disp_factor = 1.0 / vunit / (boxsize / (h0/100.0))
!           print*, f(1:5)
!           print*, local_nz, local_z_start, nz, ny, nx, disp_factor
!           print*, trim(outputfile)
!           print*, size(f)
!           call write_particles(f,local_nz,local_z_start,nz,ny,nx,outputfile,&
!                disp_factor)
!        else
        call grafic_write(f,local_nz,local_z_start,nz,ny,nx,outputfile)
!        endif
     endif ! subgrid is not implemented...

!!$  else
!!$     !  hanning=F, extract next-level subgrid and append to grafic2.dat.
!!$     !  First surround subvolume with 1/2-size buffer and wrap periodically.
!!$     sigma=0.0
!!$     do m3=1,2*m3s
!!$        i3=m3+m3off
!!$        !  Periodic boundary conditions on top grid.
!!$        if (i3.lt.1) i3=i3+np3
!!$        if (i3.gt.np3) i3=i3-np3
!!$        l3=i3
!!$        if (m3.gt.1.5*m3s) l3=l3-2*m3s
!!$        if (l3.lt.1) l3=l3+np3
!!$        do m2=1,2*m2s
!!$           i2=m2+m2off
!!$           if (i2.lt.1) i2=i2+np2
!!$           if (i2.gt.np2) i2=i2-np2
!!$           l2=i2
!!$           if (m2.gt.1.5*m2s) l2=l2-2*m2s
!!$           if (l2.lt.1) l2=l2+np2
!!$           do m1=1,2*m1s
!!$              i1=m1+m1off
!!$              if (i1.lt.1) i1=i1+np1
!!$              if (i1.gt.np1) i1=i1-np1
!!$              l1=i1
!!$              if (m1.gt.1.5*m1s) l1=l1-2*m1s
!!$              if (l1.lt.1) l1=l1+np1
!!$              slice(m1,m2)=f(l1,l2,l3)
!!$              sigma=sigma+f(l1,l2,l3)**2
!!$           end do
!!$        end do
!!$        write(10) ((slice(m1,m2),m1=1,2*m1s),m2=1,2*m2s)
!!$     end do
!!$     sigma=sqrt(sigma/(8*m1s*m2s*m3s))
!!$     if (myid==0) print*,'After extraction, component ',icomp, &
!!$          &      ' has subvolume RMS=',real(sigma)
!!$  end if
     !  End loop over icomp.
     deallocate(f)
  enddo
  !
!!$if (.not.hanning) close(10)
  call rfftwnd_mpi_destroy_plan(plan)
  call rfftwnd_mpi_destroy_plan(iplan)
  call cleanup_fftw

  call mpi_finalize(ierr)
  stop
end program grafic1

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dsigma(ak)
  !  This function calculates the variance of density with a Hanning
  !  filter at the grid Nyquist frequency.
  !
  implicit none
  double precision dsigma,ak
  real asig,dx,p
  common /dsig/ asig,dx
  external p
  !
  if (ak.le.0.0d0.or.ak.gt.3.1415926535898d0/dx) then
     dsigma=0.0d0
     return
  end if
  !  Hanning filter.
  dsigma=ak*ak*p(real(ak),asig)*cos(0.5*ak*dx)
  !
  return
end function dsigma

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function dsigmadispl(ak)
  !  This function calculates the variance of displacement with a Hanning
  !  filter at the grid Nyquist frequency.
  !
  implicit none
  double precision dsigmadispl,ak
  real asig,dx,p
  common /dsig/ asig,dx
  external p
  !
  if (ak.le.0.0d0.or.ak.gt.3.1415926535898d0/dx) then
     dsigmadispl=0.0d0
     return
  end if
  !  Hanning filter. 
  dsigmadispl=ak*ak*p(real(ak),asig)/ak/ak/3.*cos(0.5*ak*dx)
  !
  return
end function dsigmadispl
