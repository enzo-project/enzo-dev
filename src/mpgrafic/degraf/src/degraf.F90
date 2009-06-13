
!--------------------------------------------------------------------
!
! The aim of this program is to take as input a collection of 
! grafic files (either density or velocity fields), keep their 
! low frequencies and output smaller corresponding boxes in real
! space. 
! Care must be taken to properly shift the velocity fields in
! real space compared to the density field...
!
! Output file names will be constructed in the following way:
! Suppose that the inputfile name is ic_deltab containing a box
! of linear size 1024, and suppose that the current value of the
! degrading factor is 2, then the output filename will be 
! ic_deltab_512
!
! degramax must be a power of 2, smaller than boxsize/2
! The code will produce degraded versions *down to* boxsize/degramax
! if the recurse flag is set, otherwise it produces only one
! output at resolution boxsize/degramax
!  
!-------------------------------------------------------------------- 

program degraf

  use grafic_types
  use grafic_io
#ifdef ENZO
  use enzo_io
#endif
  use transform
  use lowpass
  implicit none
#ifdef ADD1US
#define  rfftw3d_f77_create_plan  rfftw3d_f77_create_plan_
#define  rfftwnd_f77_destroy_plan rfftwnd_f77_destroy_plan_
#define  rfftwnd_f77_one_real_to_complex rfftwnd_f77_one_real_to_complex_
#define  rfftw3d_f77_mpi_create_plan  rfftw3d_f77_mpi_create_plan_
#define  rfftwnd_f77_mpi_destroy_plan rfftwnd_f77_mpi_destroy_plan_
#define  rfftwnd_f77_mpi rfftwnd_f77_mpi_
#define  rfftwnd_f77_mpi_local_sizes rfftwnd_f77_mpi_local_sizes_
#endif

  ! Code identification
  character(len=*), parameter :: code='DEGRAF'
  character(len=*), parameter :: version='0.1'

  ! Namelist (parameters) declaration
  integer(i4b), parameter :: nfilesmax=20
  character(len=filenamelen) :: namelistfile
  integer(i4b) :: nfiles=1
  character(len=filenamelen), dimension(nfilesmax) :: inputfiles
  integer(i4b) :: degramax=2
  logical(lgt) :: recurse=.false.
  real(sp) :: width_in=1.0
  integer(i4b) :: buffercells=0
  real(dp), dimension(3) :: shift=0 ! Default is no translation
  logical(lgt) :: hanning=.false. ! Default is no hanning correction/application

  namelist /parameters/ nfiles, inputfiles, degramax, recurse, shift, hanning, &
       width_in, buffercells

  ! Mpi stuff
  integer(i4b) :: myid, nproc, ierr

  ! Box header variables
  type(taille) :: htin, htout
  type(cosmo) :: hcosm

  ! Miscellaneous
  integer(i4b) :: narg
  integer(i4b) :: i,j,k,ifiles
  integer(i4b) :: degrafact,degramin
!  integer(i4b), external :: iargc

  integer(i8b) :: plan
#ifdef DOUB
  real(dp), allocatable, dimension(:) :: bufferin, bufferout
#else
  real(sp), allocatable, dimension(:) :: bufferin, bufferout
#endif
  integer(i4b) :: nxin,nyin,nzin,nxout,nyout,nzout
  integer(i4b) :: local_nzin,local_nyin,total_local_sizein
  integer(i4b) :: local_z_startin,local_y_startin
  integer(i4b) :: local_z_starts, local_nzs
  integer(i4b) :: nz_trim, zstart_trim
  integer(i4b) :: topgrid_dims(3)
  real(sp) :: width, boxsize, disp_factor, vunit
  integer :: level, maxlevel, coarse_nx, n, idim, ndims=1
  integer :: border(3)
  logical :: write_paramfile=.false.
  logical, allocatable, dimension(:,:,:) :: mask
  character(len=filenamelen) :: kfile='tempo_kfile'
  character(len=filenamelen) :: outfile=''
  character(len=4) :: dummy=''
  !------------------------------------------------------------------

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,nproc,ierr)
  call init_grafic_io
  if (myid.eq.0) write_paramfile = .true.

  ! Announce program, and read namelist
  if(myid==0) then
     print *,' '
     print *,'         ****************************************************** '
     print *,'                           '//code//' '//version
     print *,'                      Compute degraded grafic file(s)           '
     print *,'                   Using currently ', nproc,' processors        '
     print *,'         ****************************************************** '
     print *,' '
  endif
  
  narg = iargc()
  if (narg /= 1) then
     if (myid==0) then
        print*,' '//code//'> Usage: degraf <namelistfile>'
     endif
     call mpi_finalize(ierr)
     stop
  endif
  
  call getarg(1,namelistfile)
  open(1,file=namelistfile)
  read(1,nml=parameters)
  close(1)

#ifdef ENZO
  ! Check to see if buffercells is a power of two, but not one!  Zero is fine.
  if (myid==0) then
     if ((buffercells.eq.1) .or. & 
          (buffercells.lt.0) .or. &
          (buffercells.ge.2 .and. iand(buffercells, buffercells-1).ne.0)) then
        print*, "buffercells must be a power of two (>=2) or zero!"
        call mpi_finalize(ierr)
        stop
     endif
  endif
  ! Open file to write enzo parameters for cut-and-pasting
  if (write_paramfile) &
       open(enzo_handle, file="enzo.params", form="formatted", status="replace")
#endif

  ! Output namelist content
  
  if (myid==0) then
     print*,' '
     print*,' DEGRAF has been called with the following parameters:'
     print*,' Number of files to process: ',nfiles
     do i=1,nfiles
        print*,' Name of file #',i,' : ',trim(inputfiles(i))
     enddo
     print*,' Value of degramax: ',degramax
     if (recurse .eqv. .true.) then
        print*,' Recursive degradation: true'
     else
        print*,' Recursive degradation: false'
     endif
     if (hanning) then
        print*,' Hanning correction on the big box will be done,'
        print*,' and hanning filter will be applied on smaller boxes'
     endif
  endif

  ! Set up loop parameters depending on recurse
  if (recurse) then
     degramin=2
  else
     degramin=degramax
  endif

  ! Main loop on files to process

  do ifiles=1,nfiles

#ifdef ENZO
     ndims = get_ndims(inputfiles(ifiles))
#endif

     do idim=1,ndims

        if (myid==0) then
           print*,' '
           print*,' **************************************************************'
           if (ndims.gt.1) then
              write(6,'("Processing file ",A,", dimension ", I1)') &
                   trim(inputfiles(ifiles)), idim
           else
              write(6,'("Processing file ",A)') &
                   trim(inputfiles(ifiles))
           endif
        endif

        ! Get file header and content
        call grafic_read_header(inputfiles(ifiles),htin,hcosm)
        nxin=htin%nx
        nyin=htin%ny
        nzin=htin%nz

#ifdef ENZO
        topgrid_dims = (/nxin, nyin, nzin/)
        topgrid_dims = topgrid_dims / degramax
        call set_topgrid_dims(topgrid_dims)
#endif

        ! In the first time through, convert shift from box length
        ! fraction units to physical units (Mpc).
        if (ifiles.eq.1 .and. idim.eq.1) then
           shift(1) = shift(1)*htin%dx*real(nxin,dp)
           shift(2) = shift(2)*htin%dx*real(nyin,dp)
           shift(3) = shift(3)*htin%dx*real(nzin,dp)
        endif

        ! Get local buffer sizes, and prepare fft forward plan
        call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,nxin,nyin&
             &,nzin,fftw_real_to_complex,fftw_estimate)
        call rfftwnd_f77_mpi_local_sizes(plan,local_nzin,local_z_startin&
             &,local_nyin,local_y_startin,total_local_sizein)

        allocate(bufferin(total_local_sizein))
        bufferin=0
        call grafic_read(bufferin,local_nzin,local_z_startin,nzin,nyin,nxin,&
             inputfiles(ifiles),idim-1)

        ! Compute forward transform (in place)
        call mpi_barrier(mpi_comm_world,ierr)
        if (myid==0) then
           print*,' After grafic_read'
        endif

        call init_fftw(total_local_sizein)
        call fft_mpi(plan,bufferin)
        call cleanup_fftw

        call rfftwnd_f77_mpi_destroy_plan(plan)
        if (myid==0) then
           print*,' After fft on input box'
        endif

        ! Now dump big box Fourier transform to disk
        ! This avoids complicated communication between processors
        ! to pass large box coefficients to the small box...
        ! Coefficients of the small box will be read plane by plane
        ! in the large box.
        if (myid==0) then
           call grafic_write_header(kfile,htin,hcosm)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        call grafic_write(bufferin,local_nzin,local_z_startin,nzin,nyin,nxin,&
             kfile,0,1,padding_in=.true.)
        deallocate(bufferin) ! Save memory...

        !
        ! For enzo, write a cropped box for the finest subgrid
        !
        maxlevel = nint(log(real(degramax)) / log(2.0))
#ifdef ENZO
        ! Adjust width to match with a cell boundary at the coarsest
        ! level and a multiple of two in cell number
        coarse_nx = ishft(nxin, -maxlevel)
        width = width_in
        do i = 2, maxlevel
           width = width + real(2*buffercells) / real(ishft(nxin, -i+1))
        enddo
        width = 2*real(ceiling(0.5 * width * coarse_nx)) / coarse_nx

        ! Create output filename
        write(dummy,'(i1)') maxlevel
        outfile=trim(inputfiles(ifiles)) // '.' // trim(adjustl(dummy))

        ! Shift and crop the box
        htout = htin
        call fft_cut(kfile,bufferout,htout,local_nzs,local_z_starts, &
             shift,hanning)
        nz_trim = local_nzin
        zstart_trim = local_z_startin
        call create_mask(mask, htin, htout, degramax, maxlevel, width, &
             buffercells, nz_trim, zstart_trim, write_paramfile)

        ! Write it
        if (myid==0 .and. idim==1) then
           call grafic_write_header(outfile,htout,hcosm)
        endif
        call mpi_barrier(mpi_comm_world,ierr)
        call grafic_write(bufferout,nz_trim,zstart_trim,&
             htout%nz,htout%ny,htout%nx,outfile,idim-1,ndims,mask=mask,&
             level=maxlevel)
        deallocate(bufferout)
#endif

        degrafact = degramin
        level = maxlevel-1
        do while (degrafact <= degramax) ! Loop on resolutions
           htout=htin
           htout%nx=nxin/degrafact
           htout%ny=nyin/degrafact 
           htout%nz=nzin/degrafact
           htout%dx=htin%dx*real(degrafact,sp)
           htout%lx=htin%lx*real(degrafact,sp) ! Case htout%lx=xoff
           htout%ly=htin%ly*real(degrafact,sp) ! Case htout%ly=xoff
           htout%lz=htin%lz*real(degrafact,sp) ! Case htout%lz=xoff

           nxout = htout%nx
           nyout = htout%ny
           nzout = htout%nz

           call fft_cut(kfile,bufferout,htout,local_nzs,local_z_starts, &
                shift,hanning)

           ! Dump to file
#ifdef ENZO
           write(dummy,'(i1)') level
           outfile=trim(inputfiles(ifiles)) // '.' // trim(adjustl(dummy))

           ! Create mask to write the necessary data
           call create_mask(mask, htin, htout, degramax, level, width, &
                buffercells, local_nzs, local_z_starts, write_paramfile)
#else        
           write(dummy,'(i4)') htout%nx
           outfile=trim(inputfiles(ifiles))//'_'//trim(adjustl(dummy))
#endif

           if (myid==0 .and. idim==1) then
              call grafic_write_header(outfile,htout,hcosm)
           endif
           call mpi_barrier(mpi_comm_world,ierr)
           call grafic_write(bufferout,local_nzs,&
                local_z_starts,htout%nz,htout%ny,htout%nx,&
                outfile,idim-1,ndims,mask=mask,level=level)
           if (allocated(bufferout)) deallocate(bufferout)

           ! Go to next resolution
           level = level - 1
           degrafact = degrafact*2
        enddo ! resolution
        write_paramfile = .false.
     enddo ! dimensions
  enddo ! files

  ! Get rid of temporary mode file (kfile)
  if (myid==0) then
     open(99,file=kfile,form='unformatted')
     close(99,status='delete')
  endif

#ifdef ENZO
  ! Write the rest of the parameters and close the file
  if (myid==0) then
     write(enzo_handle, *) ''
     write(enzo_handle, '(A,F10.6)') &
          "CosmologySimulationOmegaBaryonNow = ", hcosm%omegab
     write(enzo_handle, '(A,F10.6)') &
          "CosmologySimulationOmegaCDMNow = ", hcosm%omegam - hcosm%omegab
     write(enzo_handle, '(A,F10.6)') "CosmologyOmegaMatterNow = ", hcosm%omegam
     write(enzo_handle, '(A,F10.6)') "CosmologyOmegaLambdaNow = ", hcosm%omegav
     write(enzo_handle, '(A,F10.6)') "CosmologyHubbleConstantNow = ", hcosm%h0 / 100.0
     write(enzo_handle, '(A,F10.6)') "CosmologyComovingBoxSize = ", &
          htin%dx * htin%nx * hcosm%h0/100.0
     write(enzo_handle, '(A,F10.6)') "CosmologyInitialRedshift = ", &
          1.0/hcosm%astart - 1.0
     close(enzo_handle)
  endif
#endif
  
  call mpi_finalize(ierr)

end program degraf
