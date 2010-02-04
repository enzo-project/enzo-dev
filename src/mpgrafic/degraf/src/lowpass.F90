module lowpass

  use grafic_types
  use grafic_io
  use transform
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

contains

  subroutine fft_cut(kfile,bufout,htout,local_nzs,local_z_starts&
       &,shift,hanning)
    
    ! Arguments
    character(len=filenamelen), intent(in) :: kfile
#ifdef DOUB
    real(dp), allocatable, dimension(:), intent(out) :: bufout
#else
    real(sp), allocatable, dimension(:), intent(out) :: bufout
#endif
    type(taille), intent(in) :: htout
    integer, intent(out) :: local_nzs, local_z_starts
    real(dp), dimension(3) :: shift
    logical(lgt), intent(in) :: hanning

    ! Local variables
    type(taille) :: htin
    type(cosmo) :: hcosm
    integer :: total_local_sizes
    integer :: local_y_starts, local_nys
    integer :: nread, nstart
    integer :: nx,ny,nz,n2x,nxs,nys,nzs,n2xs
    integer :: i,j,k,kk
    integer :: myid, ierr
    integer(i8b) :: plan,index_big,index_small
#ifdef DOUB
    real(dp), allocatable, dimension(:) :: bufin
#else
    real(sp), allocatable, dimension(:) :: bufin
#endif
    logical(lgt) :: debut, fin
    real(dp) :: xoff, dx, xoffs, dxs, dk1, dk2, dk3, akmax
    real(dp) :: xtrs, ytrs, ztrs
    real(dp) :: theta_big, theta_small, ak1, ak2, ak3, ak, ak33, ak23
    real(dp) :: ak1b, ak2b, ak3b, dk1b, dk2b, dk3b, akmaxb
    real(dp) :: norm
    complex(dpc) :: ctemp,z
    integer :: nreadtot
    logical :: ignorenyq
    real(dp) :: halfgridshift
    real(dp), parameter :: tol=1.0d-10

    call mpi_comm_rank(mpi_comm_world,myid,ierr)

    nxs=htout%nx
    nys=htout%ny
    nzs=htout%nz
    n2xs=2*(nxs/2+1)
    dxs=real(htout%dx,kind=dp)
    xoffs=real(htout%lx,kind=dp)
    xtrs=xoffs+0.5*dxs+shift(1) ! Add extra translation
    ytrs=xoffs+0.5*dxs+shift(2) ! Add extra translation
    ztrs=xoffs+0.5*dxs+shift(3) ! Add extra translation

    call rfftw3d_f77_mpi_create_plan(plan,mpi_comm_world,nxs,nys,nzs&
         &,fftw_complex_to_real,fftw_estimate)
    call rfftwnd_f77_mpi_local_sizes(plan,local_nzs,local_z_starts&
         &,local_nys,local_y_starts,total_local_sizes)
    if (.not. allocated(bufout)) allocate(bufout(total_local_sizes))

    call grafic_read_header(kfile,htin,hcosm)
    nx = htin%nx
    ny = htin%ny
    nz = htin%nz
    n2x= 2*(nx/2+1)
    dx=real(htin%dx,kind=dp)
    xoff=real(htin%lx,kind=dp)
    if (xoff.eq.xoffs) then
       xoff=0
       xoffs=0
    endif

    dk1=2.0*PI/(nxs*dxs)
    dk2=2.0*PI/(nys*dxs)
    dk3=2.0*PI/(nzs*dxs)
    akmax=2.0*PI/dxs ! Dimensioned k vectors are those of the
    ! degraded box...

    ! Check if xoff corresponds to half a mesh size of the big box
    ! modulo an integer number of mesh sizes. If this is the case,
    ! the big nyquist modes have been killed, do not try to resuscitate them...
    ignorenyq=.false.
    !!!!!halfgridshift = mod(xoff/dx-0.5,1.0)
    !halfgridshift=(xoff/dx-0.5)-nint(xoff/dx-0.5)
    !if (abs(halfgridshift).gt.tol) ignorenyq=.true.

    ! We have to separate the positive and negative frequency parts
    ! in kx because they are at the opposite borders of the big cube

    debut=.false.
    fin=.false.
    ! First element is a positive (small) kx frequency (including Nyquist)
    if (local_z_starts+1 <= nzs/2+1) debut=.true.
    ! Last element of segment is a (small) kx negative frequency
    if (local_z_starts+local_nzs >= nzs/2+2) fin=.true.

    nreadtot=0
    if (debut) then
       ! Compute number of planes to read and starting point in big box
       nread = min(local_z_starts+local_nzs,nzs/2+1)-local_z_starts
       nreadtot=nread
       nstart = local_z_starts
       allocate(bufin(nread*ny*n2x))
       call grafic_read(bufin,nread,nstart,nz,ny,nx,kfile,0,padding_in=.true., &
            serial_in=.true.)
       do k=1,nread
          ak3=(k+local_z_starts-1)*dk3
          if (k+local_z_starts.gt.nzs/2) ak3=ak3-akmax
          ak33=ak3*ak3
          do j=1,nys/2+1 ! Includes Nyquist, positive ky only
             ak2=(j-1)*dk2
             !if (j.gt.nys/2) ak2=ak2-akmax
             ak23=ak2*ak2+ak33
             do i=1,nxs/2+1 ! Includes Nyquist. Complex.
                ak1=(i-1)*dk1
                ak=sqrt(ak1*ak1+ak23)
                index_big = int((k-1)*ny+j-1,i8b)*n2x+2*i-1
                index_small = int((k-1)*nys+j-1,i8b)*n2xs+2*i-1
#ifdef DOUB
                ctemp=cmplx(bufin(index_big),bufin(index_big+1),kind=dpc)
#else
                ctemp=cmplx(real(bufin(index_big),kind=dp),real(bufin(index_big+1),kind=dp))
#endif
                ! Correct for hanning of big box, and add hanning of
                ! small box
                if (hanning .and. nxs<nx) then
                   if (ak .ge. akmax/2) then
                      ctemp = 0.0
                   else
                      ctemp = ctemp *cos(0.5*ak*dxs)/cos(0.5*ak*dx)
                   endif
                endif
                ! Compute shift back
                theta_big=0.0
                if (k+local_z_starts.ne.nz/2+1) theta_big=theta_big+ak3*xoff
                if (j.ne.ny/2+1) theta_big=theta_big+ak2*xoff
                if (i.ne.nx/2+1) theta_big=theta_big+ak1*xoff
                z=cmplx(cos(theta_big),sin(-theta_big)) ! Shift back
                if (k+local_z_starts.eq.nz/2+1 .and. .not.ignorenyq) z=z/cos(ak3*xoff)
                if (j.eq.ny/2+1 .and. .not.ignorenyq) z=z/cos(ak2*xoff)
                if (i.eq.nx/2+1 .and. .not.ignorenyq) z=z/cos(ak1*xoff)
                ! Compute shift forward
                theta_small=0
                if (k+local_z_starts.ne.nzs/2+1) theta_small=theta_small+ak3*ztrs
                if (j.ne.nys/2+1) theta_small=theta_small+ak2*ytrs
                if (i.ne.nxs/2+1) theta_small=theta_small+ak1*xtrs
                z=z*cmplx(cos(theta_small),sin(theta_small))
                if (k+local_z_starts.eq.nzs/2+1) z=z*cos(ak3*ztrs)
                if (j.eq.nys/2+1) z=z*cos(ak2*ytrs)
                if (i.eq.nxs/2+1) z=z*cos(ak1*xtrs)
                ctemp = ctemp*z ! Correct phases
#ifdef DOUB
                bufout(index_small) = real(ctemp)
                bufout(index_small+1) = aimag(ctemp)
#else
                bufout(index_small) = real(ctemp,kind=sp)
                bufout(index_small+1) = real(aimag(ctemp),kind=sp)
#endif
             enddo
          enddo
          do j=nys/2+2,nys ! Negative ky modes
             ak2=(j-1)*dk2 ! j > nys/2
             ak2=ak2-akmax
             ak23=ak2*ak2+ak33
             do i=1,nxs/2+1 ! Includes Nyquist. Complex.
                ak1=(i-1)*dk1
                ak=sqrt(ak1*ak1+ak23)
                index_big = int((k-1)*ny+j+ny-nys-1,i8b)*n2x+2*i-1
                index_small = int((k-1)*nys+j-1,i8b)*n2xs+2*i-1
#ifdef DOUB
                ctemp=cmplx(bufin(index_big),bufin(index_big+1),kind=dpc)
#else
                ctemp=cmplx(real(bufin(index_big),kind=dp),real(bufin(index_big+1),kind=dp))
#endif
                ! Correct for hanning of big box, and add hanning of
                ! small box
                if (hanning .and. nxs<nx) then
                   if (ak .ge. akmax/2) then
                      ctemp = 0.0
                   else
                      ctemp = ctemp *cos(0.5*ak*dxs)/cos(0.5*ak*dx)
                   endif
                endif
                ! Compute shift back
                theta_big=ak2*xoff ! Never touch Nyquist
                if (k+local_z_starts.ne.nz/2+1) theta_big=theta_big+ak3*xoff
                if (i.ne.nx/2+1) theta_big=theta_big+ak1*xoff
                z=cmplx(cos(theta_big),sin(-theta_big)) ! Shift back
                if (k+local_z_starts.eq.nz/2+1 .and. .not. ignorenyq) z=z/cos(ak3*xoff)
                if (i.eq.nx/2+1 .and. .not.ignorenyq) z=z/cos(ak1*xoff)
                ! Compute shift forward
                theta_small=ak2*ytrs ! Never touch Nyquist
                if (k+local_z_starts.ne.nzs/2+1) theta_small=theta_small+ak3*ztrs
                if (i.ne.nxs/2+1) theta_small=theta_small+ak1*xtrs
                z=z*cmplx(cos(theta_small),sin(theta_small))
                if (k+local_z_starts.eq.nzs/2+1) z=z*cos(ak3*ztrs)
                if (i.eq.nxs/2+1) z=z*cos(ak1*xtrs)
                ctemp = ctemp*z
#ifdef DOUB
                bufout(index_small) = real(ctemp)
                bufout(index_small+1) = aimag(ctemp)
#else
                bufout(index_small) = real(ctemp,kind=sp)
                bufout(index_small+1) = real(aimag(ctemp),kind=sp)
#endif
             enddo
          enddo
       enddo
       deallocate(bufin)
    endif
    
    if (fin) then ! Covers the negative kx modes
       ! Compute number of planes to read and starting point in big
       ! box
       nread = min(local_nzs,local_z_starts+local_nzs-nzs/2-1)
       nreadtot = nreadtot+nread
       nstart = max(local_z_starts+nz-nzs,nz-nzs/2+1)
       allocate(bufin(nread*ny*n2x))
       call grafic_read(bufin,nread,nstart,nz,ny,nx,kfile,0,padding_in=.true.,&
            serial_in=.true.)
       do k=1,nread
          kk=k+local_nzs-nread ! Fill second part of the small box slice...
          ak3=(kk+local_z_starts-1)*dk3
          if (kk+local_z_starts.gt.nzs/2) ak3=ak3-akmax
          ak33=ak3*ak3
          do j=1,nys/2+1
             ak2=(j-1)*dk2
             !if (j.gt.nys/2) ak2=ak2-akmax
             ak23=ak2*ak2+ak33
             do i=1,nxs/2+1 ! Complex
                ak1=(i-1)*dk1
                ak=sqrt(ak1*ak1+ak23)
                index_big = int((k-1)*ny+j-1,i8b)*n2x+2*i-1
                index_small = int((kk-1)*nys+j-1,i8b)*n2xs+2*i-1
#ifdef DOUB
                ctemp=cmplx(bufin(index_big),bufin(index_big+1),kind=dpc)
#else
                ctemp=cmplx(real(bufin(index_big),kind=dp),real(bufin(index_big+1),kind=dp))
#endif
                ! Correct for hanning of big box, and add hanning of
                ! small box
                if (hanning .and. nxs<nx) then
                   if (ak .ge. akmax/2) then
                      ctemp = 0.0
                   else
                      ctemp = ctemp *cos(0.5*ak*dxs)/cos(0.5*ak*dx)
                   endif
                endif
                ! Compute shift back
                theta_big=ak3*xoff ! Never touch Nyquist
                if (j.ne.ny/2+1) theta_big=theta_big+ak2*xoff
                if (i.ne.nx/2+1) theta_big=theta_big+ak1*xoff
                z=cmplx(cos(theta_big),sin(-theta_big)) ! Shift back
                if (j.eq.ny/2+1 .and. .not.ignorenyq) z=z/cos(ak2*xoff)
                if (i.eq.nx/2+1 .and. .not.ignorenyq) z=z/cos(ak1*xoff)
                ! Compute shift forward
                theta_small=ak3*ztrs ! Never touch Nyquist
                if (j.ne.nys/2+1) theta_small=theta_small+ak2*ytrs
                if (i.ne.nxs/2+1) theta_small=theta_small+ak1*xtrs
                z=z*cmplx(cos(theta_small),sin(theta_small))
                if (j.eq.nys/2+1) z=z*cos(ak2*ytrs)
                if (i.eq.nxs/2+1) z=z*cos(ak1*xtrs)
                ctemp = ctemp*z
#ifdef DOUB
                bufout(index_small) = real(ctemp)
                bufout(index_small+1) = aimag(ctemp)
#else
                bufout(index_small) = real(ctemp,kind=sp)
                bufout(index_small+1) = real(aimag(ctemp),kind=sp)
#endif
             enddo
          enddo
          do j=nys/2+2,nys
             ak2=(j-1)*dk2
             ak2=ak2-akmax ! j > nys/2
             ak23=ak2*ak2+ak33
             do i=1,nxs/2+1 ! Complex
                ak1=(i-1)*dk1
                ak=sqrt(ak1*ak1+ak23)
                index_big = int((k-1)*ny+j+ny-nys-1,i8b)*n2x+2*i-1
                index_small = int((kk-1)*nys+j-1,i8b)*n2xs+2*i-1
#ifdef DOUB
                ctemp=cmplx(bufin(index_big),bufin(index_big+1),kind=dpc)
#else
                ctemp=cmplx(real(bufin(index_big),kind=dp),real(bufin(index_big+1),kind=dp))
#endif
                ! Correct for hanning of big box, and add hanning of
                ! small box
                if (hanning .and. nxs<nx) then
                   if (ak .ge. akmax/2) then
                      ctemp = 0.0
                   else
                      ctemp = ctemp *cos(0.5*ak*dxs)/cos(0.5*ak*dx)
                   endif
                endif
                ! Compute shift back
                theta_big=ak3*xoff ! Never touch Nyquist
                theta_big=theta_big+ak2*xoff ! Never touch Nyquist
                if (i.ne.nx/2+1) theta_big=theta_big+ak1*xoff
                z=cmplx(cos(theta_big),sin(-theta_big)) ! Shift back
                if (i.eq.nx/2+1 .and. .not.ignorenyq) z=z/cos(ak1*xoff)
                ! Compute shift forward
                theta_small=ak3*ztrs ! Never touch Nyquist
                theta_small=theta_small+ak2*ytrs ! Never touch Nyquist
                if (i.ne.nxs/2+1) theta_small=theta_small+ak1*xtrs
                z=z*cmplx(cos(theta_small),sin(theta_small))
                if (i.eq.nxs/2+1) z=z*cos(ak1*xtrs)
                ctemp = ctemp*z
#ifdef DOUB
                bufout(index_small) = real(ctemp)
                bufout(index_small+1) = aimag(ctemp)
#else
                bufout(index_small) = real(ctemp,kind=sp)
                bufout(index_small+1) = real(aimag(ctemp),kind=sp)
#endif
             enddo
          enddo
       enddo
       deallocate(bufin)

    endif

    write(6,'("Proc #",I4," :: nreadtot, local_nzs = ",I8,I8)'), myid, &
         nreadtot, local_nzs

    call init_fftw(total_local_sizes)
    call fft_mpi(plan,bufout)
    ! Normalize the output: normalize fftw, and correct for ratio of
    ! sizes...
    norm = real(nx,dp)*real(ny,dp)*real(nz,dp)
    bufout = bufout/norm
    call cleanup_fftw
    call rfftwnd_f77_mpi_destroy_plan(plan)
  end subroutine fft_cut

end module lowpass
