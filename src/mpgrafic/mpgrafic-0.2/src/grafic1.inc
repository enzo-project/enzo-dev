
!c  grafic1.inc
!c  User must set parameters then make grafic1
!c
!c  Parameters for offsetting baryon and CDM velocity fields by
!c  0.5 grid spacing.
!c  The origin of the baryon density field is (0,0,0).  The origin
!c  of the baryon velocity field will be offvelb*(0.5,0.5,0.5)*dx
!c  and similarly for offvelc.  offvelb/c should be +/-1 or 0.
	integer offvelb,offvelc
	parameter (offvelb=0,offvelc=0)
!c
!c  Initial density fluctuation amplitude on the ultimate subgrid scale.
!c  For large grids (> 10^7) one may wish to decrease this; for small
!c  refinement grids (< 10^6) and multiple levels one may wish to increase it.
	real sigstart
	parameter (sigstart=0.2)
!c
!c  Top grid mesh size.
!c
	integer np1,np2,np3
	parameter (np1=256,np2=256,np3=256)
!c
!c np[1-3]:      particle lattice size for computing initial conditions.
!c                 All must be divisible by 2.  They need not be equal.
!c  The memory required by grafic1 is np1*np2*np3 words plus a few.
    