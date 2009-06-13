!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine splder(y,dy,n)
!       Splder fits a cubic spline to y and returns the first derivatives at
!       the grid points in dy.  dy is equivalent to a 4th-order Pade
!       difference formula for dy/di.
!       
	implicit double precision (a-h,o-z)
	dimension f(100001),gspl(100001),y(n),dy(n)
	common /spline/ gspl
!       
	n1=n-1
	if (n1.gt.100000) &
           write(*,*) 'Spline array overflow!!! n1=',n1,'>100000'
!       Quartic fit to dy/di at boundaries, assuming d3y/di3=0.
	f(1)=(-10.0d0*y(1)+15.0d0*y(2)-6.0d0*y(3)+y(4))/6.0d0
	f(n)=(10.0d0*y(n)-15.0d0*y(n1)+6.0d0*y(n-2)-y(n-3))/6.0d0
!       Solve the tridiagonal system
!       dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1,
!       with dy(1)=f(1), dy(n)=f(n).
	do 10 i=2,n1
	   f(i)=gspl(i)*(3.0d0*(y(i+1)-y(i-1))-f(i-1))
 10	continue
	dy(n)=f(n)
	do 20 i=n1,1,-1
	   dy(i)=f(i)-gspl(i)*dy(i+1)
 20	continue
	return
	end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine spintn(y,dy,z,n)
!       Spintn calculates z(i) = the integral from 1 to i of s(j)dj, for 1 <= i <= n,
!       where s is the spline fit to y, and dy = ds/di.
!       
	implicit double precision (a-h,o-z)
	dimension y(n),dy(n),z(n)

	ztemp=0.0d0
	do 10 i=2,n
	   ztemp1=ztemp+0.5d0*(y(i-1)+y(i))+(dy(i-1)-dy(i))/12.0d0
	   z(i-1)=ztemp
	   ztemp=ztemp1
 10	continue
	z(n)=ztemp
	return
	end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine splini
!       Splini must be called before splder to initialize array g in common.

	implicit double precision (a-h,o-z)
	dimension gspl(100001)
	common /spline/ gspl
	save /spline/
!       
	gspl(1)=0.0d0
	do 10 i=2,100001
	   gspl(i)=1.0d0/(4.0d0-gspl(i-1))
 10	continue
	return
	end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine matinv(a,n,np,y)
!       Invert the double precision nxn matrix a; physical dimensions are np x np.
!       Inverse is returned in y.  WARNING: a is destroyed!
	double precision a(np,np),y(np,np)
	parameter (nmax=1000)
	dimension indx(nmax)

	if (n.gt.nmax) then
	   write(*,*) 'increase storage in ludcmp!  n,nmax=',n,nmax
	end if
	do 12 j=1,n
	   do 11 i=1,n
	      y(i,j)=0.0d0
 11	   continue
	   y(j,j)=1.0d0
 12	continue
	call ludcmp(a,n,np,indx)
	do 13 j=1,n
	   call lubksb(a,n,np,indx,y(1,j))
 13	continue
	return
	end

	subroutine ludcmp(a,n,np,indx)
	parameter (tiny=0.0d-20,nmax=1000)
	double precision a(np,np),sum,temp
	double precision aamax,vv(nmax)
	dimension indx(n)
	if (n.gt.nmax) then
	   write(*,*) 'increase storage in ludcmp!  n,nmax=',n,nmax
	end if
	do 12 i=1,n
	   aamax=0.0d0
	   do 11 j=1,n
	      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11	   continue
	   if (aamax.eq.0.d0) then
	      print*, 'singular matrix.'
	      stop
	   endif
	   vv(i)=1./aamax
 12	continue
	do 19 j=1,n
	   if (j.gt.1) then
	      do 14 i=1,j-1
		 sum=a(i,j)
		 if (i.gt.1) then
		    do 13 k=1,i-1
		       sum=sum-a(i,k)*a(k,j)
 13		    continue
		    a(i,j)=sum
		 endif
 14	      continue
	   endif
	   aamax=0.0d0
	   do 16 i=j,n
	      sum=a(i,j)
	      if (j.gt.1) then
		 do 15 k=1,j-1
		    sum=sum-a(i,k)*a(k,j)
 15		 continue
		 a(i,j)=sum
	      endif
	      temp=vv(i)*abs(sum)
	      if (temp.ge.aamax) then
		 imax=i
		 aamax=temp
	      endif
 16	   continue
	   if (j.ne.imax) then
	      do 17 k=1,n
		 temp=a(imax,k)
		 a(imax,k)=a(j,k)
		 a(j,k)=temp
 17	      continue
	      vv(imax)=vv(j)
	   endif
	   indx(j)=imax
	   if (j.ne.n) then
!       if (a(j,j).eq.0.d0) a(j,j)=tiny
	      temp=1.0d0/a(j,j)
	      do 18 i=j+1,n
		 a(i,j)=a(i,j)*temp
 18	      continue
	   endif
 19	continue
!	if (a(n,n).eq.0.0d0) a(n,n)=tiny
	return
	end
!       
	subroutine lubksb(a,n,np,indx,b)
	double precision a(np,np),b(np),sum
	dimension indx(n)
	ii=0
	do 12 i=1,n
	   ll=indx(i)
	   sum=b(ll)
	   b(ll)=b(i)
	   if (ii.ne.0) then
	      do 11 j=ii,i-1
		 sum=sum-a(i,j)*b(j)
 11	      continue
	   else if (sum.ne.0.0d0) then
	      ii=i
	   endif
	   b(i)=sum
 12	continue
	do 14 i=n,1,-1
	   sum=b(i)
	   if (i.lt.n) then
	      do 13 j=i+1,n
		 sum=sum-a(i,j)*b(j)
 13	      continue
	   endif
	   b(i)=sum/a(i,i)
 14	continue
	return
	end
