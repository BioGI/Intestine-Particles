INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(13,307)	! floating point precision (double)
INTEGER, PARAMETER :: lng = KIND(10000000)					! maximum integer value ("long")

!z0=10.0_dbl
!x0=30.0_dbl
!x1=50.0_dbl
!y0=30.0_dbl
!y1=50.0_dbl
!npx=10_lng
!npy=10_lng
!
!open(50,file='particle-a.txt')
!write(50,*) npx*npy
!  do i=1,npx
!	do j=1,npy
!  		write(50,*) (x1-x0)/(npx-1)*(i-1)+x0,(y1-y0)/(npy-1)*(j-1)+y0,z0
!	enddo
!  end do

z0=10.0_dbl
x0=36.0_dbl
x1=46.0_dbl
y0=41.0_dbl
np=10_lng
open(50,file='particle-a.txt')
write(50,*) np
do i=1,np
  write(50,*) (x1-x0)/(np-1)*(i-1)+x0,y0,z0
end do


!x1=40.5
!x2=80.5
!y0=43.0
!y1=58.0
!y2=78.0
!np=10
!open(50,file='particle-a.dat')
!write(50,*) np
!write(50,*) 2*np
!do i=1,np
!  write(50,*) x1,(y1-y0)/(np-1)*(i-1)+y0,0.0
!end do
!do i=1,np
!  write(50,*) x2,(y2-y0)/(np-1)*(i-1)+y0,0.0
!end do
close(50)

stop
end
