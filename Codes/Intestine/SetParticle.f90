INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(13,307)		! floating point precision (double)
INTEGER, PARAMETER :: lng = KIND(10000000)			! maximum integer value ("long")


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

close(50)

stop
end
