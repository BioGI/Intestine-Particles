IMPLICIT NONE
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(13,307)			! floating point precision (double)
INTEGER, PARAMETER :: lng = KIND(10000000)				! maximum integer value ("long")
INTEGER, allocatable :: seed1(:)
INTEGER :: seed_size1,seed_date1(8),nbins,nbinmiddle,i,np,j,npindex,k
REAL(dbl), ALLOCATABLE	:: randnomono(:),randno(:),v0R(:),Q0R(:),Q0RdR(:),v0RdR(:),Rbins(:),Q0RdRint(:),Radlist(:)
REAl(dbl) :: R0,Rstar,sigR,sigmax,vptotal,fourbythreepi,xmax,xmin,ymax,ymin,zmax,zmin,deltaR,sumvolume	
REAL(dbl) :: x_center, y_center, z_center
REAL(dbl) :: rMax, teta1Max, teta2Max, rr, teta1, teta2
REAL(dbl) :: x_particle, y_particle, z_particle
REAL(lng) :: PI, dR, Radius
!--------------------------------------------------------------------------------------
CALL DATE_AND_TIME(VALUES=seed_date1)
CALL RANDOM_SEED(size=seed_size1)
ALLOCATE(seed1(seed_size1))
CALL RANDOM_SEED(GET=seed1)
seed1=972
CALL RANDOM_SEED(put=seed1)
DEALLOCATE(seed1)

!------------------- Monodisperse Collection ------------------------------------------
x_center= 	40.0
y_center= 	40.0
z_center=	30.0

PI=3.1415926535897932384626433832

rMax= 		25.0_dbl
teta1Max= 	2*PI
teta2Max= 	2*PI


np=500_lng

ALLOCATE(randnomono(3_lng*np))
CALL RANDOM_NUMBER(randnomono)

ALLOCATE(randno(3_lng*np))
CALL RANDOM_NUMBER(randno)


R0 = 0.0001_dbl ! cm
dR = 0.0010_dbl
open(52,file='particle.dat')
write(52,*) np

do i=1,np
  rr=    (randnomono(3*(i-1)+1))**(1.0/3.0) * rMax
  teta1=  randnomono(3*(i-1)+2) 		   * teta1Max
  teta2=  randnomono(3*(i-1)+3) 	    	   * teta2Max
  Radius= randno(3*(i-1)+1)*dR + R0
 
  x_particle = x_center + rr* sin(teta1)* cos(teta2) 
  y_particle = y_center + rr* sin(teta1)* sin(teta2)
  z_particle = z_center + rr* cos(teta1)	 

  write(52,*) i, x_particle, y_particle, z_particle, Radius
end do

close(52)
!write(*,*) np*(88.0_dbl/21.0_dbl)*(R0**3.0_dbl)








!----------------- Polydisperse Collection From Yanxing -------------------------------
zmin=1.0_dbl+2.0_dbl
xmin=0.2_dbl*240_dbl+2.0_dbl
ymin=0.2_dbl*240_dbl+2.0_dbl
zmax=240.0_dbl-2_dbl
xmax=240.0_dbl - xmin
ymax=240.0_dbl - ymin
np=250_lng
fourbythreepi = 88.0_dbl/21.0_dbl
nbins = 20_lng

ALLOCATE(Rbins(nbins))
ALLOCATE(v0R(nbins))
ALLOCATE(v0RdR(nbins))
ALLOCATE(Q0RdR(nbins))
ALLOCATE(Q0RdRint(nbins))
ALLOCATE(Radlist(np))

!open(51,file='np250-nb20.txt')
!do i=1,nbins 
!  read(51,*) Rbins(i),v0R(i),v0RdR(i),Q0RdR(i),Q0RdRint(i)
!end do
!close(51)
!
!np =sum(Q0RdRint)
!k = 0
!do i = 1,nbins
!	!sumvolume = sumvolume + INT(Q0RdRint(i))*(fourbythreepi/8.0_dbl)*((Rbins(i)*0.0001_dbl)**3.0)
!	do j = 1,INT(Q0RdRint(i))
!		k = k + 1
!		Radlist(k) = 0.5_dbl*Rbins(i)*0.0001_dbl
!		!sumvolume = sumvolume + 1.0_dbl*(fourbythreepi)*((0.5_dbl*Rbins(i)*0.0001_dbl)**3.0)
!		write(*,*) k,Radlist(k)
!	enddo
!enddo
!
!sumvolume = 0.0_dbl
!
!do i = 1,np
!	sumvolume = sumvolume +  fourbythreepi*(Radlist(i)**3.0)
!enddo 
!
!write(*,*) k,np,sumvolume
!
!ALLOCATE(randno(3_lng*np))
!CALL RANDOM_NUMBER(randno)
!
!open(50,file='particle-a-polydisperse.txt')
!write(50,*) np
!
!do i=1,np
! write(50,*) i,xmin+(xmax-xmin)*randno(3*(i-1)+1),ymin+(ymax-ymin)*randno(3*(i-1)+2),zmin+(zmax-zmin)*randno(3*(i-1)+3),Radlist(i)
!enddo
!
!--------------------------------------------------------------------
!close(50)

stop
end
