IMPLICIT NONE
INTEGER, PARAMETER :: dbl = SELECTED_REAL_KIND(13,307)			! floating point precision (double)
INTEGER, PARAMETER :: lng = KIND(10000000)				! maximum integer value ("long")
INTEGER, allocatable :: seed1(:)
INTEGER :: seed_size1,seed_date1(8),nbins,nbinmiddle,i,np,j,npindex,k
REAL(dbl), ALLOCATABLE	:: randnomono(:),randno(:),v0R(:),Q0R(:),Q0RdR(:),v0RdR(:),Rbins(:),Q0RdRint(:),Radlist(:)
REAl(dbl) :: R0,Rstar,sigR,sigmax,vptotal,fourbythreepi,xmax,xmin,ymax,ymin,zmax,zmin,deltaR,sumvolume	

!--------------------------------------------------------------------------------------
CALL DATE_AND_TIME(VALUES=seed_date1)
CALL RANDOM_SEED(size=seed_size1)
ALLOCATE(seed1(seed_size1))
CALL RANDOM_SEED(GET=seed1)
seed1=972
CALL RANDOM_SEED(put=seed1)
DEALLOCATE(seed1)

!------------------- Monodisperse Collection ------------------------------------------
zmin=1.0_dbl+2.0_dbl
xmin=0.2_dbl*240_dbl+2.0_dbl
ymin=0.2_dbl*240_dbl+2.0_dbl
zmax=240.0_dbl - zmin
xmax=240.0_dbl - xmin
ymax=240.0_dbl - ymin
np=14_lng!250_lng!200_lng

ALLOCATE(randnomono(3_lng*np))
CALL RANDOM_NUMBER(randnomono)

R0 = 0.00263008138299_dbl ! cm

open(52,file='particle-a-14.txt')
write(52,*) np

do i=1,np
  write(52,*) i,xmin+(xmax-xmin)*randnomono(3*(i-1)+1),ymin+(ymax-ymin)*randnomono(3*(i-1)+2),zmin+(zmax-zmin)*randnomono(3*(i-1)+3),R0
end do

close(52)
write(*,*) np*(88.0_dbl/21.0_dbl)*(R0**3.0_dbl)

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

open(51,file='np250-nb20.txt')
do i=1,nbins 
  read(51,*) Rbins(i),v0R(i),v0RdR(i),Q0RdR(i),Q0RdRint(i)
end do
close(51)

np =sum(Q0RdRint)
k = 0
do i = 1,nbins
	!sumvolume = sumvolume + INT(Q0RdRint(i))*(fourbythreepi/8.0_dbl)*((Rbins(i)*0.0001_dbl)**3.0)
	do j = 1,INT(Q0RdRint(i))
		k = k + 1
		Radlist(k) = 0.5_dbl*Rbins(i)*0.0001_dbl
		!sumvolume = sumvolume + 1.0_dbl*(fourbythreepi)*((0.5_dbl*Rbins(i)*0.0001_dbl)**3.0)
		write(*,*) k,Radlist(k)
	enddo
enddo

sumvolume = 0.0_dbl

do i = 1,np
	sumvolume = sumvolume +  fourbythreepi*(Radlist(i)**3.0)
enddo 

write(*,*) k,np,sumvolume

ALLOCATE(randno(3_lng*np))
CALL RANDOM_NUMBER(randno)

open(50,file='particle-a-polydisperse.txt')
write(50,*) np

do i=1,np
 write(50,*) i,xmin+(xmax-xmin)*randno(3*(i-1)+1),ymin+(ymax-ymin)*randno(3*(i-1)+2),zmin+(zmax-zmin)*randno(3*(i-1)+3),Radlist(i)
enddo

!--------------------------------------------------------------------
close(50)

stop
end
