IMPLICIT NONE

INTEGER, PARAMETER 	:: dbl = SELECTED_REAL_KIND(13,307)			! floating point precision (double)
INTEGER, PARAMETER 	:: lng = KIND(10000000)					! maximum integer value ("long")

INTEGER,   ALLOCATABLE  :: seed1(:), Q0RdRint(:)
REAL(dbl), ALLOCATABLE	:: randnomono(:),randno(:)
REAL(dbl), ALLOCATABLE  :: v0R(:),Q0R(:),Q0RdR(:),v0RdR(:),Rbins(:),Radlist(:)

INTEGER   :: seed_size1,seed_date1(8),nbins,nbinmiddle,i,np,j,npindex,k

REAl(dbl) :: R0,Rstar,sigR,sigmax,vptotal,fourbythreepi
REAL(dbl) :: xmax,xmin,ymax,ymin,zmax,zmin,deltaR,sumvolume	
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

!!---------------------------------------------------------------------------------------------------
!!----- Monodisperse Collection ---------------------------------------------------------------------
!!---------------------------------------------------------------------------------------------------
!x_center= 	80.0
!y_center= 	80.0
!z_center=	60.0
!PI	=	3.1415926535897932384626433832
!rMax	= 	50.0_dbl
!np	=	500_lng
!R0 	= 	0.000010_dbl 
!dR 	= 	0.000140_dbl
!teta1Max= 	2*PI
!teta2Max= 	2*PI
!
!ALLOCATE(randnomono(3_lng*np))
!CALL RANDOM_NUMBER(randnomono)
!
!ALLOCATE(randno(3_lng*np))
!CALL RANDOM_NUMBER(randno)
!
!OPEN(52,file='particle.dat')
!WRITE(52,*) np
!
!DO i=1,np
!  rr=    (randnomono(3*(i-1)+1))**(1.0/3.0)* rMax
!  teta1=  randnomono(3*(i-1)+2) 	   * teta1Max
!  teta2=  randnomono(3*(i-1)+3) 	   * teta2Max
!  Radius= randno(3*(i-1)+1)*dR + R0
! 
!  x_particle = x_center + rr* sin(teta1)* cos(teta2) 
!  y_particle = y_center + rr* sin(teta1)* sin(teta2)
!  z_particle = z_center + rr* cos(teta1)	 
!
!  write(52,*) i, x_particle, y_particle, z_particle, Radius
!END DO
!
!close(52)

!---------------------------------------------------------------------------------------------------
!----- Polydisperse Collection From Yanxing --------------------------------------------------------
!---------------------------------------------------------------------------------------------------
PI     		= 3.1415926535897932384626433832
fourbythreepi 	= 88.0_dbl/21.0_dbl
x_center	= 58.0
y_center	= 58.0
z_center	= 150.0
rMax   		= 50.0_dbl
teta1Max	= 2*PI
teta2Max	= 2*PI
nbins 		= 20 

ALLOCATE(Rbins(nbins))
ALLOCATE(v0R(nbins))
ALLOCATE(v0RdR(nbins))
ALLOCATE(Q0RdR(nbins))
ALLOCATE(Q0RdRint(nbins))

open(51,file='Par_Dist.dat')
DO i= 1,nbins 
  read(51,*) Rbins(i), v0R(i), v0RdR(i), Q0RdR(i), Q0RdRint(i)
  write(*,*) i,  Rbins(i), Q0RdRint(i)
END DO
close(51)


np= sum(Q0RdRint)

ALLOCATE(Radlist(np))

k= 0
DO i= 1, nbins
   DO j = 1, Q0RdRint(i)
      k = k + 1
      Radlist(k) = 0.5_dbl * Rbins(i)* 1.0e-6
      write(*,*) k, Radlist(k)
   END DO
END DO

sumvolume = 0.0_dbl
DO i= 1,np
   sumvolume = sumvolume +  fourbythreepi*(Radlist(i)**3.0)
END DO 
write(*,*) np,sumvolume

ALLOCATE(randno(3_lng*np))
CALL RANDOM_NUMBER(randno)

open(50,file='particle.dat')
write(50,*) np
DO i=1,np
  rr         = (randno(3*(i-1)+1))**(1.0/3.0)* rMax
  teta1      =  randno(3*(i-1)+2)           * teta1Max
  teta2      =  randno(3*(i-1)+3)           * teta2Max
  x_particle = x_center + rr* sin(teta1)* cos(teta2)
  y_particle = y_center + rr* sin(teta1)* sin(teta2)
  z_particle = z_center + rr* cos(teta1)
  write(50,*) i, x_particle, y_particle, z_particle, Radlist(i) 
END DO

!--------------------------------------------------------------------
close(50)

stop
END

