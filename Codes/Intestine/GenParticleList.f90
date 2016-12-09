
IMPLICIT NONE

INTEGER, PARAMETER 	:: dbl = SELECTED_REAL_KIND(13,307)			! floating point precision (double)
INTEGER, PARAMETER 	:: lng = KIND(10000000)					! maximum integer value ("long")
INTEGER   :: seed_size1,seed_date1(8),nbins,nbinmiddle,i,np,j,npindex,k

INTEGER,   ALLOCATABLE  :: seed1(:), Q0RdRint(:)
REAL(dbl), ALLOCATABLE	:: randnomono(:),randno(:)
REAL(dbl), ALLOCATABLE  :: Rlist(:)

REAl(dbl) :: R0,Rstar
REAL(dbl) :: xmax,xmin,ymax,ymin,zmax,zmin	
REAL(dbl) :: x_center, y_center, z_center
REAL(dbl) :: D, L, Dx, Dy, Dz
REAL(dbl) :: rMax, teta1Max, teta2Max, rr, teta1, teta2
REAL(dbl) :: x_particle, y_particle, z_particle
REAL(lng) :: PI 
LOGICAL   :: Falg_Couette

CALL DATE_AND_TIME(VALUES=seed_date1)
CALL RANDOM_SEED(size=seed_size1)
ALLOCATE(seed1(seed_size1))
CALL RANDOM_SEED(GET=seed1)
seed1=972
CALL RANDOM_SEED(put=seed1)
DEALLOCATE(seed1)

!---------------------------------------------------------------------------------------------------
!----- Polydisperse Collection From Yanxing --------------------------------------------------------
!---------------------------------------------------------------------------------------------------

!-----  For Couette Cases
D        = 52.0
L        = 51.0      
Dx       = 40.0  
Dy       = 40.0
Dz       = 49.0

!-----  For IntestineCases

x_center = 51.0_dbl
y_center = 51.0_dbl
z_center = 50.0_dbl
PI     	 = 3.1415926535897932384626433832
rMax   	 = 53.0_dbl
teta1Max = 2*PI
teta2Max = 2*PI

open(50,file='Particle_Sizes.dat')
read(50,*) np
ALLOCATE(Rlist(np))
DO i= 1, np 
   read(50,*) Rlist(i)
   WRITE(*,*) i,Rlist(i)
END DO
close(50)

ALLOCATE(randno(3_lng*np))
CALL RANDOM_NUMBER(randno)

open(51,file='particle.dat')
write(51,*) np
DO i=1,np
   Falg_Couette = .TRUE.
   IF (Falg_Couette) THEN
      x_particle = ((D-Dx)/2.0_dbl) + randno(3*(i-1)+1)* Dx 
      y_particle = ((D-Dy)/2.0_dbl) + randno(3*(i-1)+2)* Dy 
      z_particle = ((L-Dz)/2.0_dbl) + randno(3*(i-1)+3)* Dz 
      write(51,*) i, x_particle, y_particle, z_particle, Rlist(i) 
   ELSE
     rr         = (randno(3*(i-1)+1))**(1.0/3.0)* rMax
     teta1      =  randno(3*(i-1)+2)           * teta1Max
     teta2      =  randno(3*(i-1)+3)           * teta2Max
     x_particle = x_center + rr* sin(teta1)* cos(teta2)
     y_particle = y_center + rr* sin(teta1)* sin(teta2)
     z_particle = z_center + rr* cos(teta1)
     write(51,"(I6,3F16.8,E14.6)") i, x_particle, y_particle, z_particle, Rlist(i) 
   END IF
END DO
close(51)

stop
END

