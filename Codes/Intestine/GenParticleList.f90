
IMPLICIT NONE

INTEGER, PARAMETER 	:: dbl = SELECTED_REAL_KIND(13,307)			! floating point precision (double)
INTEGER, PARAMETER 	:: lng = KIND(10000000)					! maximum integer value ("long")
INTEGER   :: seed_size1,seed_date1(8),nbins,nbinmiddle,i,np,j,npindex,k

INTEGER,   ALLOCATABLE  :: seed1(:), Q0RdRint(:)
REAL(dbl), ALLOCATABLE	:: randnomono(:),randno(:)
REAL(dbl), ALLOCATABLE  :: v0R(:),Q0R(:),Q0RdR(:),v0RdR(:),Rbins(:),Radlist(:)

REAl(dbl) :: R0,Rstar,sigR,sigmax,vptotal
REAL(dbl) :: xmax,xmin,ymax,ymin,zmax,zmin,deltaR,sumvolume	
REAL(dbl) :: x_center, y_center, z_center
REAL(dbl) :: D, L, Dx, Dy, Dz
REAL(dbl) :: rMax, teta1Max, teta2Max, rr, teta1, teta2
REAL(dbl) :: x_particle, y_particle, z_particle
REAL(lng) :: PI, dR, Radius
LOGICAL   :: Falg_Couette, Flag_Measured_Particle_Dist

CALL DATE_AND_TIME(VALUES=seed_date1)
CALL RANDOM_SEED(size=seed_size1)
ALLOCATE(seed1(seed_size1))
CALL RANDOM_SEED(GET=seed1)
seed1=972
CALL RANDOM_SEED(put=seed1)
DEALLOCATE(seed1)


Flag_Measured_Particle_Dist= .TRUE.

IF (.NOT. Flag_Measured_Particle_Dist) THEN
   !---------------------------------------------------------------------------------------------------
   !----- Polydisperse Collection From Yanxing --------------------------------------------------------
   !---------------------------------------------------------------------------------------------------
   D        = 52.0
   L        = 51.0      
   Dx       = 40.0  
   Dy       = 40.0
   Dz       = 49.0

   x_center = 51.0_dbl
   y_center = 51.0_dbl
   z_center = 50.0_dbl
   PI     	 = 3.1415926535897932384626433832
   rMax   	 = 53.0_dbl
   teta1Max = 2*PI
   teta2Max = 2*PI
   nbins 	 = 20 

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
      sumvolume = sumvolume + (4.0_dbl/3.0_dbl)* PI * (Radlist(i)**3.0)
   END DO 
   write(*,*) np,sumvolume

   ALLOCATE(randno(3_lng*np))
   CALL RANDOM_NUMBER(randno)

   open(50,file='particle.dat')
   write(50,*) np
   DO i=1,np
      Falg_Couette = .TRUE.
      IF (Falg_Couette) THEN
         x_particle = ((D-Dx)/2.0_dbl) + randno(3*(i-1)+1)* Dx 
         y_particle = ((D-Dy)/2.0_dbl) + randno(3*(i-1)+2)* Dy 
         z_particle = ((L-Dz)/2.0_dbl) + randno(3*(i-1)+3)* Dz 
         write(50,*) i, x_particle, y_particle, z_particle, Radlist(i) 
      ELSE
         rr         = (randno(3*(i-1)+1))**(1.0/3.0)* rMax
         teta1      =  randno(3*(i-1)+2)           * teta1Max
         teta2      =  randno(3*(i-1)+3)           * teta2Max
         x_particle = x_center + rr* sin(teta1)* cos(teta2)
         y_particle = y_center + rr* sin(teta1)* sin(teta2)
         z_particle = z_center + rr* cos(teta1)
         write(50,*) i, x_particle, y_particle, z_particle, Radlist(i) 
      END IF
   END DO
   close(50)


ELSEIF (Flag_Measured_Particle_Dist) THEN
   Dx      = 3.95 
   Dy      = 13.9
   Dz      = 301.0
   x_center= 7.0
   y_center= 15.0
   z_center= 303.0

   open(52,file='Par_Rad_Measured.dat')
   read(52,*) np
   ALLOCATE(Radlist(np))
   DO i= 1,np 
      read(52,*) Radlist(i)
   END DO
   close(52)

   ALLOCATE(randno(3_lng*np))
   CALL RANDOM_NUMBER(randno)

   open(50,file='particle.dat')
   write(50,*) np
   DO i=1,np
      x_particle = (x_center-Dx) + randno(3*(i-1)+1)* 2.0_dbl* Dx 
      y_particle = (y_center-Dy) + randno(3*(i-1)+2)* 2.0_dbl* Dy 
      z_particle = (z_center-Dz) + randno(3*(i-1)+3)* 2.0_dbl* Dz 
      write(50,"(I5,3F14.6,E12.4)") i, x_particle, y_particle, z_particle, Radlist(i)/2000000.000000_dbl 
   END DO
   close(50)
END IF

stop
END

