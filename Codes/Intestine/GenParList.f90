!===================================================================================================
! This subroutine randomly creates the particle's locations  
!===================================================================================================
USE SetPrecision
USE Setup
USE Geometry

IMPLICIT NONE

INTEGER,   ALLOCATABLE  :: seed(:), Q0RdRint(:)
REAL(dbl), ALLOCATABLE	:: randno(:), Rlist(:)
REAL(dbl) :: D_Couette, L_Couette, Dx_Couette, Dy_Couette, Dz_Couette
REAL(dbl) :: x_particle, y_particle, z_particle
REAL(dbl) :: R_Particle, R_Par_Max, D_Par_Max, R_Boundary
REAL(dbl) :: R_left, R_right, dz
INTEGER   :: n,num_Par, i, j, k, z_left, z_right, clock 
LOGICAL   :: Falg_Couette

CALL ReadInput
nuL  = (2.0_dbl*tau -1.0_dbl)/6.0_dbl	! lattice kinematic viscosity
denL = 1.0_dbl		                	  ! arbitrary lattice density (1.0 for convenience)
kMin = 1
kMax = nz
nzSub= nz

CALL AllocateArrays                   ! allocate global variable arrays  
CALL Geometry_Setup
CALL BoundaryPosition

CALL RANDOM_SEED(size=n)
ALLOCATE(seed(n))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37*(/(i-1, i=1, n)/)
CALL RANDOM_SEED(put=seed)
DEALLOCATE(seed)
ALLOCATE(randno(10000000))
CALL RANDOM_NUMBER(randno)

!---------------------------------------------------------------------------------------------------
!----- Polydisperse Collection From Yanxing --------------------------------------------------------
!---------------------------------------------------------------------------------------------------

!-----  For Couette Cases
D_Couette    = 52.0
L_Couette    = 51.0      
Dx_Couette   = 40.0  
Dy_Couette   = 40.0
Dz_Couette   = 49.0
Falg_Couette = .FALSE.

R_Par_Max=0
open(50,file='Particle_Sizes.dat')
read(50,*) num_Par
ALLOCATE(Rlist(num_Par))
DO i= 1, num_Par 
   read(50,*) Rlist(i)        
   Rlist(i)= Rlist(i)/1.0e6         ! Changing units from micron to meter
   IF (Rlist(i).GT.R_Par_Max) THEN
      R_Par_Max= Rlist(i)
   ENDIF
END DO
close(50)


D_Par_Max = 2.0*R_Par_Max
write(*,*) 'D_Par_Max',D_Par_Max


open(51,file='particle.dat')
write(51,*) num_Par
j=0
DO i=1,num_Par
   IF (Falg_Couette) THEN
      x_particle = ((D-Dx_Couette)/2.0_dbl) + randno(3*(i-1)+1)* Dx_Couette 
      y_particle = ((D-Dy_Couette)/2.0_dbl) + randno(3*(i-1)+2)* Dy_Couette 
      z_particle = ((L-Dz_Couette)/2.0_dbl) + randno(3*(i-1)+3)* Dz_Couette 
      write(51,*) i, x_particle, y_particle, z_particle, Rlist(i) 
   ELSE
100  j=j+1
     x_particle = 1.0 + randno(3*(j-1)+1) * (nx-1) 
     y_particle = 1.0 + randno(3*(j-1)+2) * (ny-1) 
     z_particle = 1.0 + randno(3*(j-1)+3) * (nz-1)
     R_particle = SQRT( (x_particle-((nx+1)/2.0))**2.0 + (y_particle-((ny+1)/2.0))**2.0 ) 
     R_particle = R_Particle* xcf
     z_left     = Floor(z_particle)
     z_right    = CEILING(z_particle)
     R_left     = r(z_left) 
     R_right    = r(z_right)
     dz         = z_right-z_left
     IF (dz.NE. 0) THEN
        R_Boundary = R_left + ((z_particle-z_left)/(z_right-z_left)) * (R_right-R_left)  
     ELSE
        R_Boundary = R_left
     END IF

     write(*,*) i,x_particle,y_particle,z_particle,z_left,z_right,R_left,R_right,R_Boundary,R_particle,D_Par_Max

     IF (R_particle .LT. (R_Boundary-2.0*D_Par_Max) ) THEN   
        write(51,"(I6,3F16.8,E14.6)") i, x_particle, y_particle, z_particle, Rlist(i) 
     ELSE
        GOTO 100
     END IF
   END IF
END DO
close(51)

stop
END

