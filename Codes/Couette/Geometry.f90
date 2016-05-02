!==================================================================================================
MODULE Geometry	! Defines the geometry for the simulation
						! Subroutines (NodeFlags, BoundaryVelocity)
!==================================================================================================
USE SetPrecision      
USE Setup
USE LBM
USE IC
USE BClbm
USE BCscalar
USE MPI

IMPLICIT NONE 

CONTAINS


!--------------------------------------------------------------------------------------------------
SUBROUTINE Geometry_Setup				! sets up the geometry
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER :: isize,idate(8)				! size of seed array for the random number genreator, array for output of DATE_AND_TIME
INTEGER,ALLOCATABLE  :: iseed(:)			! seeds for random number generator
INTEGER(lng) :: i,j,k,kk,iCon,it,iPer,nPers_INT		! index variables
INTEGER(lng) :: nvz,nvt,n,g								! index variables
INTEGER(lng) :: mpierr										! MPI standard error variable 
REAL(dbl) :: macroFreq										! macroscopic contraction frequency
INTEGER(lng) :: xaxis,yaxis								! axes index variables

REAL(dbl) :: D_X, D_Y
D_X = 20.000_dbl * D 
D_Y= 0.50_dbl *D

! Define the lattice <=> physical conversion factors
IF(domaintype .EQ. 0) THEN
        xcf 		= (0.5_lng*D_x)/(nx-1_lng)		! length conversion factor: x-direction
        ycf 		= (0.5_lng*D_y)/(ny-1_lng)		! length conversion factor: y-direction
ELSE
        ! begin Balaji added
        xcf 		= (1.0_lng*D_x)/(nx-1_lng)		! length conversion factor: x-direction
        ycf 		= (1.0_lng*D_y)/(ny-1_lng)		! length conversion factor: y-direction
        ! end Balaji added
ENDIF

zcf 		= L/nz					! length conversion factor: z-direction
tcf 		= nuL*((xcf*xcf)/nu)			! time conversion factor
dcf 		= den/denL				! density conversion factor
vcf 		= xcf/tcf				! velocity conversion factor
pcf 		= cs*cs*vcf*vcf				! pressure conversion factor

! Determine the number of time steps to run
nt = ANINT((nPers*Tmix)/tcf)

! Initialize arrays
node		= -99_lng				! node flag array
rDomIn		= 0.0_dbl				! radius at each z-location
rIn		= 0.0_dbl				! temporary radius array for entire computational domain
rDomOut		= 0.0_dbl				! radius at each z-location
rOut		= 0.0_dbl				! temporary radius array for entire computational domain
velDomIn	= 0.0_dbl				! wall velocity at each z-location (global)
velIn		= 0.0_dbl				! wall velocity at each z-location (local)
velDomOut	= 0.0_dbl				! wall velocity at each z-location (global)
velOut		= 0.0_dbl				! wall velocity at each z-location (local)

! Check to ensure xcf=ycf=zcf (LBM grid must be cubic)
IF((ABS(xcf-ycf) .GE. 1E-8) .OR. (ABS(xcf-zcf) .GE. 1E-8) .OR. (ABS(ycf-zcf) .GE. 1E-8)) THEN
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) "Conversion factors not equal... Geometry_Setup.f90: Line 93."
  WRITE(1000,*) "xcf=", xcf, "ycf=", ycf, "zcf=", zcf
  WRITE(1000,*) "L=", L, "D=", D
  WRITE(1000,*) "nx=", nx, "ny=", ny, "nz=", nz
  CLOSE(1000)
  STOP
END IF


! IF CONDITION TO CHECK IF THE DOMAIN TO BE MODELLED IS FULL CYLINDER OR JUST A QUARTER OF A CYLINDER
IF(domaintype .EQ. 0) THEN 
      ! Fill out x,y,z arrays (local)
      DO i=0,nxSub+1
        x(i) = ((iMin - 1_lng) + (i-1_lng))*xcf
      END DO
      
      DO j=0,nySub+1
        y(j) = ((jMin - 1_lng) + (j-1_lng))*ycf
      END DO
      
      DO k=0,nzSub+1
        z(k) = (((kMin - 1_lng) + k) - 0.5_dbl)*zcf
      END DO
      
      ! Fill out xx,yy,zz arrays (global)
      DO i=0,nx+1
        xx(i) = (i-1_lng)*xcf
      END DO
      
      DO j=0,ny+1
        yy(j) = (j-1_lng)*ycf
      END DO
      
      DO k=0,nz+1
        zz(k) = (k - 0.5_dbl)*zcf
      END DO
      
      ! Center node locations
      Ci = 1	
      Cj = 1
      Ck = ANINT(0.5_dbl*nz)

ELSE
      xaxis=ANINT(0.5_dbl*(nx+1))
      yaxis=ANINT(0.5_dbl*(ny+1))
      
      ! Fill out x,y,z arrays (local)
      DO i=0,nxSub+1
        x(i) = ((iMin - 1_lng - (xaxis-1_lng)) + (i-1_lng))*xcf
      END DO
      
      DO j=0,nySub+1
        y(j) = ((jMin - 1_lng - (yaxis-1_lng)) + (j-1_lng))*ycf
      END DO
      
      DO k=0,nzSub+1
        z(k) = (((kMin - 1_lng) + k) - 0.5_dbl)*zcf
      END DO
      
      ! Fill out xx,yy,zz arrays (global)
      DO i=0,nx+1
        xx(i) = (i-1_lng-(xaxis-1_lng))*xcf
      END DO
      
      DO j=0,ny+1
        yy(j) = (j-1_lng-(yaxis-1_lng))*ycf
      END DO
      
      DO k=0,nz+1
        zz(k) = (k - 0.5_dbl)*zcf
      END DO
      
      ! Center node locations
      Ci = xaxis
      Cj = yaxis
      Ck = ANINT(0.5_dbl*nz)
ENDIF	

! Mode 1 - Peristalsis
a1		= (0.5_dbl*D)/(2.0_dbl - epsOVERa1)				! mean half-width of wave1
eps1 		= epsOVERa1*a1							! occlusional distance
lambda1		= L/numw1							! wavelength
aOVERlam1	= a1/lambda1							! ratio of mean half-width to wavelength 
kw1		= (2.0_dbl*PI)/lambda1						! wave number
amp1		= 0.5_dbl*((0.5_dbl*D)-eps1)					! amplitude of the wave
Tp		= 1.0_dbl 		! lambda1/s1							! peristaltic period
Re1		= ((s1*(0.5_dbl*D))/nu)*((0.5_dbl*D)/lambda1)			! Reynolds number based on mode 1

! Mode 2 - Segmental Contractions
a2		= (0.5_dbl*D)/(2.0_dbl - epsOVERa2)				! mean half-width of wave1 (based on peristalsis definition)
eps2 		= epsOVERa2*a2							! occlusional distance
lambda2		= L/numw2							! wavelength (physical units)
nlambda2	= nz/numw2							! wavelength (nodes)
aOVERlam2	= a2/lambda2							! ratio of mean half-width to wavelength 
amp2		= 0.5_dbl*((0.5_dbl*D)-eps2)					! amplitude of the wave
shift2		= 0.5_dbl*((0.5_dbl*D)+eps2)					! amplitude of the wave
segment		= nlambda2/6_lng						! length of each segment of the segmental wave   !!!!! CAREFUL HERE WITH SYMMETRY!
seg1L		= 1_lng + segment						! left point of sloped segement 1
seg1R		= 1_lng + 2_lng*segment						! right point of sloped segement 1
seg2R		= nlambda2 - segment						! right point of sloped segement 2
seg2L		= nlambda2 - (2_lng*segment)					! left point of sloped segement 2
s2		= (0.5_dbl*D)/Ts						! speed of collapse fo segmental contraction
Re2		= (s2*(0.5_dbl*D))/nu						! Reynolds number based on mode 2

! Allocate and initialize the villi arrays
numVilli				= numVilliZ*numVilliTheta		! determine the total number of villi
IF(numVilliGroups .GT. 1_lng) THEN
  numVilliActual	= numVilli - numVilliGroups*numVilliTheta		! determine the actual number of villi (subracting those skipping in grouping)
ELSE
  numVilliActual = numVilli							! if there is only 1 group, numVilli is numVilliActual
END IF
ALLOCATE(villiLoc(numVilli,10))							! location and other information of the villi
villiLoc = 0.0_dbl

! Set up the villi groups
ALLOCATE(villiGroup(numVilli))							! array of which group the villi are in
villiGroup = 0.0_dbl
n = 0_lng
DO nvz=1,numVilliZ
  DO nvt=1,numVilliTheta

    n = n + 1_lng								! villus number

    DO g=1,numVilliGroups

      IF((n .GT. ((g-1_lng)*numVilliTheta*(numVilliZ/numVilliGroups))) .AND.		&
         (n .LE. (g*numVilliTheta*(numVilliZ/numVilliGroups)))) THEN
        villiGroup(n) = g
      END IF

    END DO

  END DO
END DO

! Convert villous length and radius of the villi to meters
Lv = Lv*(0.000001_dbl)
Rv = Rv*(0.000001_dbl)	

! Determine villous frequency
macroFreq = 1.0_dbl/Tmix
vFreqT = freqRatioT*macroFreq
vFreqZ = freqRatioZ*macroFreq

IF(freqRatioT .LT. 0.00000001_dbl) THEN						! if the frequency ratio = 0 (no active villous motion) then zero out the contribution
  activeVflagT = 0.0_dbl
ELSE
  activeVflagT = 1.0_dbl
END IF

IF(freqRatioZ .LT. 0.00000001_dbl) THEN						! if the frequency ratio = 0 (no active villous motion) then zero out the contribution
  activeVflagZ = 0.0_dbl
ELSE
  activeVflagZ = 1.0_dbl
END IF

! Convert villiAngle from degrees to radians
villiAngle = (villiAngle/180.0_dbl)*PI
IF(restart .EQV. .FALSE.) THEN
  IF(randORord .EQ. RANDOM) THEN
    ALLOCATE(rnd(2_lng*numVilli))						! allocate the array of random numbers for random villi phase angles
    IF(myid .EQ. master) THEN
      CALL DATE_AND_TIME(VALUES=idate)						! get the date and time (for more different seeds each time)
      CALL RANDOM_SEED(SIZE=isize)						! get the size of the seed array
      ALLOCATE(iseed(isize))							! allocate the seed array
      CALL RANDOM_SEED(GET=iseed)						! get the seed array
      iseed = iseed*(idate(8)-500_lng)    					! idate(8) contains millisecond
      CALL RANDOM_SEED(PUT=iseed)						! use the seed array
      CALL RANDOM_NUMBER(rnd)							! get the actual random numbers
      DEALLOCATE(iseed)
      ! print the rnd array for restarting
      OPEN(1777,FILE='rnd.dat')
      DO i=1,2_lng*numVilli
        WRITE(1777,*) rnd(i)
      END DO
      CLOSE(1777)
    END IF
    CALL MPI_BCAST(rnd,2_lng*numVilli,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)! send/recv rnd on all processing units
  END IF

  CALL AdvanceGeometry

END IF

!------------------------------------------------
END SUBROUTINE Geometry_Setup
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE AdvanceGeometry	! advances the geometry in time
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

CALL BoundaryPosition		! Calculate the radius at the current time step
CALL VilliPosition
CALL BoundaryVelocity		! Calculate the velocity at boundary point
CALL SetNodes			! Flag the fluid/solid nodes based on the new geometry
!------------------------------------------------
END SUBROUTINE AdvanceGeometry
!------------------------------------------------











!--------------------------------------------------------------------------------------------------
SUBROUTINE BoundaryPosition		! Calculates the position of the wall at the current time step
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

REAL(dbl) :: h1(0:nz+1)				! Mode 1 (peristalsis)
REAL(dbl) :: h1In(0:nz+1),h1Out(0:nz+1)		! Mode 1 (Couette device)
REAL(dbl) :: h2(0:nz+1)				! Mode 2	(segmental)
REAL(dbl) :: Ac, lambdaC, shiftC		! temporary variables for the cos slopes
REAL(dbl) :: time				! time
INTEGER(lng) :: i,j,ii,k			! indices
REAL(dbl) :: D_X, D_Y

!----- Initialize Variables
time 	= 0.0_dbl				! time					
h1In 	= 0.0_dbl				! mode 1 height
h1Out 	= 0.0_dbl				! mode 2 height
rDomIn	= 0.0_dbl				! summed height
rDomOut	= 0.0_dbl				! summed height

D_X = 20*D
D_Y= 0.50_dbl *D

time	= iter*tcf
DO i=0,nz-1
   h1Out(i) = -0.38_dbl*D_x + 5.0000e-5 + s1*time 	   
   h1In(i)  = -0.48_dbl*D_x + 5.0000e-5 +s1*time 	 
END DO

!----- since PI cannot be stored exactly, the wavelength(s) does/do not EXACTLY span the domain...
!----- set h1(nz) to h1(0) and h1(nz+1) to h(1) to ensure periodicity
h1In(nz)   = h1In(0)
h1In(nz+1) = h1In(1)
h1Out(nz)  = h1Out(0)
h1Out(nz+1)= h1Out(1)

DO i=0,nz+1
  rDomIn(i) = h1In(i)
  rDomOut(i)= h1Out(i)
END DO

rIn(0:nzSub+1) = rDomIn(kMin-1:kMax+1)
rOut(0:nzSub+1)= rDomOut(kMin-1:kMax+1)

IF (myid .EQ. master) THEN
   CALL SurfaceArea
END IF

!------------------------------------------------
END SUBROUTINE BoundaryPosition
!------------------------------------------------









!--------------------------------------------------------------------------------------------------
SUBROUTINE BoundaryVelocity			! defines the velocity of the solid boundaries (fills "ub", "vb", and "wb" arrays)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

REAL(dbl) :: v1(0:nz+1), v2(0:nz+1)		! velocity arrays for each mode
REAL(dbl) :: v1In(0:nz+1), v1Out(0:nz+1)	! velocity arrays for each mode
REAL(dbl) :: lambdaC				! wavelength of the cos segments (mode 2)
REAL(dbl) :: time				! time
INTEGER(lng) :: i,j,ii				! indices

!----- Initialize Variables
time	  = 0.0_dbl				! time
velDomIn  = 0.0_dbl				! summed velocity
velDomOut = 0.0_dbl				! summed velocity
v1In	  = 0.0_dbl				! mode 1 velocity
v1Out	  = 0.0_dbl				! mode 1 velocity

time = iter*tcf

DO i=0,nz-1  					! Balaji added to ensure periodicity just like in h1. 
   v1In(i) = s1 	! s1* 0.5_dbl   
   v1Out(i)= s1 	!-s1* 0.5_dbl
END DO

v1In(nz)=   v1In(0)
v1In(nz+1)= v1In(1)
v1Out(nz)=  v1Out(0)
v1Out(nz+1)=v1Out(1)

!----- Sum the modes in a weighted linear combination
DO i=0,nz+1
   velDomIn(i) = v1In(i)
   velDomOut(i) = v1Out(i)
END DO

velIn(0:nzSub+1) = velDomIn(kMin-1:kMax+1)/vcf
velOut(0:nzSub+1) = velDomOut(kMin-1:kMax+1)/vcf

!------------------------------------------------
END SUBROUTINE BoundaryVelocity
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE SetNodes					! defines the geometry via the "node" array of flags
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 
INTEGER(lng)	:: i,j,k,m,iComm	! index variables
REAL(dbl)		:: rijk				! radius of the current node
REAL(dbl)      :: ubx,uby,ubz		! boundary velocity
INTEGER(lng) :: mpierr										! MPI standard error variable 
INTEGER(lng) :: numFluid
REAL(dbl) :: phiInTemp,phiOutTemp,phiTotalOld,phiTotalNew,rhoInTemp,rhoOutTemp,zcf3

zcf3 = zcf*zcf*zcf 

phiInTemp   = 0.0_dbl
phiOutTemp  = 0.0_dbl
phiTotalOld = 0.0_dbl
phiTotalNew = 0.0_dbl
rhoInTemp   = 0.0_dbl
rhoOutTemp  = 0.0_dbl

DO k=1,nzSub
  DO j=1,nySub
    DO i=1,nxSub
       IF (node(i,j,k) .EQ. FLUID) THEN
          phiTotalOld = phiTotalOld + phi(i,j,k)*zcf3
       END IF
    END DO
  END DO
END DO

!----- Flag the interior nodes and give values to nodes that just came in
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         rijk = x(i)							! height at current location
         IF ((rijk .LT. rOut(k)).AND.(rijk .GT. rIn(k))) THEN
            IF (node(i,j,k) .EQ. SOLID) THEN				! just came into the domain from interior solid
               ubx = velIn(k) 	! 0.0_dbl
               uby = 0.0_dbl
               ubz = 0.0_dbl 	! velIn(k)
               CALL SetProperties(i,j,k,ubx,uby,ubz)
               phiInTemp = phiInTemp + phi(i,j,k)*zcf3
	       rhoInTemp = rhoInTemp + rho(i,j,k)*zcf3
            ELSE IF(node(i,j,k) .EQ. SOLID2) THEN 			! Just came into the domain from exterior solid
               ubx = velOut(k) 	!0.0_dbl
               uby = 0.0_dbl
               ubz = 0.0_dbl 	! velOut(k)
               CALL SetProperties(i,j,k,ubx,uby,ubz)
               phiInTemp = phiInTemp + phi(i,j,k)*zcf3
	       rhoInTemp = rhoInTemp + rho(i,j,k)*zcf3
            END IF
	
            node(i,j,k)	= FLUID						! reset the SOLID node that just came in to FLUID
         ELSE
	    IF (node(i,j,k) .EQ. FLUID) THEN
	       phiOutTemp = phiOutTemp + phi(i,j,k)*zcf3
               rhoOutTemp = rhoOutTemp + rho(i,j,k)*zcf3
	    END IF

            IF (rijk .GE. rOut(k)) THEN
	        node(i,j,k) = SOLID2					! if rijk is GE rOut(k) then it's a SOLID2 node (Exterior)
	    ELSE IF (rijk .LE. rIn(k)) THEN		
        	node(i,j,k) = SOLID					! if rijk is LE rIn(k) then it's a SOLID node (Interior)
	    ENDIF
		
         END IF
       END DO
   END DO
END DO

numFluid = 0_lng
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
	 IF (node(i,j,k) .EQ. FLUID) THEN
	    numFluid = numFluid + 1_lng
	    phiTotalNew =phiTotalNew + phi(i,j,k)*zcf3
	 END IF
      END DO
   END DO
END DO

phiInNodes = phiInNodes + phiInTemp
phiOutNodes = phiOutNodes + phiOutTemp

DO iComm=1,2
   i = YZ_RecvIndex(OppCommDir(iComm))					! i index of the phantom nodes
   DO j=0,nySub+1_lng
      rijk = x(i)							! height at current location
      DO k=0,nzSub+1_lng
         IF ((rijk .LT. rOut(k)).AND.(rijk .GT. rIn(k))) THEN
            node(i,j,k) = FLUID						! set the SOLID node that just came in to FLUID
         ELSE
	    IF (rijk .GE. rOut(k)) THEN
	       node(i,j,k) = SOLID2					! if rijk is GE rOut(k) then it's a SOLID2 node (Exterior)
	    ELSE IF (rijk .LE. rIn(k)) THEN
               node(i,j,k) = SOLID					! if rijk is LE rIn(k) then it's a SOLID node (Interior)
	    ENDIF
         END IF
     END DO
  END DO
END DO


!----- ZX Faces
DO iComm=3,4
   j = ZX_RecvIndex(OppCommDir(iComm))					! j index of the phantom nodes
   DO i=0,nxSub+1_lng
      rijk = x(i)							! height at current location
      DO k=0,nzSub+1_lng
         IF ((rijk .LT. rOut(k)).AND.(rijk .GT. rIn(k))) THEN
            node(i,j,k) = FLUID						! set the SOLID node that just came in to FLUID
         ELSE
	    IF (rijk .GE. rOut(k)) THEN
	       node(i,j,k) = SOLID2					! if rijk is GE rOut(k) then it's a SOLID2 node (Exterior)
	    ELSE IF (rijk .LE. rIn(k)) THEN
               node(i,j,k) = SOLID					! if rijk is LE rIn(k) then it's a SOLID node (Interior)
	    ENDIF
         END IF
      END DO
   END DO
END DO


!----- XY Faces
DO iComm=5,6
   k = XY_RecvIndex(OppCommDir(iComm))					! k index of the phantom nodes
   DO j=0,nySub+1_lng
      DO i=0,nxSub+1_lng
         rijk = x(i)							! height at current location
         IF ((rijk .LT. rOut(k)).AND.(rijk .GT. rIn(k))) THEN
            node(i,j,k) = FLUID						! set the SOLID node that just came in to FLUID
         ELSE
	    IF (rijk .GE. rOut(k)) THEN
	       node(i,j,k) = SOLID2					! if rijk is GE rOut(k) then it's a SOLID2 node (Exterior)
	    ELSE IF (rijk .LE. rIn(k)) THEN
               node(i,j,k) = SOLID					! if rijk is LE rIn(k) then it's a SOLID node (Interior)
	    ENDIF
          END IF
    END DO
  END DO
END DO


CALL SetNodesVilli							! set the villi node flags

IF (domaintype .EQ. 0) THEN  						! only needed when planes of symmetry exist
   CALL SymmetryBC							! ensure symmetric node placement
ENDIF

!------------------------------------------------
END SUBROUTINE SetNodes
!------------------------------------------------








!===================================================================================================
SUBROUTINE SetProperties(i,j,k,ubx,uby,ubz)    
!===================================================================================================
! Set properties to nodes that just came into the fluid domain (uncovered)
IMPLICIT NONE

INTEGER(lng),INTENT(IN) :: i,j,k       			! current node location
REAL(dbl)   ,INTENT(IN) :: ubx,uby,ubz 			! velocity of the boundary
INTEGER(lng)    :: m 		          		! Direction index
INTEGER(lng)    :: ip1,jp1,kp1,iB,jB,kB  		! First neighboring node location
INTEGER(lng)    :: ip2,jp2,kp2,iC,jC,kC     	 	! Second neighboring node location
CHARACTER(7)    :: iter_char            		! iteration stored as a character
REAL(dbl)       :: feq                  		! equilibrium distribution function
REAL(dbl)	:: Geom_norm_x,Geom_norm_y,Geom_norm_z	! geometry normal vector
REAL(dbl)	:: q,n_prod,n_prod_max

!----- Enforcing velocity and denisty for the uncovered nodes ignoring the averaged value ----------
rho(i,j,k)= denL
u(i,j,k)  = ubx                                         
v(i,j,k)  = uby
w(i,j,k)  = ubz

!----- Estimating phi for uncovered node based on the prescribed BC --------------------------------
Geom_norm_x= 1.0
Geom_norm_y= 0.0
Geom_norm_z= 0.0

n_prod_max= 0.0_dbl

!---------------------------------------------------------------------------------------------------
!----- Estimating  phi for uncovered node in case of Dirichlet BC ----------------------------------
!----- One Fluid neighboring node is needed for interpolation --------------------------------------
! ---- Uncoverd node is A and neighbor is B --------------------------------------------------------
!---------------------------------------------------------------------------------------------------
IF (coeffGrad .eq. 0) then
   !----- Finding the mth direction closest to normal vector
   !----- Only one fluid neighboring node is needed for interpolation
   DO m=0,NumDistDirs
      ip1= i+ ex(m) 		
      jp1= j+ ey(m)			
      kp1= k+ ez(m)
      IF (node(ip1,jp1,kp1) .EQ. FLUID) THEN
         n_prod= abs( Geom_norm_x * (ip1-i) + Geom_norm_y * (jp1-j) + Geom_norm_z * (kp1-k))	
         IF (n_prod_max .LT. n_prod)THEN
            n_prod_max= n_prod
            iB= ip1
            jB= jp1
            kB= kp1
         END IF
      END IF
   END DO
   CALL qCalcFarhad(i,q)
   phi(i,j,k)= (phi(iB,jB,kB)-phiWall)*q/(1.0_dbl+q) + phiWall
END IF

!---------------------------------------------------------------------------------------------------
!----- Estimating  phi for the uncovered node in case of Neumann BC --------------------------------
!----- Two Fluid neighboring nodes are needed for extrapolation ------------------------------------
!----- Uncovered node is A, first neighbor is B and second neighbor is C ---------------------------
!---------------------------------------------------------------------------------------------------
IF (coeffGrad .ne. 0) then
   DO m=0,NumDistDirs
      !----- Finding the two neighboring nodes in mth direction -----------------------------------
      ip1= i+ ex(m)  
      jp1= j+ ey(m)  
      kp1= k+ ez(m) 
      ip2= i+ 2* ex(m)
      jp2= j+ 2* ey(m)
      kp2= k+ 2* ez(m)

      !------ Taking care of y-dir and z-dir periodicity for two neighboring nodes -----------------
      IF (jp1 .GT. ny) THEN
         jp1= jp1- ny
      ELSE IF (jp1 .LT. 1) THEN
         jp1= jp1+ny
      END IF
      IF (jp2 .GT. ny) THEN
         jp2= jp2- ny
      ELSE IF (jp2 .LT. 1) THEN
         jp2= jp2+ ny
      END IF
      IF (kp1 .GT. nz) THEN
         kp1= kp1- nz
      ELSE IF (kp1 .LT. 1) THEN
         kp1= kp1+ nz
      END IF
      IF (kp2 .GT. nz) THEN
         kp2= kp2- nz
      ELSE IF (kp2 .LT. 1) THEN
         kp2= kp2+ nz
      END IF

      !----- Finding the mth direction vetor closest to normal vector ------------------------------
      IF (node(ip2,jp1,kp1) .EQ. FLUID) THEN
         IF (node(ip2,jp2,kp2) .EQ. FLUID) THEN
            n_prod= abs( Geom_norm_x*(ip1-i) + Geom_norm_y*(jp1-j) + Geom_norm_z*(kp1-k))
            IF (n_prod_max .LT. n_prod)THEN
               n_prod_max= n_prod
               iB= ip1
               jB= jp1
               kB= kp1
               iC= ip2
               jC= jp2
               kC= kp2
            END IF
         END IF
      END IF  
   END DO
   
   q= 1.0_dbl	! approximating that uncovered node is at the wall (it is close enough ...)
   phi(i,j,k)= ( (phi(iB,jB,kB)*(1.0+q)*(1.0+q)/(1.0+2.0*q))	&
 	       - (phi(iC,jC,kC)*q*q/(1.0+2.0*q)) 		&
               - (q*(1+q)/(1+2.0*q))*(coeffConst/coeffGrad) ) 	&
	       / (1.0 - (q*(1+q)/(1+2.0*q))*(coeffPhi/coeffGrad))  
END IF

!---------------------------------------------------------------------------------------------------
!----- distribution functions set to equilibrium at the uncovered node -----------------------------
!---------------------------------------------------------------------------------------------------
DO m=0,NumDistDirs
  CALL Equilibrium_LOCAL(m,rho(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),feq)                   
  f(m,i,j,k) = feq
END DO

!===================================================================================================
END SUBROUTINE SetProperties
!===================================================================================================










!--------------------------------------------------------------------------------------------------
SUBROUTINE SurfaceArea			! calculate the surface area at the current time and write it to a file
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

!REAL(dbl) :: SA					! surface area
!REAL(dbl) :: r2,r1,z2,z1		! radius and z-location for each set of consecutive points
!INTEGER(lng) :: kk				! index variable
REAL(dbl) :: SA,SAIn,SAOut					! surface area
REAL(dbl) :: r2,r1,z2,z1		! radius and z-location for each set of consecutive points
INTEGER(lng) :: kk				! index variable

SAIn = 0.0_dbl
SAIn = (zz(1)-zz(0))*(yy(ny)-yy(1))
SAIn = SAIn + (zz(nz+1)-zz(nz))*(yy(ny)-yy(1))

!----- interior domain nodes
DO kk=1,nz-1
  SAIn = SAIn + (zz(kk+1)-zz(kk))*(yy(ny)-yy(1))
END DO

SAOut = 0.0_dbl
SAOut = (zz(1)-zz(0))*(yy(ny)-yy(1))
SAOut = SAOut + (zz(nz+1)-zz(nz))*(yy(ny)-yy(1))

!----- interior domain nodes
DO kk=1,nz-1
  SAIn = SAIn + (zz(kk+1)-zz(kk))*(yy(ny)-yy(1))
END DO

SA = SAIn + SAOut ! compute total surface area

!----- account for the villi
SA = SA - numVilliActual*(PI*Rv*Rv)						! subtract the cross sectional area of the villous bases from the total outer surface area
SA = SA + numVilliActual*(2.0_dbl*PI*Rv*(Lv-Rv))				! add the surface area from the villous cylinders
SA = SA + numVilliActual*(2.0_dbl*PI*Rv*Rv)					! add the surface area from the villous tips

!----- open and write to a file
IF(iter .GT. 0) THEN
  WRITE(2474,'(2E25.15)') REAL(iter/(nt/nPers)), SA,SAIn,SAOut			! write surface area to file
  CALL FLUSH(2474)
END IF

!------------------------------------------------
END SUBROUTINE SurfaceArea
!------------------------------------------------










































!--------------------------------------------------------------------------------------------------
SUBROUTINE VilliPosition		! Calculates the position of the villi at the current time step
!--------------------------------------------------------------------------------------------------
!!!!!! NOTE: FOR ALL THE VILLI STUFF, Output stuff WE USE rDomOut for the time being instead of rDom. Need to modify to include rDomIn if needed 
!!!!!! This is temporary to prevent the code from blowing up as we are not interested in using the villi at this time. 


IMPLICIT NONE

CALL VilliBase						! determine the x,y,z location of each villus base
CALL VilliMove						! specify the movement of the villi
CALL VilliTip						! determine the x,y,z locatin of each villus tip

!------------------------------------------------
END SUBROUTINE VilliPosition
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE VilliBase				! Calculates the x,y,z location of each villus base
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: nvz,nvt,n		! index variables
INTEGER(lng) :: ivz				! k location of each villus
REAL(dbl) :: vx,vy,vz,vt		! x,y,z, and theta locations of each villus

! initialize the villi counter
n = 0_lng

DO nvz=1,numVilliZ

  vz 	= L*(REAL(nvz - 0.5_dbl)/(REAL(numVilliZ)))							! z location of the villus (real)
  ivz	= NINT(vz/zcf)																	! k location of the villus (integer)												

  DO nvt=1,numVilliTheta 
      
    n = n + 1_lng																		! count the number of villi
 
    vt = (REAL(nvt - 0.5_dbl))*((0.5_dbl*PI)/(REAL(numVilliTheta)))	! theta location of the villus

    vx = rDomOut(ivz)*COS(vt)															! x location of the villus
    vy = rDomOut(ivz)*SIN(vt)															! y location of the villus

    ! store the real location of each villus
    villiLoc(n,1) = vx
    villiLoc(n,2) = vy
    villiLoc(n,3) = vz

  END DO
    
END DO

!------------------------------------------------
END SUBROUTINE VilliBase
!------------------------------------------------











!--------------------------------------------------------------------------------------------------
SUBROUTINE VilliMove						! specifies the movement of each villus
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: n							! index variable
INTEGER(lng) :: ivz						! k location of each villus
REAL(dbl) :: vx,vy,vz					! x,y,z, and theta locations of each villus
REAL(dbl) :: rL,rR						! radii of the nodes to the left and right of the current villus location
REAL(dbl) :: zL,zR						! axial locations of the nodes to the left and right of the current villus location
REAL(dbl) :: thetaR, thetaX			! angles of each villus with respect to the radius and x axis
REAL(dbl) :: time							! time

DO n=1,numVilli

  ! current time
  time = iter*tcf

  ! store thetaR and thetaX from the previous iteration
  villiLoc(n,9) = villiLoc(n,5)		! thetaR at previous timestep
  villiLoc(n,10) = villiLoc(n,4)		! thetaX at previous timestep

  ! x,y,z location of the villous base
  vx = villiLoc(n,1)						
  vy = villiLoc(n,2)
  vz = villiLoc(n,3)

  ! ---------------- passive movement from "riding" on the intestinal wall ------------------------
  ivz = vz/zcf									! k node location (along axial direction)
  rL = rDomOut(ivz-1)							! radius of the node to the left of the current villus (k-1)
  rR = rDomOut(ivz+1)							! radius of the node to the right of the current villus (k+1)
  zL = zz(ivz-1)								! axial distance of the node to the left of the current villus (k-1)
  zR = zz(ivz+1)								! axial distance of the node to the right of the current villus (k+1)
  thetaR = ATAN((rR - rL)/(zR - zL))	! angle of the villus with respect to the radius
  thetaX = ATAN(vy/vx)						! angle of the villus with respect to the x-axis
  ! -----------------------------------------------------------------------------------------------

  ! ---------------- active movement from specified villous motion --------------------------------
  IF(MOD(villiGroup(n),2_lng) .EQ. 0) THEN		! even groups
    IF(randORord .EQ. RANDOM) THEN
      thetaX = thetaX - (activeVflagT)*(villiAngle*SIN(2.0_dbl*PI*vFreqT*time + rnd(n)*2.0_dbl*PI))						! azimuthal direction (random)
      thetaR = thetaR - (activeVflagZ)*(villiAngle*SIN(2.0_dbl*PI*vFreqZ*time + (rnd(n+numVilli)*2.0_dbl*PI)))			! axial direction (random)
    ELSE IF(randORord .EQ. ORDERED) THEN
      thetaX = thetaX - (activeVflagT)*(villiAngle*SIN(2.0_dbl*PI*vFreqT*time))													! azimuthal direction (ordered) 
      thetaR = thetaR - (activeVflagZ)*(villiAngle*SIN(2.0_dbl*PI*vFreqZ*time))													! axial direction (ordered)
    ELSE
      OPEN(1000,FILE="error.txt")
      WRITE(1000,*) "Error in VilliMove in Geometry.f90 at line 535: randORord is not RANDOM(1) or ORDERED(2)..."
      WRITE(1000,*) "randORord=", randORord
      CLOSE(1000)
      STOP
    END IF
  ELSE
    IF(randORord .EQ. RANDOM) THEN		! odd groups
      thetaX = thetaX + (activeVflagT)*(villiAngle*SIN(2.0_dbl*PI*vFreqT*time + rnd(n)*2.0_dbl*PI))						! azimuthal direction (random)
      thetaR = thetaR + (activeVflagZ)*(villiAngle*SIN(2.0_dbl*PI*vFreqZ*time + (rnd(n+numVilli)*2.0_dbl*PI)))			! axial direction (random)
    ELSE IF(randORord .EQ. ORDERED) THEN
      thetaX = thetaX + (activeVflagT)*(villiAngle*SIN(2.0_dbl*PI*vFreqT*time))													! azimuthal direction (ordered) 
      thetaR = thetaR + (activeVflagZ)*(villiAngle*SIN(2.0_dbl*PI*vFreqZ*time))													! axial direction (ordered)
    ELSE
      OPEN(1000,FILE="error.txt")
      WRITE(1000,*) "Error in VilliMove in Geometry.f90 at line 535: randORord is not RANDOM(1) or ORDERED(2)..."
      WRITE(1000,*) "randORord=", randORord
      CLOSE(1000)
      STOP
    END IF
  END IF
  ! -----------------------------------------------------------------------------------------------

  ! store the angles for each villus
  villiLoc(n,4) = thetaX
  villiLoc(n,5) = thetaR
    
END DO

!------------------------------------------------
END SUBROUTINE VilliMove
!------------------------------------------------











!--------------------------------------------------------------------------------------------------
SUBROUTINE VilliTip						! Calculates the x,y,z location of each villus tip (minus the half hemisphere)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: n							! index variables
REAL(dbl) :: vx,vy,vz					! x,y,z, and theta locations of each villus
REAL(dbl) :: thetaR, thetaX			! angles of each villus with respect to the radius and x axis
REAL(dbl) :: vx2,vy2,vz2				! x,y,z location of each villus endpoint

DO n=1,numVilli

  ! x,y,z location and angles of the villus
  vx = villiLoc(n,1)						
  vy = villiLoc(n,2)
  vz = villiLoc(n,3)
  thetaX = villiLoc(n,4) 
  thetaR = villiLoc(n,5)

  ! calculate the end point of the villus (minus the hemisphere)
  vx2 = vx - (Lv-Rv)*COS(thetaR)*COS(thetaX)
  vy2 = vy - (Lv-Rv)*COS(thetaR)*SIN(thetaX)
  vz2 = vz + (Lv-Rv)*SIN(thetaR)

  ! store the location of the endpoints for each villus
  villiLoc(n,6) = vx2
  villiLoc(n,7) = vy2
  villiLoc(n,8) = vz2
    
END DO

!------------------------------------------------
END SUBROUTINE VilliTip
!------------------------------------------------











!--------------------------------------------------------------------------------------------------
SUBROUTINE SetNodesVilli					! defines the geometry of the villi via the "node" array of flags
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

INTEGER(lng)	:: n,i,j,k					! index variables
INTEGER(lng)	:: ii,jj,kk					! index variables
INTEGER(lng)   :: nvz,nvt					! index variables
INTEGER(lng)	:: iminV,imaxV				! i-indices of the current block to be checked
INTEGER(lng)	:: jminV,jmaxV				! j-indices of the current block to be checked
INTEGER(lng)	:: kminV,kmaxV				! k-indices of the current block to be checked
REAL(dbl)		:: villiVec(3)				! vector from each villus base to the top of the villus cylinder
REAL(dbl)		:: pointVec(3)				! vector from each villus base to the current point
REAL(dbl)		:: dotProd, crossProd	! dot and cross products between the current point and the current villus
REAL(dbl)		:: sinTheta, cosTheta	! sine and cosine of the angles between the two vectors
REAL(dbl)		:: magVilli, magPoint	! magnitudes of the two vectors
REAL(dbl)		:: dist						! distance from the current point to the villus
REAL(dbl)		:: Cx,Cy,Cz					! vector between villous base and point on the villus closest to the current point
REAL(dbl)		:: uV,vV,wV					! velocity used at the uncovered nodes

DO nvz=1,numVilliZ

  IF((MOD(nvz,(numVilliZ/numVilliGroups)) .NE. 0) .OR. (numVilliGroups .EQ. 1)) THEN		! skip a row of villi between groups (unless only 1 group)

    DO nvt=1,numVilliTheta

      n = (nvz-1_lng)*numVilliTheta + nvt										! villus number

      iminV = MIN(villiLoc(n,1)/xcf, villiLoc(n,6)/xcf) - INT(ANINT(1.5_dbl*Rv/xcf))
      imaxV = MAX(villiLoc(n,1)/xcf, villiLoc(n,6)/xcf) + INT(ANINT(1.5_dbl*Rv/xcf))
      jminV = MIN(villiLoc(n,2)/ycf, villiLoc(n,7)/ycf) - INT(ANINT(1.5_dbl*Rv/ycf))
      jmaxV = MAX(villiLoc(n,2)/ycf, villiLoc(n,7)/ycf) + INT(ANINT(1.5_dbl*Rv/ycf))
      kminV = MIN((villiLoc(n,3)/zcf+0.5_dbl), (villiLoc(n,8)/zcf+0.5_dbl)) - INT(ANINT(1.5_dbl*Rv/zcf))
      kmaxV = MAX((villiLoc(n,3)/zcf+0.5_dbl), (villiLoc(n,8)/zcf+0.5_dbl)) + INT(ANINT(1.5_dbl*Rv/zcf))

      DO kk=kminV,kmaxV
        DO jj=jminV,jmaxV
          DO ii=iminV,imaxV

            ! check to if the point is in the subdomain
            IF(((ii .GE. iMin-1_lng) .AND. (ii .LE. iMax+1_lng)) .AND.	&
               ((jj .GE. jMin-1_lng) .AND. (jj .LE. jMax+1_lng)) .AND.	&
               ((kk .GE. kMin-1_lng) .AND. (kk .LE. kMax+1_lng))) THEN

              ! transform into local subdomain coordinates
              i = ii - (iMin - 1_lng)
              j = jj - (jMin - 1_lng)
              k = kk - (kMin - 1_lng)

              ! ignore the solid nodes
              IF(node(i,j,k) .NE. SOLID) THEN

                ! define a vector between the villus base and the current point
                pointVec(1) = (xx(ii)-villiLoc(n,1))					! x-coordinate
                pointVec(2) = (yy(jj)-villiLoc(n,2))					! y-coordinate
                pointVec(3) = (zz(kk)-villiLoc(n,3))					! z-coordinate

                ! define a vector between the villus base and the top of the villus cylinder
                villiVec(1) = (villiLoc(n,6)-villiLoc(n,1))			! x-coordinate
                villiVec(2) = (villiLoc(n,7)-villiLoc(n,2))			! y-coordinate
                villiVec(3) = (villiLoc(n,8)-villiLoc(n,3))			! z-coordinate

                ! compute the dot product of villiVec and pointVec
                dotProd = villiVec(1)*pointVec(1) + villiVec(2)*pointVec(2) + villiVec(3)*pointVec(3)
          
                ! calculate the magnitudes of villiVec and pointVec
                magVilli = SQRT(villiVec(1)*villiVec(1) + villiVec(2)*villiVec(2) + villiVec(3)*villiVec(3))
                magPoint = SQRT(pointVec(1)*pointVec(1) + pointVec(2)*pointVec(2) + pointVec(3)*pointVec(3))

                ! get the cosine of the angle between the two vectors
                cosTheta = dotProd/(magVilli*magPoint)

                ! check to see if the point is above or below the top of the villus cylinder and calculate the proper distance between the point and the villus
                IF(magPoint*cosTheta .LE. magVilli) THEN				! below
                  ! calcualte the shortest distance between the point and the centerline of the villus cylinder
                  sinTheta = SQRT(1.0_dbl-cosTheta*cosTheta)		! sine of the angle between the two vectors
                  dist = magPoint*sinTheta								! distance between the point and the CL
                ELSE																! above
                  ! calculate the distance between the current point and the top of the villus cylinder (base of hemisphere)
                  dist = SQRT((xx(ii)-villiLoc(n,6))**2 + (yy(jj)-villiLoc(n,7))**2 + (zz(kk)-villiLoc(n,8))**2)
                END IF

                ! check to see if the villus is covering the node
                IF(dist .LE. Rv) THEN				! covers
                  node(i,j,k) = -n					! flag the node as being covered by the nth villus (-n)
                ELSE
                  IF(node(i,j,k) .EQ. -n) THEN
                    ! find the influence of villous velocity on the current point
                    CALL CalcC(i,j,k,n,Cx,Cy,Cz)
                    CALL VilliVelocity(n,Cx,Cy,Cz,uV,vV,wV)
!                    CALL NeighborVelocity(i,j,k,uV,vV,wV)
                    CALL SetProperties(i,j,k,uV,vV,wV)
                    node(i,j,k) = FLUID				! fluid node that was covered last time step
                  END IF
                END IF

              END IF
 
            END IF

          END DO
        END DO
      END DO

    END DO 

  END IF

END DO

!------------------------------------------------
END SUBROUTINE SetNodesVilli
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE CalcC(i,j,k,vNum,Cx,Cy,Cz)				! calculates the point, C, at which the line of influence intersects the boundary (using "ray tracing" - see wikipedia article)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k					! current node
INTEGER(lng), INTENT(IN) :: vNum						! number of the current villus
REAL(dbl), INTENT(OUT) :: Cx,Cy,Cz					! point of intersection
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Vx,Vy,Vz									! villous vector (base to tip)
REAL(dbl) :: Ix,Iy,Iz									! vector between base of the villus to the intersection point beween the vector from A to the centerline hits the centerline
REAL(dbl) :: AHx,AHy,AHz								! vector between current node and center of the villous hemisphere
REAL(dbl) :: AIx,AIy,AIz								! vector between current node and centerline of the villous
REAL(dbl) :: AHMagn,AIMagn								! magnitudes of AH and AI vectors
REAL(dbl) :: dxV,dyV,dzV								! unit vector pointing in the direction of the villous vector
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: AP											! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: dotAV,AMag,cosThetaAV					! dot product of A and V, magnitudes of A and V, cosine of the angle between A and V
REAL(dbl) :: VMag											! magnitudes of V

! find the vector between the villous base and the current node (shift coordinate system)
Ax = (x(i) - villiLoc(vNum,1))						! x-coordinate
Ay = (y(j) - villiLoc(vNum,2))						! y-coordinate
Az = (z(k) - villiLoc(vNum,3))						! z-coordinate

! find the vector between the villous base and the villous tip
Vx = (villiLoc(vNum,6)-villiLoc(vNum,1))			! x-coordinate
Vy = (villiLoc(vNum,7)-villiLoc(vNum,2))			! y-coordinate
Vz = (villiLoc(vNum,8)-villiLoc(vNum,3))			! z-coordinate

! determine if the current point is above or below the top of the villous cylinder (in the villous frame)
!		and calculate the distance from the current node to the wall accordingly
dotAV	= Ax*Vx + Ay*Vy + Az*Vz							! dot product of A and V
AMag	= SQRT(Ax*Ax + Ay*Ay + Az*Az)					! magnitude of A
VMag	= SQRT(Vx*Vx + Vy*Vy + Vz*Vz)					! magnitude of V
cosThetaAV = dotAV/(AMag*VMag)						! cosine of angle between A and V
IF(AMag*cosThetaAV .GE. VMag) THEN					! the node is above the cylindrical part of the villus (in the villous frame)
  ! find the vector and unit vector from point A to center of hemisphere
  AHx 	= Vx - Ax										! x-coordinate
  AHy 	= Vy - Ay										! y-coordinate
  AHz 	= Vz - Az										! z-coordinate
  AHMagn	= SQRT(AHx*AHx + AHy*AHy + AHz*AHz)		! magnitude of AH
  dx 		= AHx/AHMagn									! x-coordinate
  dy		= AHy/AHMagn									! y-coordinate
  dz		= AHz/AHMagn									! z-coordinate
  ! find the distance between the the current node and the villus  
  AP = AHMagn - Rv										! distance from the current node (point A) to the villus (point P)		
ELSE															! solid node is in the cylinder
  ! find the point of intersection between the vector pointing from the current node (A) to the villous centerline, and the vector pointing from the villous base to villous tip
  ! unit vector in the direction from the villous base to villous tip
  dxV = Vx/VMag	
  dyV = Vy/VMag
  dzV = Vz/VMag
  ! intersection point
  Ix = villiLoc(vNum,1) + dxV*(AMag*cosThetaAV)
  Iy = villiLoc(vNum,2) + dyV*(AMag*cosThetaAV)
  Iz = villiLoc(vNum,3) + dzV*(AMag*cosThetaAV)
  ! find the vector between A and I and its magnitude
  AIx = Ix - Ax
  AIy = Iy - Ay
  AIz = Iz - Az
  AIMagn = SQRT(AIx*AIx + AIy*AIy + AIz*AIz)
  ! find the unit vector between A and I
  dx = AIx/AIMagn
  dy = AIy/AIMagn
  dz = AIz/AIMagn			
  ! find the distance between the the current node and the villus  
  AP = AIMagn - Rv		
END IF

! find the vector from the base of the villus to the point P		
Cx = Ax + AP*dx											! x-coordinate
Cy = Ay + AP*dy											! y-coordinate
Cz = Az + AP*dz											! z-coordinate	
	
!------------------------------------------------
END SUBROUTINE CalcC
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE NeighborVelocity(i,j,k,ubx,uby,ubz)	! calculate the average neighboring node velocity
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

INTEGER(lng), INTENT(IN) :: i,j,k				! current node location
REAL(dbl), INTENT(OUT) :: ubx,uby,ubz			! average velocity of the neighboring nodes
INTEGER(lng)	:: m,ii,jj,kk						! index variables
INTEGER(lng)	:: numFLUIDs						! number of fluid nodes
REAL(dbl)		:: uSum,vSum,wSum					! sum of the velocities of the neighboring fluid nodes
CHARACTER(7)	:: iter_char						! iteration stored as a character

! initialize the quantities
uSum = 0.0_dbl
vSum = 0.0_dbl
wSum = 0.0_dbl
numFLUIDs = 0_lng

! calculate the average velocity of the current node's neighbors
DO m=1,NumDistDirs

  ii = i + ex(m)
  jj = j + ey(m)
  kk = k + ez(m)

  IF(((ii .GE. 0) .AND. (ii .LE. nxSub+1_lng)) .AND.	&
     ((jj .GE. 0) .AND. (jj .LE. nySub+1_lng)) .AND.	&
     ((kk .GE. 0) .AND. (kk .LE. nzSub+1_lng))) THEN

    IF(node(ii,jj,kk) .EQ. FLUID) THEN
      uSum = uSum + u(ii,jj,kk)
      vSum = vSum + v(ii,jj,kk)
      wSum = wSum + w(ii,jj,kk)    
      numFLUIDs = numFLUIDs + 1_lng     
    END IF       

  END IF

END DO

IF(numFLUIDs .NE. 0_lng) THEN

  ubx = uSum/numFluids								
  uby = vSum/numFluids														
  ubz = wSum/numFluids	

ELSE

  WRITE(iter_char(1:7),'(I7.7)') iter

  OPEN(6679,FILE='errorG2-'//iter_char//'-'//sub//'.txt')

  WRITE(6679,*) 'iter', iter
  WRITE(6679,*) 'i,j,k:', i,j,k
  WRITE(6679,*) 'node(i,j,k)', node(i,j,k)
  WRITE(6679,*) 'u(i,j,k)', u(i,j,k)  
  WRITE(6679,*) 'v(i,j,k)', v(i,j,k)  
  WRITE(6679,*) 'w(i,j,k)', w(i,j,k)  
  WRITE(6679,*) 'numFLUIDs', numFLUIDs
  WRITE(6679,*)
  WRITE(6679,*)

  DO m=1,NumDistDirs

    ii = i + ex(m)
    jj = j + ey(m)
    kk = k + ez(m)  

    IF(((ii .GE. 0) .AND. (ii .LE. nxSub+1_lng)) .AND.	&
       ((jj .GE. 0) .AND. (jj .LE. nySub+1_lng)) .AND.	&
       ((kk .GE. 0) .AND. (kk .LE. nzSub+1_lng))) THEN
   
      WRITE(6679,*) 'ii,jj,kk:', ii,jj,kk
      WRITE(6679,*) 'node(ii,jj,kk)', node(ii,jj,kk)
      WRITE(6679,*) 'rho(ii,jj,kk)', rho(ii,jj,kk)
      WRITE(6679,*)

    ELSE

      WRITE(6679,*) '(ii,jj,kk) is out of bounds'
      WRITE(6679,*) 'ii,jj,kk:', ii,jj,kk
      WRITE(6679,*) 'imin',imin
      WRITE(6679,*) 'imax',imax
      WRITE(6679,*) 'jmin',jmin
      WRITE(6679,*) 'jmax',jmax
      WRITE(6679,*) 'kmin',kmin
      WRITE(6679,*) 'kmax',kmax
      WRITE(6679,*) 'node(ii,jj,kk)', node(ii,jj,kk)
      WRITE(6679,*) 'rho(ii,jj,kk)', rho(ii,jj,kk)
      WRITE(6679,*)

    END IF
  
  END DO

  ubx = 0.0_dbl								
  uby = 0.0_dbl														
  ubz = 0.0_dbl	

  CLOSE(6679)

END IF

!------------------------------------------------
END SUBROUTINE NeighborVelocity
!------------------------------------------------





!================================================
END MODULE Geometry
!================================================
