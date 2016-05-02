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

! Define the lattice <=> physical conversion factors
IF(domaintype .EQ. 0) THEN
        xcf 		= (0.5_lng*D_x)/(nx-1_lng)	! length conversion factor: x-direction
        ycf 		= (0.5_lng*D_y)/ny		! length conversion factor: y-direction
ELSE
        ! begin Balaji added
        xcf 		= (1.0_lng*D_x)/(nx-1_lng)	! length conversion factor: x-direction
        ycf 		= (1.0_lng*D_y)/ny		! length conversion factor: y-direction
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

!----- Initialize Variables
time 	= 0.0_dbl				! time					
h1In 	= 0.0_dbl				! mode 1 height
h1Out 	= 0.0_dbl				! mode 2 height
rDomIn	= 0.0_dbl				! summed height
rDomOut	= 0.0_dbl				! summed height

time	= iter*tcf
DO i=1,nz
   h1Out(i) = -0.38_dbl*D_x + L*0.125*(1.0_dbl + cos(2.0_dbl*z(i)*pi/L) ) + 5.0000e-5 + s1*time 	   
   h1In(i)  = -0.48_dbl*D_x + L*0.125*(1.0_dbl + cos(2.0_dbl*z(i)*pi/L) ) + 5.0000e-5 + s1*time
END DO

!----- since PI cannot be stored exactly, the wavelength(s) does/do not EXACTLY span the domain...
!----- set h1(nz) to h1(0) and h1(nz+1) to h(1) to ensure periodicity
h1In(nz+1)   = h1In(1)
h1In(0) = h1In(nz)
h1Out(nz+1)   = h1Out(1)
h1Out(0) = h1Out(nz)

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

  INTEGER(lng),INTENT(IN) :: i,j,k       ! current node location
  REAL(dbl)   ,INTENT(IN) :: ubx,uby,ubz ! velocity of the boundary
  INTEGER(lng)    :: m           ! Direction index
  INTEGER(lng)    :: im1, jm1, km1         ! Node inside the wall
  INTEGER(lng)    :: ip1,jp1,kp1,iB,jB,kB  ! First neighboring node location
  INTEGER(lng)    :: ip2,jp2,kp2,iC,jC,kC      ! Second neighboring node location
  CHARACTER(7)    :: iter_char            ! iteration stored as a character
  REAL(dbl)       :: feq                  ! equilibrium distribution function
  REAL(dbl):: Geom_norm_x,Geom_norm_y,Geom_norm_z! geometry normal vector
  REAL(dbl):: q,n_prod,n_prod_max

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
        ip1= i + ex(m)
        jp1= j + ey(m)
        kp1= k + ez(m)
        im1= i - ex(m)
        jm1= j - ey(m)
        km1= k - ez(m)
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
     
     CALL qCalc(m,i,j,k,im1,jm1,km1,q)
     
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

     q= 1.0_dbl! approximating that uncovered node is at the wall (it is close enough ...)
     phi(i,j,k)= ( (phi(iB,jB,kB)*(1.0+q)*(1.0+q)/(1.0+2.0*q))&
          - (phi(iC,jC,kC)*q*q/(1.0+2.0*q)) &
          - (q*(1+q)/(1+2.0*q))*(coeffConst/coeffGrad) ) &
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
