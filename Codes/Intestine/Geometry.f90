!==================================================================================================
MODULE Geometry	! Defines the geometry for the simulation
						! Subroutines (NodeFlags, BoundaryVelocity)
!==================================================================================================
USE SetPrecision      
USE Setup
USE LBM
USE ICBC
USE MPI

IMPLICIT NONE 

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE Geometry_Setup					! sets up the geometry
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER :: isize,idate(8)									! size of seed array for the random number genreator, array for output of DATE_AND_TIME
INTEGER,ALLOCATABLE  :: iseed(:)							! seeds for random number generator
INTEGER(lng) :: i,j,k,kk,iCon,it,iPer,nPers_INT		! index variables
INTEGER(lng) :: nvz,nvt,n,g								! index variables
INTEGER(lng) :: mpierr										! MPI standard error variable 
REAL(dbl) :: macroFreq										! macroscopic contraction frequency
INTEGER(lng) :: xaxis,yaxis								! axes index variables

! Define the lattice <=> physical conversion factors
IF(domaintype .EQ. 0) THEN
        xcf 		= (0.5_lng*D)/(nx-1_lng)		! length conversion factor: x-direction
        ycf 		= (0.5_lng*D)/(ny-1_lng)		! length conversion factor: y-direction
ELSE
        ! begin Balaji added
        xcf 		= (1.0_lng*D)/(nx-1_lng)		! length conversion factor: x-direction
        ycf 		= (1.0_lng*D)/(ny-1_lng)		! length conversion factor: y-direction
        ! end Balaji added
ENDIF

zcf 		= L/nz								! length conversion factor: z-direction
tcf 		= nuL*((xcf*xcf)/nu)				! time conversion factor
dcf 		= den/denL							! density conversion factor
vcf 		= xcf/tcf							! velocity conversion factor
pcf 		= cs*cs*vcf*vcf					! pressure conversion factor

! Determine the number of time steps to run
nt = ANINT((nPers*Tmix)/tcf)

! Initialize arrays
node		= -99_lng							! node flag array
rDom		= 0.0_dbl							! radius at each z-location
r			= 0.0_dbl							! temporary radius array for entire computational domain
velDom	= 0.0_dbl							! wall velocity at each z-location (global)
vel		= 0.0_dbl							! wall velocity at each z-location (local)

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
      ! begin Balaji added 
      !INTEGER(lng) :: xaxis,yaxis								! axes index variables
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
      ! end Balaji added 
ENDIF
! Mode 1 - Peristalsis
a1				= (0.5_dbl*D)/(2.0_dbl - epsOVERa1)					! mean half-width of wave1
eps1 			= epsOVERa1*a1												! occlusional distance
lambda1		= L/numw1													! wavelength
aOVERlam1	= a1/lambda1												! ratio of mean half-width to wavelength 
kw1			= (2.0_dbl*PI)/lambda1									! wave number
amp1			= 0.5_dbl*((0.5_dbl*D)-eps1)							! amplitude of the wave
Tp				= lambda1/s1												! peristaltic period
Re1			= ((s1*(0.5_dbl*D))/nu)*((0.5_dbl*D)/lambda1)	! Reynolds number based on mode 1

! Mode 2 - Segmental Contractions
a2				= (0.5_dbl*D)/(2.0_dbl - epsOVERa2)					! mean half-width of wave1 (based on peristalsis definition)
eps2 			= epsOVERa2*a2												! occlusional distance
lambda2		= L/numw2													! wavelength (physical units)
nlambda2		= nz/numw2													! wavelength (nodes)
aOVERlam2	= a2/lambda2												! ratio of mean half-width to wavelength 
amp2			= 0.5_dbl*((0.5_dbl*D)-eps2)							! amplitude of the wave
shift2		= 0.5_dbl*((0.5_dbl*D)+eps2)							! amplitude of the wave
segment		= nlambda2/6_lng											! length of each segment of the segmental wave   !!!!! CAREFUL HERE WITH SYMMETRY!
seg1L			= 1_lng + segment											! left point of sloped segement 1
seg1R			= 1_lng + 2_lng*segment									! right point of sloped segement 1
seg2R			= nlambda2 - segment										! right point of sloped segement 2
seg2L			= nlambda2 - (2_lng*segment)							! left point of sloped segement 2
s2				= (0.5_dbl*D)/Ts											! speed of collapse fo segmental contraction
Re2			= (s2*(0.5_dbl*D))/nu									! Reynolds number based on mode 2

! Allocate and initialize the villi arrays
numVilli				= numVilliZ*numVilliTheta						! determine the total number of villi
IF(numVilliGroups .GT. 1_lng) THEN
  numVilliActual	= numVilli - numVilliGroups*numVilliTheta	! determine the actual number of villi (subracting those skipping in grouping)
ELSE
  numVilliActual = numVilli											! if there is only 1 group, numVilli is numVilliActual
END IF
ALLOCATE(villiLoc(numVilli,10))										! location and other information of the villi
villiLoc = 0.0_dbl

! Set up the villi groups
ALLOCATE(villiGroup(numVilli))										! array of which group the villi are in
villiGroup = 0.0_dbl
n = 0_lng
DO nvz=1,numVilliZ
  DO nvt=1,numVilliTheta

    n = n + 1_lng															! villus number

    DO g=1,numVilliGroups

      IF((n .GT. ((g-1_lng)*numVilliTheta*(numVilliZ/numVilliGroups))) .AND.		&
         (n .LE. (g*numVilliTheta*(numVilliZ/numVilliGroups)))) THEN
        villiGroup(n) = g
      END IF

    END DO

  END DO
END DO

! Write villiGroup to a test file
!OPEN(173,FILE='villiGroup-'//sub//'.dat')
!  DO n=1,numVilli
!    WRITE(173,*) 'n =', n, 'villiGroup(n)=', villiGroup(n)
!  END DO
!CLOSE(173)
!STOP

! Convert villous length and radius of the villi to meters
Lv = Lv*(0.000001_dbl)
Rv = Rv*(0.000001_dbl)	

! Determine villous frequency
macroFreq = 1.0_dbl/Tmix
vFreqT = freqRatioT*macroFreq
vFreqZ = freqRatioZ*macroFreq

IF(freqRatioT .LT. 0.00000001_dbl) THEN		! if the frequency ratio = 0 (no active villous motion) then zero out the contribution
  activeVflagT = 0.0_dbl
ELSE
  activeVflagT = 1.0_dbl
END IF

IF(freqRatioZ .LT. 0.00000001_dbl) THEN		! if the frequency ratio = 0 (no active villous motion) then zero out the contribution
  activeVflagZ = 0.0_dbl
ELSE
  activeVflagZ = 1.0_dbl
END IF

! Convert villiAngle from degrees to radians
villiAngle = (villiAngle/180.0_dbl)*PI

!IF(restart .EQ. .FALSE.) THEN
IF(restart .EQV. .FALSE.) THEN
!IF(restart .eq. FALSE) THEN

  IF(randORord .EQ. RANDOM) THEN
    ! Fill out the random array for the villous oscilliatory phases
    ALLOCATE(rnd(2_lng*numVilli))																		! allocate the array of random numbers for random villi phase angles
    IF(myid .EQ. master) THEN

      CALL DATE_AND_TIME(VALUES=idate)																	! get the date and time (for more different seeds each time)
      CALL RANDOM_SEED(SIZE=isize)																		! get the size of the seed array
      ALLOCATE(iseed(isize))																				! allocate the seed array
      CALL RANDOM_SEED(GET=iseed)																		! get the seed array
      iseed = iseed*(idate(8)-500_lng)    															! idate(8) contains millisecond
      CALL RANDOM_SEED(PUT=iseed)																		! use the seed array
      CALL RANDOM_NUMBER(rnd)																				! get the actual random numbers
      DEALLOCATE(iseed)

      ! print the rnd array for restarting
      OPEN(1777,FILE='rnd.dat')
      DO i=1,2_lng*numVilli
        WRITE(1777,*) rnd(i)
      END DO
      CLOSE(1777)

    END IF

    CALL MPI_BCAST(rnd,2_lng*numVilli,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpierr)		! send/recv rnd on all processing units

  END IF

  ! Initialize the Geometry
  CALL AdvanceGeometry

END IF

!------------------------------------------------
END SUBROUTINE Geometry_Setup
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE AdvanceGeometry												! advances the geometry in time
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

! Calculate the radius at the current time step
CALL BoundaryPosition
CALL VilliPosition

! Calculate the velocity at boundary point
CALL BoundaryVelocity

! Flag the fluid/solid nodes based on the new geometry
CALL SetNodes

!------------------------------------------------
END SUBROUTINE AdvanceGeometry
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE BoundaryPosition		! Calculates the position of the wall at the current time step
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

REAL(dbl) :: h1(0:nz+1)				! Mode 1 (peristalsis)
REAL(dbl) :: h2(0:nz+1)				! Mode 2	(segmental)
REAL(dbl) :: Ac, lambdaC, shiftC	! temporary variables for the cos slopes
REAL(dbl) :: time						! time
INTEGER(lng) :: i,j,ii,k			! indices

! Initialize Variables
time 	= 0.0_dbl						! time					
!rDom	= 0.5_dbl*D						! summed height
!h1 	= 0.5_dbl*D						! mode 1 height
!h2 	= 0.5_dbl*D						! mode 2 height
h1 	= 0.0_dbl						! mode 1 height
h2 	= 0.0_dbl						! mode 2 height
rDom	= 0.0_dbl						! summed height

! Current Physical Time
time	= iter*tcf

!------------------------- Mode 1 - peristalsis -----------------------------
DO i=0,nz-1

  h1(i) 	= amp1*(COS(kw1*(zz(i) - (s1*time)))) + (0.5_dbl*D - amp1)
!! Yanxing's expression
!  h1(i)         = amp1*sin(2.0_dbl*PI*((real(i,dbl)-0.5_dbl)/real(nz,dbl)-0.1_dbl*iter/real(nz,dbl))+pi/2.0_dbl)+ (0.5_dbl*D - amp1)

END DO

! since PI cannot be stored exactly, the wavelength(s) does/do not EXACTLY span the domain...
! set h1(nz) to h1(0) and h1(nz+1) to h(1) to ensure periodicity
h1(nz) 	= h1(0)
h1(nz+1)= h1(1)

!IF((h1(0) .NE. h1(nz)) .OR. (h1(nz+1) .NE. h1(1))) THEN
!  WRITE(6678,*) 'h1(0)   ', h1(0)
!  WRITE(6678,*) 'h1(nz)  ', h1(nz)
!  WRITE(6678,*) 'h1(1)   ', h1(1)
!  WRITE(6678,*) 'h1(nz+1)', h1(nz+1)
!  WRITE(6678,*) 'start', (zz(0) + zz(1))/2.0_dbl
!  WRITE(6678,*) 'end', (zz(nz) + zz(nz+1))/2.0_dbl
!  WRITE(6678,*) 'zz(0)   ', zz(0)
!  WRITE(6678,*) 'zz(1)   ', zz(1)
!  WRITE(6678,*) 'zz(nz)  ', zz(nz)
!  WRITE(6678,*) 'zz(nz+1)', zz(nz+1)
!  STOP
!END IF
!----------------------------------------------------------------------------

!------------------- Mode 2 - segmental contractions ------------------------

! Calculate the geometry for the first wave
! First Straight Piece
DO i=0,seg1L

  h2(i) = amp2*(COS(((2.0_dbl*PI)/Ts)*time)) + shift2
  
END DO

! Second Straight Piece
DO i=seg1R,seg2L

  h2(i) = amp2*(COS(((2.0_dbl*PI)/Ts)*(time-(Ts/2.0_dbl)))) + shift2
  
END DO

! Third Straight Piece
DO i=seg2R,nlambda2+1

  h2(i) = amp2*(COS(((2.0_dbl*PI)/Ts)*time)) + shift2
  
END DO

! First Cos Piece
Ac	= 0.5_dbl*(h2(seg1L)-h2(seg1R))
lambdaC	= 2.0_dbl*(zz(seg1L)-zz(seg1R))
shiftC	= 0.5_dbl*(h2(seg1L)+h2(seg1R))
DO i=seg1L+1,seg1R-1

  h2(i) = Ac*COS((2.0_dbl*PI/lambdaC)*(zz(i)-zz(seg1L))) + shiftC
  
END DO

! Second Cos Piece
Ac			= 0.5_dbl*(h2(seg2L)-h2(seg2R))
lambdaC	= 2.0_dbl*(zz(seg2L)-zz(seg2R))
shiftC	= 0.5_dbl*(h2(seg2L)+h2(seg2R))
DO i=seg2L+1,seg2R-1

  h2(i) = Ac*COS((2.0_dbl*PI/lambdaC)*(zz(i)-zz(seg2L))) + shiftC
  
END DO

!IF((h2(0) .NE. h2(nz)) .OR. (h2(nz+1) .NE. h2(1))) THEN
!  WRITE(6678,*) 'line 214'
!  WRITE(6678,*) 'h2(0)', h2(0)
!  WRITE(6678,*) 'h2(nz)', h2(nz)
!  WRITE(6678,*) 'h2(1)', h2(1)
!  WRITE(6678,*) 'h2(nz+1)', h2(nz+1)
!  STOP
!END IF

! Repeat for the rest of the waves
DO j=1,(numw2-1)
  DO i=0,nlambda2+1

    ii = i + j*nlambda2
    h2(ii) = h2(i)

  END DO
END DO

!IF((h2(0) .NE. h2(nz)) .OR. (h2(nz+1) .NE. h2(1))) THEN
!  WRITE(6678,*) 'line 233'
!  WRITE(6678,*) 'h2(0)', h2(0)
!  WRITE(6678,*) 'h2(nz)', h2(nz)
!  WRITE(6678,*) 'h2(1)', h2(1)
!  WRITE(6678,*) 'h2(nz+1)', h2(nz+1)
!  STOP
!END IF

! "fudging" to make sure that the whole domain is filled (and periodic) - more logic (and computational expense would be
! necessary to do this correctly: ideally, one would determine if an even or odd number of waves was specified
! and then work from either end, and meet in the middle to ensure a symetric domain...
!h2(nz-1:nz+1) = h2(1)

!----------------------------------------------------------------------------

!-------------------------------- Mode Sum  ---------------------------------

! Sum the modes in a weighted linear combination
DO i=0,nz+1
  rDom(i) = wc1*h1(i) + wc2*h2(i)
END DO

!IF((rDom(0) .NE. rDom(nz)) .OR. (rDom(nz+1) .NE. rDom(1))) THEN
!  WRITE(6678,*) 'rDom(0)', rDom(0)
!  WRITE(6678,*) 'rDom(nz)', rDom(nz)
!  WRITE(6678,*) 'rDom(1)', rDom(1)
!  WRITE(6678,*) 'rDom(nz+1)', rDom(nz+1)
!  STOP
!END IF

!----------------------------------------------------------------------------

! Fill out the local radius array
r(0:nzSub+1) = rDom(kMin-1:kMax+1)

!IF(iter .EQ. 1) THEN
!  OPEN(697,FILE='r-'//sub//'.dat')
!  WRITE(697,*) 'VARIABLES = z, "r(z)"'
!  CLOSE(697)
!END IF
!
!OPEN(697,FILE='r-'//sub//'.dat',POSITION='APPEND')
!WRITE(697,*) 'ZONE T="', (iter*tcf)/Tmix, '" I=', nzSub+2,' F=POINT'
!
!DO k=0,nzSub+1
!  WRITE(697,*) z(k), r(k)
!END DO
!CLOSE(697)

IF(myid .EQ. master) THEN

!  ! print the radius as a function of time in 3 locations
!  IF(iter .EQ. iter0) THEN
!    OPEN(648,FILE='rLocs.dat')
!    WRITE(648,'(A50)') 'VARIABLES = "period", "rL", "rC", "rR"'
!    WRITE(648,'(A20)') 'ZONE F=POINT'
!    CLOSE(648)
!  END IF
!
!  OPEN(648,FILE='rLocs.dat',POSITION='APPEND')
!  WRITE(648,'(4E25.15)') REAL(iter/(nt/nPers)), rDom(1), rDom(nz/2), rDom(nz)
!  CLOSE(648)

!  ! print the radius along the z-axis periodically in time.
!  IF((MOD(iter,(nt/numOuts)) .EQ. 0) .OR. (iter .EQ. phiStart) .OR. (iter .EQ. nt)) THEN
!
!    IF(iter .EQ. iter0) THEN
!      OPEN(697,FILE='rZones.dat')
!      WRITE(697,*) 'VARIABLES = z, r'
!      CLOSE(697)
!    END IF
!
!    OPEN(697,FILE='rZones.dat',POSITION='APPEND')
!    WRITE(697,*) 'ZONE T="', (iter*tcf)/Tmix, '" I=', nzSub+2,' F=POINT'
!
!    DO kk=0,nzSub+1
!      WRITE(697,*) zz(kk), rDom(kk)
!    END DO
!    CLOSE(697)
!
!  END IF

  ! calculate the surface area
  CALL SurfaceArea

END IF

!------------------------------------------------
END SUBROUTINE BoundaryPosition
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SurfaceArea			! calculate the surface area at the current time and write it to a file
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

REAL(dbl) :: SA					! surface area
REAL(dbl) :: r2,r1,z2,z1		! radius and z-location for each set of consecutive points
INTEGER(lng) :: kk				! index variable

! initialize the surface area
SA = 0.0_dbl

! approximate the surface area as if the nodes were connected linearly with the neighboring nodes
! surface area between left phantom node and 1st domain node
r1 = 0.5_dbl*(rDom(0) + rDom(1))
z1 = 0.5_dbl*(zz(0) + zz(1))
r2 = rDom(1)
z2 = zz(1)
SA = -PI*(r1 + r2)*SQRT(1.0_dbl + ((r1-r2)/(z1-z2))**2)*(z1 - z2)

! surface area between right phantom node and last domain node
r1 = rDom(nz)
z1 = zz(nz)
r2 = 0.5_dbl*(rDom(nz) + rDom(nz+1))
z2 = 0.5_dbl*(zz(nz) + zz(nz+1))
SA = SA - PI*(r1 + r2)*SQRT(1.0_dbl + ((r1-r2)/(z1-z2))**2)*(z1 - z2)

! interior domain nodes
DO kk=1,nz-1
  r1 = rDom(kk)
  r2 = rDom(kk+1)
  z1 = zz(kk)
  z2 = zz(kk+1)
  SA = SA - PI*(r1 + r2)*SQRT(1.0_dbl + ((r1-r2)/(z1-z2))**2)*(z1 - z2)
END DO

! account for the villi
SA = SA - numVilliActual*(PI*Rv*Rv)						! subtract the cross sectional area of the villous bases from the total outer surface area
SA = SA + numVilliActual*(2.0_dbl*PI*Rv*(Lv-Rv))				! add the surface area from the villous cylinders
SA = SA + numVilliActual*(2.0_dbl*PI*Rv*Rv)					! add the surface area from the villous tips

! open and write to a file
IF(iter .GT. 0) THEN
  WRITE(2474,'(2E25.15)') REAL(iter/(nt/nPers)), SA			! write surface area to file
  CALL FLUSH(2474)
END IF

!------------------------------------------------
END SUBROUTINE SurfaceArea
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE VilliPosition		! Calculates the position of the villi at the current time step
!--------------------------------------------------------------------------------------------------
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

    vx = rDom(ivz)*COS(vt)															! x location of the villus
    vy = rDom(ivz)*SIN(vt)															! y location of the villus

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
  rL = rDom(ivz-1)							! radius of the node to the left of the current villus (k-1)
  rR = rDom(ivz+1)							! radius of the node to the right of the current villus (k+1)
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
SUBROUTINE BoundaryVelocity			! defines the velocity of the solid boundaries (fills "ub", "vb", and "wb" arrays)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

REAL(dbl) :: v1(0:nz+1), v2(0:nz+1)	! velocity arrays for each mode
REAL(dbl) :: lambdaC						! wavelength of the cos segments (mode 2)
REAL(dbl) :: time							! time
INTEGER(lng) :: i,j,ii					! indices

! Initialize Variables
time		= 0.0_dbl						! time
velDom	= 0.0_dbl						! summed velocity
v1			= 0.0_dbl						! mode 1 velocity
v2 		= 0.0_dbl						! mode 2 velocity				

! Current Physical Time
time = iter*tcf

!------------------------- Mode 1 - peristalsis -----------------------------
!DO i=1,nz
DO i=0,nz-1 ! Balaji added to ensure periodicity just like in h1. 

  v1(i)	= kw1*s1*amp1*(SIN(kw1*(zz(i) - (s1*time))))
  !v1(i)	= -kw1*s1*amp1*(SIN(kw1*(zz(i) - (s1*time))))
!! Yanxing's expression
!  v1(i)         = -kw1*s1*amp1*cos(2.0_dbl*PI*((real(i,dbl)-0.5_dbl)/real(nz,dbl)-0.1_dbl*iter/real(nz,dbl))+pi/2.0_dbl)

END DO

! Balaji added
!v1(0)=v1(nz)
!v1(nz+1)=v1(1)
v1(nz)=v1(0)
v1(nz+1)=v1(1)
!----------------------------------------------------------------------------

!------------------- Mode 2 - segmental contractions  -----------------------

! Calculate the wall velocity for the first wave
! First Straight Piece
DO i=0,seg1L

  v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*time))*((2.0_dbl*PI)/Ts)
  
END DO

! Second Straight Piece
DO i=seg1R,seg2L

  v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*(time-(Ts/2.0_dbl))))*((2.0_dbl*PI)/Ts)
  
END DO

! Third Straight Piece
DO i=seg2R,nlambda2

  v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*time))*((2.0_dbl*PI)/Ts)
  
END DO

! First Cos Piece
lambdaC	= 2.0_dbl*(zz(seg1L)-zz(seg1R))
DO i=seg1L+1,seg1R-1

  v2(i) = (0.5_dbl*(v2(seg1L)-v2(seg1R)))*COS((2.0_dbl*PI/lambdaC)*(zz(i)-zz(seg1L))) &
        + (0.5_dbl*(v2(seg1L)+v2(seg1R)))
    
END DO

! Second Cos Piece
lambdaC	= 2.0_dbl*(zz(seg2L)-zz(seg2R))
DO i=seg2L+1,seg2R-1

  v2(i) = (0.5_dbl*(v2(seg2L)-v2(seg2R)))*COS((2.0_dbl*PI/lambdaC)*(zz(i)-zz(seg2L))) &
        + (0.5_dbl*(v2(seg2L)+v2(seg2R)))

END DO

! Repeat for the rest of the waves
DO j=1,(numw2-1)
  DO i=1,nlambda2+1

    ii = i + j*nlambda2
    v2(ii) = v2(i)

  END DO
END DO

! "fudging" to make sure that the whole domain is filled (and periodic) - more logic (and computational expense would be
! necessary to do this correctly: ideally, one would determine if an even or odd number of waves was specified
! and then work from either end, and meet in the middle to ensure a symetric domain...
v2(nz-1:nz+1) = v2(1)

!----------------------------------------------------------------------------

!-------------------------------- Mode Sum  ---------------------------------

! Sum the modes in a weighted linear combination
DO i=0,nz+1
  velDom(i) = wc1*v1(i) + wc2*v2(i)
END DO

!----------------------------------------------------------------------------

! Fill out the local velocity array
vel(0:nzSub+1) = velDom(kMin-1:kMax+1)/vcf

!IF(iter .EQ. 1) THEN
!  OPEN(698,FILE='vel-'//sub//'.dat')
!  WRITE(698,*) 'VARIABLES = z, "vel(z)"'
!  CLOSE(698)
!END IF
!
!OPEN(698,FILE='vel-'//sub//'.dat',POSITION='APPEND')
!WRITE(698,*) 'ZONE T="', (iter*tcf)/T, '" I=', nzSub+2,' F=POINT'
!
!DO kk=0,nzSub+1
!  WRITE(698,*) z(kk), vel(kk)
!END DO
!CLOSE(698)  

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

! Flag the interior nodes and give values to nodes that just came in
DO k=1,nzSub
  DO j=1,nySub
    DO i=1,nxSub

      rijk = SQRT(x(i)*x(i) + y(j)*y(j))

      IF(rijk .LT. r(k)) THEN

        IF(node(i,j,k) .EQ. SOLID) THEN														! just came into the domain
          
          ! calculate the wall velocity (boundary)

          ubx = vel(k)*(x(i)/rijk)
          uby = vel(k)*(y(j)/rijk)
          ubz = 0.0_dbl
          
	  !! Balaji added
	  !CALL NeighborVelocity(i,j,k,ubx,uby,ubz)
	  !IF (ubx.EQ.0.0_dbl .AND. uby.EQ.0.0_dbl) THEN
          !ubx = vel(k)*(x(i)/rijk)
          !uby = vel(k)*(y(j)/rijk)
          !ubz = 0.0_dbl
	  !ENDIF

          CALL SetProperties(i,j,k,ubx,uby,ubz)

        END IF
	
        node(i,j,k)	= FLUID																		! reset the SOLID node that just came in to FLUID

      ELSE

        node(i,j,k) = SOLID																		! if rijk is GT r(k) then it's a SOLID node

	!! Balaji added
	!rho(i,j,k)=0.0_dbl
	!u(i,j,k)=0.0_dbl
	!v(i,j,k)=0.0_dbl
	!w(i,j,k)=0.0_dbl


      END IF

    END DO
  END DO
END DO

! Loop through the phantom nodes, and set the entity, but do not give values
! YZ Faces
DO iComm=1,2

  i = YZ_RecvIndex(OppCommDir(iComm))															! i index of the phantom nodes
 	
  DO j=0,nySub+1_lng

    rijk = SQRT(x(i)*x(i) + y(j)*y(j))

    DO k=0,nzSub+1_lng

      IF(rijk .LT. r(k)) THEN
        node(i,j,k) = FLUID																		! set the SOLID node that just came in to FLUID

	  !! Balaji added
	  !ubx = vel(k)*(x(i)/rijk)
          !uby = vel(k)*(y(j)/rijk)
          !ubz = 0.0_dbl
          !CALL SetProperties(i,j,k,ubx,uby,ubz)
      ELSE
        node(i,j,k) = SOLID																		! if rijk is GT r(k) then it's a SOLID node


	!! Balaji added
	!rho(i,j,k)=0.0_dbl
	!u(i,j,k)=0.0_dbl
	!v(i,j,k)=0.0_dbl
	!w(i,j,k)=0.0_dbl
      END IF
        
    END DO
  END DO

END DO

! ZX Faces
DO iComm=3,4

  j = ZX_RecvIndex(OppCommDir(iComm))															! j index of the phantom nodes

  DO i=0,nxSub+1_lng

    rijk = SQRT(x(i)*x(i) + y(j)*y(j))

    DO k=0,nzSub+1_lng

      IF(rijk .LT. r(k)) THEN
        node(i,j,k) = FLUID																		! set the SOLID node that just came in to FLUID
	  
	  !! Balaji added
	  !ubx = vel(k)*(x(i)/rijk)
          !uby = vel(k)*(y(j)/rijk)
          !ubz = 0.0_dbl
          !CALL SetProperties(i,j,k,ubx,uby,ubz)
      ELSE
        node(i,j,k) = SOLID																		! if rijk is GT r(k) then it's a SOLID node


	!! Balaji added
	!rho(i,j,k)=0.0_dbl
	!u(i,j,k)=0.0_dbl
	!v(i,j,k)=0.0_dbl
	!w(i,j,k)=0.0_dbl
      END IF
        
    END DO

  END DO

END DO

! XY Faces
DO iComm=5,6

  k = XY_RecvIndex(OppCommDir(iComm))															! k index of the phantom nodes

  DO j=0,nySub+1_lng
    DO i=0,nxSub+1_lng

      rijk = SQRT(x(i)*x(i) + y(j)*y(j))

      IF(rijk .LT. r(k)) THEN
        node(i,j,k) = FLUID																		! set the SOLID node that just came in to FLUID

	  !! Balaji added
	  !ubx = vel(k)*(x(i)/rijk)
          !uby = vel(k)*(y(j)/rijk)
          !ubz = 0.0_dbl
          !CALL SetProperties(i,j,k,ubx,uby,ubz)
      ELSE
        node(i,j,k) = SOLID																		! if rijk is GT r(k) then it's a SOLID node


	!! Balaji added
	!rho(i,j,k)=0.0_dbl
	!u(i,j,k)=0.0_dbl
	!v(i,j,k)=0.0_dbl
	!w(i,j,k)=0.0_dbl
      END IF

    END DO
  END DO

END DO

CALL SetNodesVilli																					!  set the villi node flags

!CALL SymmetryBC																						!	ensure symmetric node placement
!!CALL SymmetryBC_NODE																				!	ensure symmetric node placement

! Balaji added to make domain full 3D
IF(domaintype .EQ. 0) THEN  ! only needed when planes of symmetry exist
	CALL SymmetryBC																						!	ensure symmetric node placement
	!CALL SymmetryBC_NODE																				!	ensure symmetric node placement
ENDIF


!! Write the node flags to file
!IF(iter .EQ. 0) THEN
!  OPEN(699,FILE='flag-'//sub//'.dat')
!  WRITE(699,*) 'VARIABLES = "i", "j", "k", "node"'
!  CLOSE(699)
!END IF
!
!IF((MOD(iter,(nt/100))) .EQ. 0 .OR. (iter .EQ. 1)) THEN
!  OPEN(699,FILE='flag-'//sub//'.dat',POSITION='APPEND')
!  WRITE(699,'(A8,E15.5,3(A5,I4),A9)') 'ZONE T="', (iter*tcf)/Tmix, '" I=', nxSub+2,' J=', nySub+2,' K=', nzSub+2,' F=POINT'
!  DO k=0,nzSub+1
!    DO j=0,nySub+1
!      DO i=0,nxSub+1
!
!        WRITE(699,*) x(i), y(j), z(k), node(i,j,k)
!
!      END DO
!    END DO
!  END DO
!  CLOSE(699)  
!  CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)														! synchronize all processing units before next loop [Intrinsic]
!  STOP
!END IF

!------------------------------------------------
END SUBROUTINE SetNodes
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

!--------------------------------------------------------------------------------------------------
SUBROUTINE SetProperties(i,j,k,ubx,uby,ubz)	! give properties to nodes that just came into the fluid domain (uncovered)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE 

INTEGER(lng), INTENT(IN) :: i,j,k				! current node location
REAL(dbl), INTENT(IN) :: ubx,uby,ubz			! velocity of the boundary
INTEGER(lng)	:: m,ii,jj,kk						! index variables
INTEGER(lng)	:: numFLUIDs						! number of fluid nodes
REAL(dbl)		:: rhoSum, rhoTemp				! sum of the densities of the neighboring fluid nodes, pre-set density
REAL(dbl)		:: feq								! equilibrium distribution function
CHARACTER(7)	:: iter_char						! iteration stored as a character

! initialize the sum of surrounding densities
rhoSum = 0.0_dbl
numFLUIDs = 0_lng

! calculate the average density of the current node's neighbors
DO m=1,NumDistDirs

  ii = i + ex(m)
  jj = j + ey(m)
  kk = k + ez(m)

  IF(((ii .GE. 0) .AND. (ii .LE. nxSub+1_lng)) .AND.	&
     ((jj .GE. 0) .AND. (jj .LE. nySub+1_lng)) .AND.	&
     ((kk .GE. 0) .AND. (kk .LE. nzSub+1_lng))) THEN

    IF(node(ii,jj,kk) .EQ. FLUID) THEN
      rhoSum = rhoSum + rho(ii,jj,kk)
      numFLUIDs = numFLUIDs + 1_lng     
    END IF       

  END IF

END DO

! This should rarely happen...
IF(numFLUIDs .NE. 0_lng) THEN

  rho(i,j,k) = rhoSum/numFLUIDs

ELSE

!  WRITE(iter_char(1:7),'(I7.7)') iter
!
!  rhoTemp = rho(i,j,k)

  rho(i,j,k) = denL

!  OPEN(6679,FILE='errorG-'//iter_char//'-'//sub//'.txt')
!
!  WRITE(6679,*) 'iter', iter
!  WRITE(6679,*) 'i,j,k:', i,j,k
!  WRITE(6679,*) 'node(i,j,k)', node(i,j,k)
!  WRITE(6679,*) 'rhoTemp', rhoTemp
!  WRITE(6679,*) 'rho(i,j,k)', rho(i,j,k)  
!  WRITE(6679,*) 'numFLUIDs', numFLUIDs
!  WRITE(6679,*)
!  WRITE(6679,*)
!
!  DO m=1,NumDistDirs
!
!    ii = i + ex(m)
!    jj = j + ey(m)
!    kk = k + ez(m)  
!
!    IF(((ii .GE. 0) .AND. (ii .LE. nxSub+1_lng)) .AND.	&
!       ((jj .GE. 0) .AND. (jj .LE. nySub+1_lng)) .AND.	&
!       ((kk .GE. 0) .AND. (kk .LE. nzSub+1_lng))) THEN
!   
!      WRITE(6679,*) 'ii,jj,kk:', ii,jj,kk
!      WRITE(6679,*) 'node(ii,jj,kk)', node(ii,jj,kk)
!      WRITE(6679,*) 'rho(ii,jj,kk)', rho(ii,jj,kk)
!      WRITE(6679,*)
!
!    ELSE
!
!      WRITE(6679,*) '(ii,jj,kk) is out of bounds'
!      WRITE(6679,*) 'ii,jj,kk:', ii,jj,kk
!      WRITE(6679,*) 'imin',imin
!      WRITE(6679,*) 'imax',imax
!      WRITE(6679,*) 'jmin',jmin
!      WRITE(6679,*) 'jmax',jmax
!      WRITE(6679,*) 'kmin',kmin
!      WRITE(6679,*) 'kmax',kmax
!      WRITE(6679,*) 'node(ii,jj,kk)', node(ii,jj,kk)
!      WRITE(6679,*) 'rho(ii,jj,kk)', rho(ii,jj,kk)
!      WRITE(6679,*)
!
!    END IF
!      
!  END DO
!  CLOSE(6679)

END IF


! velocity and scalar (use boundary conditions)
u(i,j,k) 	= ubx									! wall velocity			
v(i,j,k) 	= uby														
w(i,j,k) 	= ubz
phi(i,j,k)	= phiWall							! scalar			

!rho(i,j,k) = denL
!! velocity and scalar (use boundary conditions)
!u(i,j,k) 	= 0									! wall velocity			
!v(i,j,k) 	= 0														
!w(i,j,k) 	= 0



! distribution functions (set to equilibrium)
DO m=0,NumDistDirs
  CALL Equilibrium_LOCAL(m,rho(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),feq)	! distribution functions
  f(m,i,j,k) = feq
END DO

!------------------------------------------------
END SUBROUTINE SetProperties
!------------------------------------------------

!================================================
END MODULE Geometry
!================================================
