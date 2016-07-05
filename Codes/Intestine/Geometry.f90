!==================================================================================================
MODULE Geometry						  ! Defines the geometry for the simulation
!==================================================================================================
USE SetPrecision      
USE Setup
USE LBM
USE IC
USE BClbm
USE MPI

IMPLICIT NONE 

CONTAINS



!==================================================================================================
SUBROUTINE Geometry_Setup				! sets up the geometry
!==================================================================================================
IMPLICIT NONE

INTEGER :: isize,idate(8)				! size of seed array for the random number genreator, array for output of DATE_AND_TIME
INTEGER,ALLOCATABLE  :: iseed(:)			! seeds for random number generator
INTEGER(lng) :: i,j,k,kk,iCon,it,iPer,nPers_INT		! index variables
INTEGER(lng) :: nvz,nvt,n,g				! index variables
INTEGER(lng) :: mpierr					! MPI standard error variable 
REAL(dbl)    :: macroFreq				! macroscopic contraction frequency
INTEGER(lng) :: xaxis,yaxis				! axes index variables

!----- Define the lattice <=> physical conversion factors
IF (domaintype .EQ. 0) THEN
   xcf	= (0.5_lng*D)/(nx-1_lng)			! length conversion factor: x-direction
   ycf	= (0.5_lng*D)/(ny-1_lng)			! length conversion factor: y-direction
ELSE
   xcf	= (1.0_lng*D)/(nx-1_lng)			! length conversion factor: x-direction
   ycf	= (1.0_lng*D)/(ny-1_lng)			! length conversion factor: y-direction
ENDIF

zcf 	= L/nz						! length conversion factor: z-direction
tcf 	= nuL*((xcf*xcf)/nu)				! time conversion factor
dcf 	= den/denL					! density conversion factor
vcf 	= xcf/tcf					! velocity conversion factor
pcf 	= cs*cs*vcf*vcf					! pressure conversion factor


!----- Determine the number of time steps to run
IF (restart) THEN
   OPEN(55,FILE='Restart-iter.dat')				! open initial iteration file
   READ(55,*) iter0					! read and set initial iteration
   CLOSE(55)
   iter = iter0-1_lng					! set the in
   nt = ANINT((nPers*Tmix)/tcf) + iter
ELSE
   nt = ANINT((nPers*Tmix)/tcf)
END IF


!----- Initialize arrays
node	= -99_lng					! node flag array
rDom	= 0.0_dbl					! radius at each z-location
r	= 0.0_dbl					! temporary radius array for entire computational domain
velDom	= 0.0_dbl					! wall velocity at each z-location (global)
vel	= 0.0_dbl					! wall velocity at each z-location (local)

!----- Check to ensure xcf=ycf=zcf (LBM grid must be cubic)
IF ((ABS(xcf-ycf) .GE. 1E-8) .OR. (ABS(xcf-zcf) .GE. 1E-8) .OR. (ABS(ycf-zcf) .GE. 1E-8)) THEN
   OPEN(1000,FILE="error.txt")
   WRITE(1000,*) "Conversion factors not equal... Geometry_Setup.f90: Line 93."
   WRITE(1000,*) "xcf=", xcf, "ycf=", ycf, "zcf=", zcf
   WRITE(1000,*) "L=", L, "D=", D
   WRITE(1000,*) "nx=", nx, "ny=", ny, "nz=", nz
   CLOSE(1000)
   STOP
END IF

!------ IF CONDITION TO CHECK IF THE DOMAIN TO BE MODELLED IS FULL CYLINDER OR JUST A QUARTER OF A CYLINDER
IF (domaintype .EQ. 0) THEN 
   !----- Fill out x,y,z arrays (local)
   DO i=0,nxSub+1
      x(i) = ((iMin - 1_lng) + (i-1_lng))*xcf
   END DO
   DO j=0,nySub+1
      y(j) = ((jMin - 1_lng) + (j-1_lng))*ycf
   END DO
   DO k=0,nzSub+1
      z(k) = (((kMin - 1_lng) + k) - 0.5_dbl)*zcf
   END DO
      
   !------ Fill out xx,yy,zz arrays (global)
   DO i=0,nx+1
      xx(i) = (i-1_lng)*xcf
   END DO
   DO j=0,ny+1
      yy(j) = (j-1_lng)*ycf
   END DO
   DO k=0,nz+1
      zz(k) = (k - 0.5_dbl)*zcf
   END DO
      
   !----- Center node locations
   Ci = 1	
   Cj = 1
   Ck = ANINT(0.5_dbl*nz)
ELSE
   !----- begin Balaji added 
   xaxis=ANINT(0.5_dbl*(nx+1))
   yaxis=ANINT(0.5_dbl*(ny+1))
      
   !----- Fill out x,y,z arrays (local)
   DO i=0,nxSub+1
      x(i) = ((iMin - 1_lng - (xaxis-1_lng)) + (i-1_lng))*xcf
   END DO
   DO j=0,nySub+1
      y(j) = ((jMin - 1_lng - (yaxis-1_lng)) + (j-1_lng))*ycf
   END DO
   DO k=0,nzSub+1
      z(k) = (((kMin - 1_lng) + k) - 0.5_dbl)*zcf
   END DO
      
   !----- Fill out xx,yy,zz arrays (global)
   DO i=0,nx+1
      xx(i) = (i-1_lng-(xaxis-1_lng))*xcf
   END DO
   DO j=0,ny+1
      yy(j) = (j-1_lng-(yaxis-1_lng))*ycf
   END DO
   DO k=0,nz+1
      zz(k) = (k - 0.5_dbl)*zcf
   END DO
      
   !----- Center node locations
   Ci = xaxis
   Cj = yaxis
   Ck = ANINT(0.5_dbl*nz)
   !------ end Balaji added 
ENDIF

!----- Mode 1 - Peristalsis ------------------------------------------------------------------------
a1	 = (0.5_dbl*D)/(2.0_dbl - epsOVERa1)					! mean half-width of wave1
eps1 	 = epsOVERa1*a1								! occlusional distance
lambda1	 = L/numw1								! wavelength
aOVERlam1= a1/lambda1								! ratio of mean half-width to wavelength 
kw1	 = (2.0_dbl*PI)/lambda1							! wave number
amp1	 = 0.5_dbl*((0.5_dbl*D)-eps1)						! amplitude of the wave
Tp	 = lambda1/s1								! peristaltic period
Re1	 = ((s1*(0.5_dbl*D))/nu)*((0.5_dbl*D)/lambda1)				! Reynolds number based on mode 1

!----- Mode 2 - Segmental Contractions -------------------------------------------------------------
a2	 = (0.5_dbl*D)/(2.0_dbl - epsOVERa2)					! mean half-width of wave1 (based on peristalsis definition)
eps2 	 = epsOVERa2*a2								! occlusional distance
lambda2	 = L/numw2								! wavelength (physical units)
nlambda2 = nz/numw2								! wavelength (nodes)
aOVERlam2= a2/lambda2								! ratio of mean half-width to wavelength 
amp2	 = 0.5_dbl*((0.5_dbl*D)-eps2)						! amplitude of the wave
shift2	 = 0.5_dbl*((0.5_dbl*D)+eps2)						! amplitude of the wave
segment	 = nlambda2/6_lng							! length of each segment of the segmental wave   !!!!! CAREFUL HERE WITH SYMMETRY!
seg1L	 = 1_lng + segment							! left point of sloped segement 1
seg1R	 = 1_lng + 2_lng*segment							! right point of sloped segement 1
seg2R	 = nlambda2 - segment							! right point of sloped segement 2
seg2L	 = nlambda2 - (2_lng*segment)						! left point of sloped segement 2
s2	 = (0.5_dbl*D)/Ts							! speed of collapse fo segmental contraction
Re2	= (s2*(0.5_dbl*D))/nu							! Reynolds number based on mode 2

IF (restart .EQV. .FALSE.) THEN
   CALL AdvanceGeometry
END IF
!==================================================================================================
END SUBROUTINE Geometry_Setup
!==================================================================================================








!==================================================================================================
SUBROUTINE AdvanceGeometry		! advances the geometry in time
!==================================================================================================
IMPLICIT NONE 

CALL BoundaryPosition  			! Calculate the radius at the current time step
CALL BoundaryVelocity			! Calculate the velocity at boundary point
CALL SetNodes				! Flag the fluid/solid nodes based on the new geometry
!==================================================================================================
END SUBROUTINE AdvanceGeometry
!==================================================================================================









!==================================================================================================
SUBROUTINE BoundaryPosition		! Calculates the position of the wall at the current time step
!==================================================================================================
IMPLICIT NONE

REAL(dbl) :: h1(0:nz+1)			! Mode 1 (peristalsis)
REAL(dbl) :: h2(0:nz+1)			! Mode 2	(segmental)
REAL(dbl) :: Ac, lambdaC, shiftC	! temporary variables for the cos slopes
REAL(dbl) :: time			! time
INTEGER(lng) :: i,j,ii,k		! indices

!----- Initialize Variables
time = 0.0_dbl				! time					
h1   = 0.0_dbl				! mode 1 height
h2   = 0.0_dbl				! mode 2 height
rDom = 0.0_dbl				! summed height

!----- Current Physical Time
time = iter*tcf

!------------------------- Mode 1 - peristalsis -----------------------------
DO i=0,nz-1
   h1(i) = amp1*( COS(PI + kw1*(zz(i)-s1*time)) ) + 0.5_dbl*D-amp1
   !Yanxing's expression
   !h1(i)= amp1*sin(2.0_dbl*PI*((real(i,dbl)-0.5_dbl)/real(nz,dbl)-0.1_dbl*iter/real(nz,dbl))+pi/2.0_dbl)+ (0.5_dbl*D - amp1)
END DO

!------ since PI cannot be stored exactly, the wavelength(s) does/do not EXACTLY span the domain...
!------ set h1(nz) to h1(0) and h1(nz+1) to h(1) to ensure periodicity
h1(nz) 	= h1(0)
h1(nz+1)= h1(1)


!------------------- Mode 2 - segmental contractions ------------------------
!----- Calculate the geometry for the first wave

!----- First Straight Piece
DO i=0,seg1L
   h2(i) = amp2*(COS(PI+((2.0_dbl*PI)/Ts)*time)) + shift2
END DO

!----- Second Straight Piece
DO i=seg1R,seg2L
   h2(i) = amp2*(COS(PI+((2.0_dbl*PI)/Ts)*(time-(Ts/2.0_dbl)))) + shift2
END DO

!----- Third Straight Piece
DO i=seg2R,nlambda2+1
   h2(i) = amp2*(COS(PI+((2.0_dbl*PI)/Ts)*time)) + shift2
END DO

!----- First Cos Piece
Ac	= 0.5_dbl*(h2(seg1L)-h2(seg1R))
lambdaC	= 2.0_dbl*(zz(seg1L)-zz(seg1R))
shiftC	= 0.5_dbl*(h2(seg1L)+h2(seg1R))

DO i=seg1L+1,seg1R-1
   h2(i) = Ac*COS((2.0_dbl*PI/lambdaC)*(zz(i)-zz(seg1L))) + shiftC
END DO

!----- Second Cos Piece
Ac	= 0.5_dbl*(h2(seg2L)-h2(seg2R))
lambdaC	= 2.0_dbl*(zz(seg2L)-zz(seg2R))
shiftC	= 0.5_dbl*(h2(seg2L)+h2(seg2R))
DO i=seg2L+1,seg2R-1
   h2(i) = Ac*COS((2.0_dbl*PI/lambdaC)*(zz(i)-zz(seg2L))) + shiftC
END DO

!---- Repeat for the rest of the waves
DO j=1,(numw2-1)
   DO i=0,nlambda2+1
      ii = i + j*nlambda2
      h2(ii) = h2(i)
   END DO
END DO

!"fudging" to make sure that the whole domain is filled (and periodic) - more logic (and computational expense would be
!necessary to do this correctly: ideally, one would determine if an even or odd number of waves was specified
!and then work from either end, and meet in the middle to ensure a symetric domain...
!h2(nz-1:nz+1) = h2(1)

!-------------------------------- Mode Sum  ---------------------------------
!----- Sum the modes in a weighted linear combination
DO i=0,nz+1
   rDom(i) = wc1*h1(i) + wc2*h2(i)
END DO

!----- Fill out the local radius array
r(0:nzSub+1) = rDom(kMin-1:kMax+1)

!IF (myid .EQ. master) THEN
!   CALL SurfaceArea				!----- calculate the surface area
!END IF

!==================================================================================================
END SUBROUTINE BoundaryPosition
!==================================================================================================








!==================================================================================================
SUBROUTINE SurfaceArea			! calculate the surface area at the current time and write it to a file
!==================================================================================================
IMPLICIT NONE

REAL(dbl) :: SA				! surface area
REAL(dbl) :: r2,r1,z2,z1		! radius and z-location for each set of consecutive points
INTEGER(lng) :: kk			! index variable

!----- initialize the surface area
SA = 0.0_dbl

!----- approximate the surface area as if the nodes were connected linearly with the neighboring nodes
!----- surface area between left phantom node and 1st domain node
r1 = 0.5_dbl*(rDom(0) + rDom(1))
z1 = 0.5_dbl*(zz(0) + zz(1))
r2 = rDom(1)
z2 = zz(1)
SA = -PI*(r1 + r2)*SQRT(1.0_dbl + ((r1-r2)/(z1-z2))**2)*(z1 - z2)

!----- surface area between right phantom node and last domain node
r1 = rDom(nz)
z1 = zz(nz)
r2 = 0.5_dbl*(rDom(nz) + rDom(nz+1))
z2 = 0.5_dbl*(zz(nz) + zz(nz+1))
SA = SA - PI*(r1 + r2)*SQRT(1.0_dbl + ((r1-r2)/(z1-z2))**2)*(z1 - z2)

!----- interior domain nodes
DO kk=1,nz-1
   r1 = rDom(kk)
   r2 = rDom(kk+1)
   z1 = zz(kk)
   z2 = zz(kk+1)
   SA = SA - PI*(r1 + r2)*SQRT(1.0_dbl + ((r1-r2)/(z1-z2))**2)*(z1 - z2)
END DO

!----- account for the villi
SA = SA - numVilliActual*(PI*Rv*Rv)					! subtract the cross sectional area of the villous bases from the total outer surface area
SA = SA + numVilliActual*(2.0_dbl*PI*Rv*(Lv-Rv))			! add the surface area from the villous cylinders
SA = SA + numVilliActual*(2.0_dbl*PI*Rv*Rv)				! add the surface area from the villous tips

!----- open and write to a file
!IF (iter .GT. 0) THEN
!   WRITE(2474,'(2E25.15)') REAL(iter/(nt/nPers)), SA			! write surface area to file
!   CALL FLUSH(2474)
!END IF

!==================================================================================================
END SUBROUTINE SurfaceArea
!==================================================================================================





!==================================================================================================
SUBROUTINE BoundaryVelocity			! defines the velocity of the solid boundaries (fills "ub", "vb", and "wb" arrays)
!==================================================================================================
IMPLICIT NONE 

REAL(dbl) :: v1(0:nz+1), v2(0:nz+1)		! velocity arrays for each mode
REAL(dbl) :: lambdaC				! wavelength of the cos segments (mode 2)
REAL(dbl) :: time				! time
INTEGER(lng) :: i,j,ii				! indices

!----- Initialize Variables
time  = 0.0_dbl				! time
velDom= 0.0_dbl				! summed velocity
v1    = 0.0_dbl				! mode 1 velocity
v2    = 0.0_dbl				! mode 2 velocity	
			
time  = iter*tcf			! Current Physical Time

!------------------------- Mode 1 - peristalsis -----------------------------
DO i= 0,nz-1 ! Balaji added to ensure periodicity just like in h1. 
   v1(i) = kw1*s1*amp1*(SIN(PI + kw1*(zz(i) - (s1*time))))
   !v1(i)= -kw1*s1*amp1*(SIN(kw1*(zz(i) - (s1*time))))
   !Yanxing's expression
   !v1(i)= -kw1*s1*amp1*cos(2.0_dbl*PI*((real(i,dbl)-0.5_dbl)/real(nz,dbl)-0.1_dbl*iter/real(nz,dbl))+pi/2.0_dbl)
END DO
v1(nz)=v1(0)
v1(nz+1)=v1(1)

!------------------- Mode 2 - segmental contractions  -----------------------

!----- Calculate the wall velocity for the first wave
!----- First Straight Piece
DO i=0,seg1L
   v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*time))*((2.0_dbl*PI)/Ts)
END DO

!----- Second Straight Piece
DO i=seg1R,seg2L
   v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*(time-(Ts/2.0_dbl))))*((2.0_dbl*PI)/Ts)
END DO

!----- Third Straight Piece
DO i=seg2R,nlambda2
   v2(i) = -amp2*(SIN(((2.0_dbl*PI)/Ts)*time))*((2.0_dbl*PI)/Ts)
END DO

!----- First Cos Piece
lambdaC	= 2.0_dbl*(zz(seg1L)-zz(seg1R))
DO i=seg1L+1,seg1R-1
   v2(i) = (0.5_dbl*(v2(seg1L)-v2(seg1R)))*COS((2.0_dbl*PI/lambdaC)*(zz(i)-zz(seg1L))) &
        + (0.5_dbl*(v2(seg1L)+v2(seg1R)))
END DO

!----- Second Cos Piece
lambdaC	= 2.0_dbl*(zz(seg2L)-zz(seg2R))
DO i=seg2L+1,seg2R-1
   v2(i) = (0.5_dbl*(v2(seg2L)-v2(seg2R)))*COS((2.0_dbl*PI/lambdaC)*(zz(i)-zz(seg2L))) &
        + (0.5_dbl*(v2(seg2L)+v2(seg2R)))
END DO

!----- Repeat for the rest of the waves
DO j=1,(numw2-1)
   DO i=1,nlambda2+1
      ii = i + j*nlambda2
      v2(ii) = v2(i)
   END DO
END DO

!"fudging" to make sure that the whole domain is filled (and periodic) - more logic (and computational expense would be
!necessary to do this correctly: ideally, one would determine if an even or odd number of waves was specified
!and then work from either end, and meet in the middle to ensure a symetric domain...

v2(nz-1:nz+1) = v2(1)

!-------------------------------- Mode Sum  ---------------------------------
!----- Sum the modes in a weighted linear combination
DO i=0,nz+1
   velDom(i) = wc1*v1(i) + wc2*v2(i)
END DO

!----------------------------------------------------------------------------
!----- Fill out the local velocity array
vel(0:nzSub+1) = velDom(kMin-1:kMax+1)/vcf

!==================================================================================================
END SUBROUTINE BoundaryVelocity
!==================================================================================================








!==================================================================================================
SUBROUTINE SetNodes					! defines the geometry via "node" array of flags
!==================================================================================================
IMPLICIT NONE 

INTEGER(lng)	:: i,j,k,m,iComm			! index variables
REAL(dbl)	:: rijk					! radius of the current node
REAL(dbl)      	:: ubx,uby,ubz				! boundary velocity
INTEGER(lng) 	:: mpierr				! MPI standard error variable 

!----- Flag the interior nodes and give values to nodes that just came in
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         rijk = SQRT(x(i)*x(i) + y(j)*y(j))
         IF (rijk .LT. r(k)) THEN
            IF (node(i,j,k) .EQ. SOLID) THEN		! just came into the domain
            !----- calculate the wall velocity (boundary)
               ubx = vel(k)*(x(i)/rijk)
               uby = vel(k)*(y(j)/rijk)
               ubz = 0.0_dbl
              CALL SetProperties(i,j,k,ubx,uby,ubz)
            END IF
            node(i,j,k)	= FLUID				! reset the SOLID node that just came in to FLUID
         ELSE
           node(i,j,k) = SOLID				! if rijk is GT r(k) then it's a SOLID node
         END IF
      END DO
   END DO
END DO

!----- Loop through the phantom nodes, and set the entity, but do not give values
!----- YZ Faces
DO iComm=1,2
   i = YZ_RecvIndex(OppCommDir(iComm))			! i index of the phantom nodes
   DO j=0,nySub+1_lng
      rijk = SQRT(x(i)*x(i) + y(j)*y(j))
      DO k=0,nzSub+1_lng
         IF (rijk .LT. r(k)) THEN
            node(i,j,k) = FLUID				! set the SOLID node that just came in to FLUID
         ELSE
            node(i,j,k) = SOLID				! if rijk is GT r(k) then it's a SOLID node
         END IF
      END DO
   END DO
END DO

!------ ZX Faces
DO iComm=3,4
   j = ZX_RecvIndex(OppCommDir(iComm))			! j index of the phantom nodes
   DO i=0,nxSub+1_lng
      rijk = SQRT(x(i)*x(i) + y(j)*y(j))
      DO k=0,nzSub+1_lng
         IF (rijk .LT. r(k)) THEN
            node(i,j,k) = FLUID				! set the SOLID node that just came in to FLUID
         ELSE
            node(i,j,k) = SOLID				! if rijk is GT r(k) then it's a SOLID node
         END IF
      END DO
   END DO
END DO

!----- XY Faces
DO iComm=5,6
   k = XY_RecvIndex(OppCommDir(iComm))			! k index of the phantom nodes
   DO j=0,nySub+1_lng
      DO i=0,nxSub+1_lng
         rijk = SQRT(x(i)*x(i) + y(j)*y(j))
         IF (rijk .LT. r(k)) THEN
            node(i,j,k) = FLUID				! set the SOLID node that just came in to FLUID
         ELSE
            node(i,j,k) = SOLID				! if rijk is GT r(k) then it's a SOLID node
         END IF
      END DO
   END DO
END DO

!----- Balaji added to make domain full 3D
IF (domaintype .EQ. 0) THEN  				! only needed when planes of symmetry exist
   CALL SymmetryBC					! ensure symmetric node placement
ENDIF
!==================================================================================================
END SUBROUTINE SetNodes
!==================================================================================================









!==================================================================================================
SUBROUTINE SetProperties(i,j,k,ubx,uby,ubz)	! set properties for nodes that just came into fluid domain (uncovered)
!==================================================================================================
IMPLICIT NONE 

INTEGER(lng), INTENT(IN) :: i,j,k				! current node location
REAL(dbl), INTENT(IN)    :: ubx,uby,ubz				! velocity of the boundary
INTEGER(lng)	:: m,ii,jj,kk					! index variables
INTEGER(lng)	:: numFLUIDs					! number of fluid nodes
REAL(dbl)	:: rhoSum, rhoTemp				! sum of the densities of the neighboring fluid nodes, pre-set density
REAL(dbl)	:: feq						! equilibrium distribution function
CHARACTER(7)	:: iter_char					! iteration stored as a character

!----- initialize the sum of surrounding densities
rhoSum = 0.0_dbl
numFLUIDs = 0_lng

!----- calculate the average density of the current node's neighbors
DO m= 1,NumDistDirs
   ii = i + ex(m)
   jj = j + ey(m)
   kk = k + ez(m)
   IF (((ii .GE. 0) .AND. (ii .LE. nxSub+1_lng)) .AND.	&
      ((jj .GE. 0) .AND. (jj .LE. nySub+1_lng)) .AND.	&
      ((kk .GE. 0) .AND. (kk .LE. nzSub+1_lng))) THEN
      IF (node(ii,jj,kk) .EQ. FLUID) THEN
         rhoSum = rhoSum + rho(ii,jj,kk)
         numFLUIDs = numFLUIDs + 1_lng     
      END IF       
  END IF
END DO

!----- This should rarely happen...
IF (numFLUIDs .NE. 0_lng) THEN
   rho(i,j,k) = rhoSum/numFLUIDs
ELSE
   rho(i,j,k) = denL
END IF

!----- velocity and scalar (use boundary conditions)
!rho(i,j,k)= denL							! Density 
u(i,j,k)  = ubx								! wall velocity			
v(i,j,k)  = uby														
w(i,j,k)  = ubz
phi(i,j,k)= phiWall							! scalar			

!----- distribution functions (set to equilibrium)
DO m= 0,NumDistDirs
   CALL Equilibrium_LOCAL(m,rho(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),feq)	! distribution functions
   f(m,i,j,k) = feq
END DO

!==================================================================================================
END SUBROUTINE SetProperties
!==================================================================================================









!==================================================================================================
END MODULE Geometry
!==================================================================================================
