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
SUBROUTINE AdvanceGeometry		! advances the geometry in time
!==================================================================================================
IMPLICIT NONE 

CALL BoundaryPosition  		         	! Calculate the radius at the current time step
CALL BoundaryVelocity			          ! Calculate the velocity at boundary point
CALL SetNodes			                	! Flag the fluid/solid nodes based on the new geometry
!==================================================================================================
END SUBROUTINE AdvanceGeometry
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

REAL(dbl) :: v0(0:nz+1), v1(0:nz+1), v2(0:nz+1)   ! velocity arrays for each mode
REAL(dbl) :: lambdaC                              ! wavelength of the cos segments (mode 2)
REAL(dbl) :: time
REAL(dbl) :: A_1, A_2, AA
REAL(dbl) :: alpha, beta, CC, DD, EE
REAL(dbl) :: A_Base, A_Change 
INTEGER(lng) :: i,j,ii				

!----- Initialize Variables
time  = 0.0_dbl				! time
velDom= 0.0_dbl				! summed velocity
v0    = 0.0_dbl       ! Mode 0 Velocity (Couette)
v1    = 0.0_dbl				! mode 1 velocity
v2    = 0.0_dbl				! mode 2 velocity	
			
time  = iter*tcf			! Current Physical Time
!------------------- Mode 0 -Couette ----------------------------------------
DO i= 0,nz+1  
   v0(i) = s1
END DO

!------------------- Mode 1 - peristalsis -----------------------------------
DO i= 0,nz-1  
   v1(i) = kw1*s1*amp1*(SIN(PI + kw1*(zz(i)+s_MovingF*time-s1*time)))
END DO
v1(nz)=   v1(0)
v1(nz+1)= v1(1)

!------------------- Mode 2 - segmental contractions  -----------------------

!----- Calculate the wall velocity for the first wave
A_base  = 204.13981e-6
A_change= 163.31300e-6
DO i= 0,nlambda2+1
   A_1  = A_Base+ A_Change*(COS(PI+(2*PI*zz(i)/L)))
   A_2  = A_Base+ A_Change*(COS(PI+(2*PI*zz(i)/L)+PI))
   alpha= time/Ts
   beta = (cos(alpha*PI))**2
   CC   = (-2.0_dbl*PI/Ts) * sin(alpha*PI) * cos(alpha*PI)
   DD   = (A_1-A_2)/PI
   EE   = ((beta*A_1+ (1-beta)*A_2)/PI)**(-0.5_dbl)
   v2(i)= CC*DD*EE
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
IF (Flag_Couette) THEN    
   DO i=0,nz+1
      velDom(i)= v0(i)  
   END DO
ELSE 
   DO i=0,nz+1
     velDom(i)= wc1*v1(i) + wc2*v2(i)
   END DO
ENDIF 


!----------------------------------------------------------------------------
!----- Fill out the local velocity array
vel(0:nzSub+1) = velDom(kMin-1:kMax+1)/vcf

!==================================================================================================
END SUBROUTINE BoundaryVelocity
!==================================================================================================








!==================================================================================================
SUBROUTINE SetNodes					    ! defines the geometry via "node" array of flags
!==================================================================================================
IMPLICIT NONE 

INTEGER(lng) :: i,j,k,m,iComm   ! index variables
REAL(dbl)    :: rijk            ! radius of the current node
REAL(dbl)    :: ubx,uby,ubz     ! boundary velocity
INTEGER(lng) :: mpierr          ! MPI standard error variable 

r(0:nzSub+1) = rDom(kMin-1:kMax+1)  ! Fill out the local radius array

!----- Flag the interior nodes and give values to nodes that just came in
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub

         IF (Flag_Couette) THEN !----- Couette Geometry ---------------------
            IF (abs(x(i)).LT.r(k)) THEN
               node(i,j,k)= FLUID 
            ELSE
               node(i,j,k)= SOLID  
            END IF
         ELSE !----------------------- Intesitne Geometry -------------------
            rijk = SQRT(x(i)*x(i) + y(j)*y(j))
            !write(*,*) 'rijk,r',rijk,r(k)
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
         END IF

      END DO
   END DO
END DO

!----- Loop through the phantom nodes, and set the entity, but do not give values
!----- YZ Faces
DO iComm=1,2
   i = YZ_RecvIndex(OppCommDir(iComm))			! i index of the phantom nodes
   DO j=0,nySub+1_lng
      DO k=0,nzSub+1_lng
         IF (abs(x(i)).LT.r(k)) THEN
            node(i,j,k)= FLUID 
         ELSE
            node(i,j,k)= SOLID  
         END IF
      END DO
   END DO
END DO

!------ ZX Faces
DO iComm=3,4
   j = ZX_RecvIndex(OppCommDir(iComm))			! j index of the phantom nodes
   DO i=0,nxSub+1_lng
      DO k=0,nzSub+1_lng
         IF (abs(x(i)).LT.r(k)) THEN
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
         IF (abs(x(i)).LT. r(k)) THEN
            node(i,j,k) = FLUID				! set the SOLID node that just came in to FLUID
         ELSE
            node(i,j,k) = SOLID				! if rijk is GT r(k) then it's a SOLID node
         END IF
      END DO
   END DO
END DO

!----- Balaji added to make domain full 3D
IF (domaintype .EQ. 0) THEN          ! only needed when planes of symmetry exist
   CALL SymmetryBC                   ! ensure symmetric node placement
ENDIF
!==================================================================================================
END SUBROUTINE SetNodes
!==================================================================================================









!==================================================================================================
SUBROUTINE SetProperties(i,j,k,ubx,uby,ubz)	! set properties for nodes that just came into fluid domain (uncovered)
!==================================================================================================
IMPLICIT NONE 

INTEGER(lng), INTENT(IN) :: i,j,k           ! current node location
REAL(dbl), INTENT(IN)    :: ubx,uby,ubz     ! velocity of the boundary
INTEGER(lng)             :: m,ii,jj,kk      ! index variables
INTEGER(lng)             :: numFLUIDs       ! number of fluid nodes
REAL(dbl)                :: rhoSum, rhoTemp ! sum of the densities of the neighboring fluid nodes, pre-set density
REAL(dbl)                :: feq             ! equilibrium distribution function
CHARACTER(7)             :: iter_char       ! iteration stored as a character

!----- initialize the sum of surrounding densities
rhoSum   = 0.0_dbl
numFLUIDs= 0_lng

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
