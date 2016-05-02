!==================================================================================================
MODULE BCscalar		! Sets Scalar Boundary Conditions
!==================================================================================================
USE SetPrecision
USE Setup  
USE MPI
USE IC
USE BClbm

CONTAINS


!--------------------------------------------------------------------------------------------------
SUBROUTINE BC_Scalar(m,i,j,k,im1,jm1,km1,phiBC)				! implements the scalar BCs 
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1				! index variables
REAL(dbl),    INTENT(OUT):: phiBC     					! scalar contribution from the boundary condition
INTEGER(lng) :: mm,ip1,jp1,kp1,iB,jB,kB    				! First neighboring node location
REAL(dbl)    :: rhoAstar, phiAstar, feq_Astar,  PkAstar	 		! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl)    :: rhoBstar, phiBstar, fPlusBstar, PkBstar 		! Values interpolated to Bstar location
REAL(dbl)    :: cosTheta, sinTheta					! COS(theta), SIN(theta)
REAL(dbl)    :: ub, vb, wb						! wall velocity (x-, y-, z- components)
REAL(dbl)    :: rijk 							! radius of the solid node
REAL(dbl)    :: Geom_norm_x,Geom_norm_y,Geom_norm_z
REAL(dbl)    :: q, n_prod, n_prod_max

!===========================================================================
! HELP: How to set different boundary conditions
! BC Scalar-zero:       coeffPhi=1      coeffGrad= 0    coeffConst= 0
! BC Scalar-Non-zero:   coeffPhi=1      coeffGrad= 0    coeffConst= phi_BC
! BC-Flux-zero:         coeffPhi=0      coeffGrad= 1    coeffConst= 0
! BC-Flux-Non-zero:     coeffPhi=0      coeffGrad= 1    coeffConst= dphi/dn
! BC-Permeability:      coeffPhi=Pw/Dm  coeffGrad=-1    coeffConst= 0
!===========================================================================

CALL qCalcFarhad(i,q)		

IF (rijk .GE. rOut(k)) THEN
   ub = velOut(km1) 			
   vb = 0.0_dbl 				
   wb = 0.0_dbl 							
ELSE IF (rijk .LE. rIn(k)) THEN
   ub = velIn(km1) 				
   vb = 0.0_dbl					
   wb = 0.0_dbl 					
END IF		

!---------------------------------------------------------------------------------------------------
!----- Computing phi at the wall in case of Dirichlet BC -------------------------------------------
!---------------------------------------------------------------------------------------------------
IF (coeffGrad .EQ. 0.0) then   
   phiWall= coeffConst/coeffPhi
END IF      

!---------------------------------------------------------------------------------------------------
!----- Estimating  phi at the wall in case of Neumann or Mixed BC ----------------------------------
!----- Two Fluid nodes are needed for extrapolation ------------------------------------------------
!----- 1st node; A, is adjacent to the boundary ----------------------------------------------------
!----- 2nd node; B, is a fluid neighbor of A, in the direction closest to geometry-normal ----------
!---------------------------------------------------------------------------------------------------
IF (coeffGrad .NE. 0.0) then
   Geom_norm_x= 1.0
   Geom_norm_y= 0.0
   Geom_norm_z= 0.0
   n_prod_max= 0.0_dbl
   !----- Finding the mth direction closest to normal vector
   !----- Only one fluid neighboring node is needed for interpolation
   DO mm=0,NumDistDirs
      ip1= i+ ex(mm)
      jp1= j+ ey(mm)
      kp1= k+ ez(mm)
      IF (node(ip1,jp1,kp1) .EQ. FLUID) THEN
         n_prod= abs( Geom_norm_x * (ip1-i) + Geom_norm_y * (jp1-j) + Geom_norm_z * (kp1-k))
         IF (n_prod_max .LT. n_prod) THEN
            n_prod_max= n_prod
            iB= ip1
            jB= jp1
            kB= kp1
         END IF
      END IF
   END DO
   phiWall= ((phiTemp(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) -    	&
             (phiTemp(iB,jB,kB)*q*q/(1.0+2.0*q)) -          	&
             (q*(1+q)/(1+2.0*q))* (coeffConst/coeffGrad)    )  	&
           /(1.0- (q*(1+q)/(1+2.0*q))*(coeffPhi/coeffGrad))
END IF

!----- neighboring node (fluid side) ---------------------------------------------------------------
ip1 = i + ex(m)
jp1 = j + ey(m)
kp1 = k + ez(m)
!------ This rarely happens (both neighboring nodes over a line are solid)
!------  use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF

!----- Computing values at A* & the scalar streamed from A* (Chpter 3 paper) -----------------------
rhoAstar= (rho(i,j,k)-rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)	! Extrapolate density
CALL Equilibrium_LOCAL(m,rhoAstar,ub,vb,wb,feq_Astar)    	! f_eq in mth direction
phiAstar= phiWall						! phi at solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar		! Contribution from A* to B*  

!----- Computing values at B* & the scalar streamed from B* (Chpter 3 paper) -----------------------
rhoBstar=   (1-q)*rho(ip1,jp1,kp1)     + q* rho(i,j,k)
phiBstar=   (1-q)*phiTemp(ip1,jp1,kp1) + q* phiTemp(i,j,k)
fPlusBstar= (1-q)*fplus(m,ip1,jp1,kp1) + q* fplus(m,i,j,k)
PkBstar=    (fplusBstar/rhoBstar - wt(m)*Delta)*phiBstar

!----- Scalar contribution from wall to the node----------------------------------------------------
phiBC=      PkAstar+ (PkAstar- PkBstar)*(1-q)
!------------------------------------------------
END SUBROUTINE BC_Scalar
!------------------------------------------------
















!================================================
END MODULE BCscalar
!================================================
