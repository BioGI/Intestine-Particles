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
REAL(dbl), INTENT(OUT)   :: phiBC     					! scalar contribution from the boundary condition

INTEGER(lng) :: ip1,jp1,kp1 						! neighboring nodes (2 away from the wall)
REAL(dbl)    :: q
REAL(dbl)    :: rhoAstar, phiAstar, feq_Astar,  PkAstar	 		! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl)    :: rhoBstar, phiBstar, fPlusBstar, PkBstar 		! Values interpolated to Bstar location
REAL(dbl)    :: cosTheta, sinTheta					! COS(theta), SIN(theta)
REAL(dbl)    :: ub, vb, wb						! wall velocity (x-, y-, z- components)
REAL(dbl)    :: rijk 							! radius of the solid node

!===========================================================================
! HELP: How to set different boundary conditions
! BC Scalar-zero:       coeffPhi=1      coeffGrad= 0    coeffConst= 0
! BC Scalar-Non-zero:   coeffPhi=1      coeffGrad= 0    coeffConst= phi_BC
! BC-Flux-zero:         coeffPhi=0      coeffGrad= 1    coeffConst= 0
! BC-Flux-Non-zero:     coeffPhi=0      coeffGrad= 1    coeffConst= dphi/dn
! BC-Permeability:      coeffPhi=Pw/Dm  coeffGrad=-1    coeffConst= 0
!===========================================================================

CALL qCalc(m,i,j,k,im1,jm1,km1,q)

IF (rijk .GE. rOut(k)) THEN
   ub = velOut(km1) 			
   vb = 0.0_dbl 				
   wb = 0.0_dbl 							
ELSE IF (rijk .LE. rIn(k)) THEN
   ub = velIn(km1) 				
   vb = 0.0_dbl					
   wb = 0.0_dbl 					
END IF		

!----- neighboring node (fluid side)	
ip1 = i + ex(m) 			
jp1 = j + ey(m)			
kp1 = k + ez(m)		

!------ if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF	

IF (coeffGrad .EQ. 0) then						! Dirichlet BC
   phiWall= coeffConst/coeffPhi
ELSE 			     						! Calculates phiWall for flux BC (eq. 28 in paper)				
   phiWall= ( (phiTemp(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) - 	&
	      (phiTemp(ip1,j,k)*q*q/(1.0+2.0*q)) -		&
	      (q*(1+q)/(1+2.0*q))* (coeffConst/coeffGrad)    ) 	&
            / (1.0- (q*(1+q)/(1+2.0*q))*(coeffPhi/coeffGrad))
END IF

!----- Computing values at A* & scalar streamed from A* (Chpter 3 paper)
rhoAstar= (rho(i,j,k)- rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)	! Extrapolate density
CALL Equilibrium_LOCAL(m,rhoAstar,ub,vb,wb,feq_Astar)    		! Equibrium distribution function in mth direction
phiAstar= phiWall							! phi at solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar			! Contribution from wall in mth direction (0 if phiWall=0)

!------ Computing values at B* & scalar streamed from B* (Chpter 3 paper)
rhoBstar=   (1-q)*rho(ip1,jp1,kp1)     + q* rho(i,j,k)
phiBstar=   (1-q)*phiTemp(ip1,jp1,kp1) + q* phiTemp(i,j,k)
fPlusBstar= (1-q)*fplus(m,ip1,jp1,kp1) + q* fplus(m,i,j,k)
PkBstar=    (fplusBstar/rhoBstar - wt(m)*Delta)*phiBstar

phiBC=      PkAstar+ (PkAstar- PkBstar)*(1-q)
!------------------------------------------------
END SUBROUTINE BC_Scalar
!------------------------------------------------

!================================================
END MODULE BCscalar
!================================================
