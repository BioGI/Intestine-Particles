!===================================================================================================
MODULE PassiveScalar	
!===================================================================================================
USE SetPrecision
USE Setup
USE IC
USE BClbm
USE BCscalar

IMPLICIT NONE 

CONTAINS


!===================================================================================================
SUBROUTINE Scalar_Setup					      ! sets up the passive scalar component
!===================================================================================================
IMPLICIT NONE

phi    = 0.0_dbl							
phiTemp= 0.0_dbl							
phiInNodes = 0.0_dbl
phiOutNodes = 0.0_dbl

Dm     = nuL/Sc							! binary molecular diffusivity (scalar in fluid)
Dmcf   = (zcf*zcf)/tcf						! conversion factor for diffusivity
Delta  = 1.0_dbl- 6.0_dbl*Dm					! scalar diffusion parameter

!---- scalar standard devation for gaussian distributions
sigma  = 0.1_dbl*D						! 1/10 of Diameter

!---- determine scalar starting iteration
phiStart= NINT((phiPer*Tmix)/tcf)
IF (phiPer.EQ.0.0) THEN
   phiStart= NINT((phiPer*Tmix)/tcf)+1 				! Balaji: to add 1 as for phiPer=0, phiSTart=0. But iter never has a value 0.
ENDIF

!===================================================================================================
END SUBROUTINE Scalar_Setup
!===================================================================================================








!===================================================================================================
SUBROUTINE Scalar				  ! calculates the evolution of scalar in the domain
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1				! index variables
REAL(dbl) :: phiBC						! scalar contribution from boundary
REAL(dbl) :: phiOutSurf,phiInSurf				! scalar contribution coming from and going into the boundary
REAL(dbl) :: tausgs						! contribution form tau_sgs term from particle closure
REAL(dbl) :: zcf3						! Cell volume

CALL ScalarDistribution						! sets/maintains initial distributions of scalar [MODULE: ICBC.f90]

!----- store the previous scalar values
phiTemp 	    = phi
Negative_phi_Counter= 0
Negative_phi_Worst  = 0.0_dbl

!----- Stream the scalar
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         IF (node(i,j,k) .EQ. FLUID) THEN
            phi(i,j,k) = Delta*phiTemp(i,j,k)  
	    phi(i,j,k) = phi(i,j,k)+ delphi_particle(i,j,k) 	! Effects of drug release
!           tausgs = ((tausgs_particle_x(i+1,j,k)-tausgs_particle_x(i-1,j,k)) + &
!	              (tausgs_particle_y(i,j+1,k)-tausgs_particle_y(i,j-1,k)) + &
!	     	      (tausgs_particle_z(i,j,k+1)-tausgs_particle_z(i,j,k-1)))*0.5_dbl
!	    phi(i,j,k) = phi(i,j,k)+ tausgs 			.
            DO m=0,NumDistDirs
               !-----  neighboring node --------------------------------------------------------------
               im1= i- ex(m)
               jm1= j- ey(m)
               km1= k- ez(m)
               IF (node(im1,jm1,km1) .EQ. FLUID) THEN 
                  phi(i,j,k) = phi(i,j,k) + (fplus(m,im1,jm1,km1)/rho(im1,jm1,km1) - wt(m)*Delta)*phiTemp(im1,jm1,km1)
               ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN									! macro- boundary
                  CALL BC_Scalar(m,i,j,k,im1,jm1,km1,phiBC) 
                  phi(i,j,k) = phi(i,j,k) + phiBC     
                  CALL AbsorbedScalarS(i,j,k,m,im1,jm1,km1,phiBC) 								! measure the absorption rate
               ELSE IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN				! villi
          
               ELSE
                  OPEN(1000,FILE="error.txt")
                  WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node(im1,jm1,km1) is out of range"
                  WRITE(1000,*) "iter",iter,"m=",m,"i=",i,"j=",j,"k=",k
                  WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
                  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
                  WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
                  WRITE(1000,*) "node(i,j,k)=",node(i,j,k),"node(im1,jm1,km1)=",node(im1,jm1,km1)
                  CLOSE(1000)
                  STOP
               END IF
            END DO

!---------- node volume in physical units (cm^3) so when printing the drung units are "mole
            zcf3 = zcf*zcf*zcf

!---------- Monitoring the negative phi
            IF (phi(i,j,k) .LT. 0.0_dbl) THEN
               Negative_phi_Counter = Negative_phi_Counter +1.0
               Negative_phi_Total   = Negative_phi_Total + phi(i,j,k) * zcf3 
               IF (phi(i,j,k) .LT. Negative_phi_Worst) THEN
                  Negative_phi_Worst = phi(i,j,k)
               ENDIF
	       phi(i,j,k) = 0.0_dbl
            END IF
         END IF
      END DO
   END DO
END DO

IF (Negative_phi_Counter.LT. 1.0) THEN
   Negative_phi_Counter = 1.0
ENDIF 

!----- Monitorin the Negative phi issue
write(2118,*) iter, Negative_phi_Counter, Negative_phi_Total, Negative_phi_Worst, Negative_phi_Total/Negative_phi_Counter

!----- Add the amount of scalar absorbed through the outer surfaces
phiAbsorbed = phiAbsorbedS 						
      
!===================================================================================================
END SUBROUTINE Scalar
!===================================================================================================







!===================================================================================================
SUBROUTINE AbsorbedScalarS(i,j,k,m,im1,jm1,km1,phiBC) 	  ! Monitoring the abosrption at boundaries	
!===================================================================================================
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m,im1,jm1,km1					! index variables
REAL(dbl),    INTENT(IN) :: phiBC   		  				! scalar contribution from the boundary condition

INTEGER(lng) :: ip1,jp1,kp1
REAL(dbl)    :: phiOUT, phiIN							! scalar values exchanged with the wall
REAL(dbl)    :: phiAbsorbedSleft,phiAbsorbedSright
REAL(dbl)    :: phiINleft, phiINright, phiOUTright,phiOUTleft
REAL(dbl)    :: feq_AO_u0
REAL(dbl)    :: rhoAstar,phiAstar, PkAstar,feq_Astar,feq_Bstar
REAL(dbl)    :: rhoA, PkA, feq_A
REAL(dbl)    :: fPlusBstar, rhoBstar, phiBstar, PkBstar
REAL(dbl)    :: ub,vb,wb, ubb,vbb,wbb
REAL(dbl)    :: q

!CALL qCalc_iter(m,i,j,k,im1,jm1,km1,q)
!ubb= 0.0_dbl
!vbb= 0.0_dbl
!wbb= 0.0_dbl
!---------------------------------------------------------------------------------------------------
!----- Computing phiOUT ----------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!CALL Equilibrium_LOCAL(bb(m),rho(i,j,k),ubb,vbb,wbb,feq_AO_u0)
!phiOUT= (feq_AO_u0/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)
!---------------------------------------------------------------------------------------------------
!---- Conmputing phiIN------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!----- neighboring node (fluid side)	
!ip1 = i + ex(m) 			
!jp1 = j + ey(m)			
!kp1 = k + ez(m)		
!IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
!  ip1 = i
!  jp1 = j
!  kp1 = k
!END IF	
!----- Computing values at A* & scalar streamed from A* (Chpter 3 paper)
!rhoAstar= (rho(i,j,k)- rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)	! extrapolate the density
!CALL Equilibrium_LOCAL(m,rhoAstar,ubb,vbb,wbb,feq_Astar)		! calculate the equibrium distribution function in the mth direction
!phiAstar= phiWall							! getting phi at the solid surface
!PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar			! contribution from the wall in mth direction (0 if phiWall=0)
!------ Computing values at B* & scalar streamed from B* (Chpter 3 paper)
!rhoBstar=   (1-q)*rho(ip1,jp1,kp1)     + q*rho(i,j,k)
!CALL Equilibrium_LOCAL(m,rhoBstar,ubb,vbb,wbb,feq_Bstar)
!phiBstar=   (1-q)*phiTemp(ip1,jp1,kp1) + q*phiTemp(i,j,k)
!PkBstar=    (feq_Bstar/rhoBstar - wt(m)*Delta)*phiBstar
!
!phiIN= PkAstar+ (PkAstar- PkBstar)*(1-q)
!!---- Modification for moving boundary in case of using only A and A* for BC
!rhoA= rho(i,j,k)
!CALL Equilibrium_LOCAL(m,rhoA,ubb,vbb,wbb,feq_A) 
!PkA= (feq_A/rhoA - wt(m)*Delta)*phiTemp(i,j,k) 
!IF(q .LT. 0.25) THEN
!  q = 0.25_dbl
!END IF 
!phiIN   = ((PkAstar - PkA)/q) + PkAstar

!--- No Modifications for moving boundaries
phiIN= phiBC                                                    	 ! contribution from wall to crrent node (in)
phiOUT= (fplus(bb(m),i,j,k)/rho(i,j,k)-wt(bb(m))*Delta)*phiTemp(i,j,k)
phiAbsorbedS = phiAbsorbedS + (phiOUT-phiIN)				! scalar absorbed at current location in mth direction
!===================================================================================================
END SUBROUTINE AbsorbedScalarS
!===================================================================================================









!===================================================================================================
SUBROUTINE ScalarInOut	   ! measures scalar left or entered the domain through the inlet or outlet
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1,iComm			! index variables
REAL(dbl) :: phiOUT, phiIN					! scalar values exchanged with the wall

!------ XY Faces (z-direction)
DO iComm=5,6
   IF (SubID(iComm) .EQ. 0) THEN					
      k = XY_SendIndex(iComm)					! k index
      DO j=1,nySub
         DO i=1,nxSub
            IF ((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN
               DO m=1,NumFs_face
                  !----- i,j,k location of neighboring node
                  im1 = i - ex(bb(f_Comps(iComm,m)))
                  jm1 = j - ey(bb(f_Comps(iComm,m)))
                  km1 = k - ez(bb(f_Comps(iComm,m)))
                  IF (node(im1,jm1,km1) .EQ. SOLID) THEN
                     phiIN= (fplus(bb(f_Comps(iComm,m)),im1,jm1,km1)/rho(im1,jm1,km1) - wt(bb(f_Comps(iComm,m)))*Delta)	&	! scalar contribution from inlet/outlet to current node
                            *phiTemp(im1,jm1,km1)		
                     phiOUT= (fplus(f_Comps(iComm,m),i,j,k)/rho(i,j,k) - wt(f_Comps(iComm,m))*Delta)*phiTemp(i,j,k)		! scalar contribution from current node to inlet/outlet
                     phiInOut= phiInOut + (phiOUT - phiIN)
                  END IF
               END DO
            END IF
         END DO
      END DO
   END IF
END DO
!===================================================================================================
END SUBROUTINE ScalarInOut
!===================================================================================================








!===================================================================================================
END MODULE PassiveScalar
!===================================================================================================
