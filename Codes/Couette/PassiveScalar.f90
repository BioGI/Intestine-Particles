!==================================================================================================
MODULE PassiveScalar		! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs)
!==================================================================================================
USE SetPrecision
USE Setup
USE IC
USE BClbm
USE BCscalar

IMPLICIT NONE 

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar_Setup	! sets up the passive scalar component
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

!-----initialize arrays
phi= 0.0_dbl								! scalar
phiTemp= 0.0_dbl							! temporary scalar

!-----scalar parameters
Dm= nuL/Sc								! binary molecular diffusivity (scalar in fluid)
Dmcf= (zcf*zcf)/tcf							! conversion factor for diffusivity
Delta= 1.0_dbl - 6.0_dbl*Dm						! scalar diffusion parameter

!-----set scalar values at sources/sinks
phiWall= 0.0_dbl							! value of scalar at the boundary

!-----set the scalar standard devation for gaussian distributions
sigma= 0.1_dbl*D							! 1/10 of the Diameter

!-----determine scalar starting iteration
phiStart= NINT((phiPer*Tmix)/tcf)
IF (phiPer.EQ.0.0) THEN
   phiStart= NINT((phiPer*Tmix)/tcf)+1 					! Balaji changed this to add 1 as for phiPer=0, phiSTart=0. But iter never has a value 0.
ENDIF

phiInNodes = 0.0_dbl
phiOutNodes = 0.0_dbl
!------------------------------------------------
END SUBROUTINE Scalar_Setup
!------------------------------------------------





!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar								! calculates the evolution of scalar in the domain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1				! index variables
REAL(dbl) :: phiBC						! scalar contribution from boundary
REAL(dbl) :: phiOutSurf,phiInSurf				! scalar contribution coming from and going into the boundary
REAL(dbl) :: tausgs						! contribution form tau_sgs term from particle closure
REAL(dbl) :: zcf3						! Cell volume

CALL ScalarDistribution						! sets/maintains initial distributions of scalar [MODULE: ICBC.f90]

!----- store the previous scalar values
phiTemp = phi
Negative_phi_Counter = 0
Negative_phi_Worst= 0

!----- Stream the scalar
DO k=1,nzSub
  DO j=1,nySub
    DO i=1,nxSub
      
      IF(node(i,j,k) .EQ. FLUID) THEN
        phi(i,j,k) = Delta*phiTemp(i,j,k)  
	phi(i,j,k) = phi(i,j,k)+ delphi_particle(i,j,k) 	! Balaji added to introduce drug concentration release
! 	Farhad Removed SGS effects
!	tausgs = ((tausgs_particle_x(i+1,j,k)-tausgs_particle_x(i-1,j,k)) + &
!		 (tausgs_particle_y(i,j+1,k)-tausgs_particle_y(i,j-1,k)) + &
!		 (tausgs_particle_z(i,j,k+1)-tausgs_particle_z(i,j,k-1)))*0.5_dbl
!
!	phi(i,j,k) = phi(i,j,k)+ tausgs 			! Balaji added - to handle SGS particle effects.
        DO m=0,NumDistDirs
          !-----i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)

         IF (km1.eq.0) then
             km1=nz ! Balaji added
          ELSE IF(km1.eq.nz+1) then
             km1=1 ! Balaji added
          END IF

          IF (jm1.eq.0) then
             jm1=ny ! Balaji added
          ELSE IF(jm1.eq.ny+1) then
             jm1=1 ! Balaji added
          END IF


          IF (node(im1,jm1,km1) .EQ. FLUID) THEN 
	     phi(i,j,k) = phi(i,j,k) + (fplus(m,im1,jm1,km1)/rho(im1,jm1,km1) - wt(m)*Delta)*phiTemp(im1,jm1,km1)
          ELSE IF((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN				! macro- boundary
!            CALL ScalarBC2(m,i,j,k,im1,jm1,km1,phiBC,phiOutSurf,phiInSurf)	! Wang: scalar flux  BC (using BB f's)	
             CALL Scalar_Fixed_BC(m,i,j,k,im1,jm1,km1,phiBC)				! Wang: scalar BC
             phi(i,j,k) = phi(i,j,k) + phiBC     
!            CALL AbsorbedScalarS2(i,j,k,m,phiOutSurf,phiInSurf)		! measure the absorption rate
             CALL AbsorbedScalarS(i,j,k,m,phiBC)				! measure the absorption rate
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN				! villi
             CALL ScalarBCV(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),phiBC)						! implement scalar boundary condition (using BB f's)	[MODULE: ICBC]
             phi(i,j,k) = phi(i,j,k) + phiBC     
             CALL AbsorbedScalarV(i,j,k,m,phiBC)									! measure the absorption rate
          ELSE
             OPEN(1000,FILE="error.txt")
             WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node(im1,jm1,km1) is out of range"
             WRITE(1000,*) "iter",iter
             WRITE(1000,*) "m=",m
             WRITE(1000,*) "i=",i,"j=",j,"k=",k
             WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
             WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
             WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
             WRITE(1000,*) "node(i,j,k)=",node(i,j,k)
             WRITE(1000,*) "node(im1,jm1,km1)=",node(im1,jm1,km1)
             CLOSE(1000)
             STOP
          END IF

        END DO

	!phi(i,j,k) = phi(i,j,k) + delphi_particle(i,j,k) ! Balaji added to introduce drug concentration release
       	!fix spurious oscillations in moment propagation method for high Sc #s

!------ node volume in physical units (cm^3) so when printing the drung units are "mole
        zcf3 = zcf*zcf*zcf

!------ Monitoring the negative phi
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

!----- Add the amount of scalar absorbed through the outer and villous surfaces
phiAbsorbed = 	phiAbsorbedS + phiAbsorbedV										! total amount of scalar absorbed up to current time
      
!------------------------------------------------
END SUBROUTINE Scalar
!------------------------------------------------




!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarS(i,j,k,m,phiBC)		! measures the total absorbed scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m						! index variables
REAL(dbl), INTENT(IN) :: phiBC   		  				! scalar contribution from the boundary condition
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

phiIN= phiBC									! contribution from the wall to the crrent node (in)
phiOUT= (fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)	! contribution to the wall from the current node (out)
phiAbsorbedS = phiAbsorbedS + (phiOUT - phiIN)					!- wt(m)*Delta*phiWall	! add the amount of scalar that has been absorbed at the current location in the current direction
!------------------------------------------------
END SUBROUTINE AbsorbedScalarS
!------------------------------------------------





!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarS2(i,j,k,m,phiOutSurf,phiInSurf)		! measures the total absorbed scalar ! New method by Balaji, Dec 2014
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m					! index variables
REAL(dbl), INTENT(IN) :: phiOutSurf,phiInSurf     			! scalar contribution from the boundary condition
REAL(dbl) :: phiOUT, phiIN						! scalar values exchanged with the wall

phiIN 	= phiInSurf							! contribution from the wall to the crrent node (in)
phiOUT	= phiOutSurf							! contribution to the wall from the current node (out)
phiAbsorbedS = phiAbsorbedS + (phiOUT - phiIN) 				!- wt(m)*Delta*phiWall	! add the amount of scalar that has been absorbed at the current location in the current direction
!------------------------------------------------
END SUBROUTINE AbsorbedScalarS2
!------------------------------------------------




!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarV(i,j,k,m,phiBC)		! measures the total absorbed scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m						! index variables
REAL(dbl), INTENT(IN) :: phiBC     						! scalar contribution from the boundary condition
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

phiIN 	= phiBC									! contribution from the wall to the crrent node (in)
phiOUT	= (fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)	! contribution to the wall from the current node (out)
phiAbsorbedV = phiAbsorbedV + (phiOUT - phiIN)					! add the amount of scalar that has been absorbed at the current location in the current direction
!------------------------------------------------
END SUBROUTINE AbsorbedScalarV
!------------------------------------------------





!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarInOut	   ! measures scalar left or entered the domain through the inlet or outlet
!--------------------------------------------------------------------------------------------------
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
               IF ((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN
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
!------------------------------------------------
END SUBROUTINE ScalarInOut
!------------------------------------------------


!================================================
END MODULE PassiveScalar
!================================================
