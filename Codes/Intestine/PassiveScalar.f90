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

INTEGER(lng):: i,j,k,m,im1,jm1,km1,mpierr							! index variables
INTEGER(lng):: Over_Sat_Counter, Over_Sat_Counter_Global
REAL(dbl)   :: Over_sat_Total,   Over_Sat_Total_Global
REAL(dbl)   :: Largest_phi, Largest_phi_Global							! OverSaturation issue monitoring
REAL(dbl)   :: phiBC 										! scalar contribution from boundary
REAL(dbl)   :: phiOutSurf,phiInSurf								! scalar contribution coming from and going into the boundary
REAL(dbl)   :: tausgs										! contribution form tau_sgs term from particle closure
REAL(dbl)   :: zcf3										! Cell volume

CALL IC_Drug_Distribution 									! sets/maintains initial distributions of scalar [MODULE: IC.f90]

IF (iter .EQ. phiStart) THEN  						 			! Calculate the intial amount of scalar
   phiTotal = 0.0_dbl
   DO k=1,nzSub
      DO j=1,nySub
         DO i=1,nxSub
            IF (node(i,j,k) .EQ. FLUID) THEN
               phiTotal = phiTotal + phi(i,j,k)
            END IF
         END DO
      END DO
   END DO
END IF


!----- store the previous scalar values
phiTemp 	    = phi + delphi_particle

Negative_phi_Counter= 0
Negative_phi_Worst  = 0.0_dbl
Negative_phi_Total  = 0.0_dbl

Over_Sat_Counter    = 0
Largest_phi	    = 0.0_dbl
Over_Sat_Total      = 0.0_dbl

!----- Stream the scalar
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         IF (node(i,j,k) .EQ. FLUID) THEN
            phi(i,j,k) = Delta * phiTemp(i,j,k)  

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

!---------- Monitoring the negative phi
            zcf3 = zcf*zcf*zcf* 1000000.0_dbl  											! node volume in physical units (cm^3) so drung units are "mole
            IF (phi(i,j,k) .LT. 0.0_dbl) THEN
               Negative_phi_Counter = Negative_phi_Counter + 1
               Negative_phi_Total   = Negative_phi_Total   + phi(i,j,k)  
               IF (phi(i,j,k) .LT. Negative_phi_Worst) THEN
                  Negative_phi_Worst = phi(i,j,k)
               END IF
               IF ( Flag_Rectify_Neg_phi) THEN
                  phi(i,j,k) = 0.0_dbl
               END IF                    
            END IF

!---------- Monitoring over saturation
            IF (phi(i,j,k) .GT. Cs_mol) THEN
               Over_Sat_Counter= Over_Sat_Counter + 1
               Over_Sat_Total  = Over_Sat_Total   + phi(i,j,k)
               IF (Largest_phi .LT. phi(i,j,k) ) THEN
                  Largest_phi= phi(i,j,k)
               END IF
            END IF

         END IF
      END DO
   END DO
END DO


!----- Printing out the outputs for Monitoring Negative-phi issue
CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
CALL MPI_ALLREDUCE(Negative_phi_Counter, Negative_phi_Counter_Global, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(Negative_phi_Total,   Negative_phi_Total_Global,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(Negative_phi_Worst,   Negative_phi_Worst_Global,   1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)

IF (myid .EQ. master) THEN
   IF (Negative_phi_Counter_Global .GE. 1) THEN
       write(2118,*) iter, Negative_phi_Counter_Global, Negative_phi_Total_Global/Cs_mol, Negative_phi_Worst_Global/Cs_mol, Negative_phi_Total_Global/(Cs_mol*Negative_phi_Counter_Global)
       CALL FLUSH(2118)
   END IF
END IF


!----- Monitoring the Over Saturation problem
CALL MPI_ALLREDUCE(Over_Sat_Counter, Over_Sat_Counter_Global, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(Over_Sat_Total,   Over_Sat_Total_Global,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(Largest_phi,      Largest_phi_Global,      1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)

IF (myid .EQ. master) THEN
   IF (Over_Sat_Counter_Global .GE. 1) THEN 
       write(2119,*) iter, Over_Sat_Counter_Global, Largest_phi_Global/Cs_mol, Over_Sat_Total_Global/(Over_Sat_Counter_Global*Cs_mol)
       CALL FLUSH(2119)
   END IF
END IF

!----- Add the amount of scalar absorbed through the outer surfaces
phiAbsorbed = phiAbsorbedS 						
!===================================================================================================
END SUBROUTINE Scalar
!===================================================================================================





!===================================================================================================
END MODULE PassiveScalar
!===================================================================================================
