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
Dm   = nuL/Sc							! binary molecular diffusivity (scalar in fluid)
Dmcf = (zcf*zcf)/tcf						! conversion factor for diffusivity
Delta= 1.0_dbl- 6.0_dbl*Dm					! scalar diffusion parameter

!---- scalar standard devation for gaussian distributions
sigma  = 0.1_dbl*D						! 1/10 of Diameter
!===================================================================================================
END SUBROUTINE Scalar_Setup
!===================================================================================================




!===================================================================================================
SUBROUTINE Scalar				! calculates the evolution of scalar in the domain
!===================================================================================================
IMPLICIT NONE

INTEGER(lng):: i,j,k,m,im1,jm1,km1,mpierr
INTEGER(lng):: iamBoundary(136,136,217)
REAL(dbl)   :: dV,dM         
INTEGER(lng):: Over_Sat_Counter, Over_Sat_Counter_Global
REAL(dbl)   :: Over_sat_Total,   Over_Sat_Total_Global
REAL(dbl)   :: Largest_phi, Largest_phi_Global							! OverSaturation issue monitoring
REAL(dbl)   :: phiBC,P_Astar_Bstar                          ! scalar contribution from boundary
REAL(dbl)   :: phiOutSurf,phiInSurf                         ! scalar contribution coming from and going into the boundary
REAL(dbl)   :: tausgs                                       ! contribution form tau_sgs term from particle closure
REAL(dbl)   :: zcf3                                         ! Cell volume
REAL(dbl)   :: phiIN,phiOUT
REAL(dbl)   :: P_A_O,P_A_B
REAL(dbl)   :: Sum_numerator,Sum_denominator,Sum_numerator_g,Sum_denominator_g
REAL(dbl)   :: alpha_BC
REAL(dbl)   :: q

CALL IC_Drug_Distribution 									! sets/maintains initial distributions of scalar [MODULE: IC.f90]

!----- store the previous scalar values
phiTemp             = phi + delphi_particle
Negative_phi_Counter= 0
Negative_phi_Worst  = 0.0_dbl
Negative_phi_Total  = 0.0_dbl
Over_Sat_Counter    = 0
Largest_phi	        = 0.0_dbl
Over_Sat_Total      = 0.0_dbl
!phiTotal            = 0.0_dbl
iamBoundary=0
dV= (100.0*zcf)**3.0_dbl

!----- Stream the scalar
IF ((coeffGrad .EQ. 1.0) .AND. (coeffPhi .EQ. 0.0) .AND. (coeffConst .EQ. 0.0) ) THEN ! No Flux or Permeability
   Sum_numerator  =0.0_dbl
   Sum_denominator=0.0_dbl
   alpha_BC=1.0_dbl
   DO k=1,nzSub
      DO j=1,nySub
         DO i=1,nxSub
            IF (node(i,j,k) .EQ. FLUID) THEN
               DO m=0,NumDistDirs       !  neighboring node 
                  im1= i- ex(m)
                  jm1= j- ey(m)
                  km1= k- ez(m)
                  IF (node(im1,jm1,km1) .EQ. SOLID) THEN		! Communicating with a solid node acroos the boundary
                     P_A_O          =(fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta) * phiTemp(i,j,k)
                     P_A_B          =(fplus(m,i,j,k)    /rho(i,j,k) - wt(m)    *Delta) * phiTemp(i,j,k)
                     CALL BC_Zero_Flux(m,i,j,k,im1,jm1,km1,q,phiBC,P_Astar_Bstar,alpha_BC) 
                     Sum_numerator  =Sum_numerator +  P_A_O + (1-q/q)*P_A_B
                     Sum_denominator=Sum_denominator+ P_Astar_Bstar/q
                  END IF
               END DO
            END IF   
         END DO 
      END DO
   END DO
   CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   CALL MPI_ALLREDUCE(Sum_numerator,  Sum_numerator_g,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
   CALL MPI_ALLREDUCE(Sum_denominator,Sum_denominator_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
   alpha_BC= Sum_numerator_g/Sum_denominator_g
   IF (myid.EQ.master) THEN
       WRITE(*,*) 'iter,up,down,alpha',iter,Sum_numerator_g,Sum_denominator_g,alpha_BC
   ENDIF
END IF 


!----- Stream the scalar
N0=0
N1=0
N2=0
N3=0
N4=0
N5=0
N6=0
N7=0
N8=0
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         IF (node(i,j,k) .EQ. FLUID) THEN
            IF (iter .EQ. iter_Start_phi) THEN
               phiTotal = phiTotal + phi(i,j,k)
            ENDIF
            phi(i,j,k) = Delta * phiTemp(i,j,k) 

            DO m=0,NumDistDirs
               !-----  neighboring node --------------------------------------------------------------
               im1= i- ex(m)
               jm1= j- ey(m)
               km1= k- ez(m)
               IF (node(im1,jm1,km1) .EQ. FLUID) THEN 
                  phi(i,j,k) = phi(i,j,k) + (fplus(m,im1,jm1,km1)/rho(im1,jm1,km1) - wt(m)*Delta)*phiTemp(im1,jm1,km1)
               ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN							                                                    		! iCommunicating with a solid node acroos the boundary
                  IF ((coeffGrad .EQ. 1.0) .AND. (coeffPhi .EQ. 0.0) .AND. (coeffConst .EQ. 0.0) ) THEN ! No Flux or Permeability
                     iamBoundary(i,j,k) = 1

                     CALL BC_Zero_Flux(m,i,j,k,im1,jm1,km1,q,phiBC,P_Astar_Bstar,alpha_BC) 
                     phi(i,j,k) = phi(i,j,k) + phiBC    
                     !phiIN = phiBC
                     !phiOUT= (fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)
                     !phiAbsorbedS = phiAbsorbedS + (phiOUT-phiIN)		                                                       		! scalar absorbed at current location in mth direction
                  ELSE                                                                                                        ! Immidiate uptake 
                     CALL BC_Scalar(m,i,j,k,im1,jm1,km1,phiBC) 
                     phi(i,j,k) = phi(i,j,k) + phiBC     
                     CALL AbsorbedScalarS(i,j,k,m,im1,jm1,km1,phiBC) 								! measure the absorption rate
                  ENDIF   
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
            IF (Flag_Rectify_Neg_phi) THEN
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
            ENDIF
         END IF
      END DO
   END DO
END DO

IF((iter-iter0) .LE. 2)THEN
   WRITE(*,*) 'CPU:N0,N1,N2,N3,N4,N5,N6,N7,N8',myid,N0,N1,N2,N3,N4,N5,N6,N7,N8
ENDIF

IF (Flag_2step_Permeability) THEN    
   DO k=1,nzSub
      DO j=1,nySub
         DO i=1,nxSub
            IF ((node(i,j,k) .EQ. FLUID).AND.(iamBoundary(i,j,k).EQ.1)) THEN
               dM=  phiTemp(i,j,k) * tcf * Pw * dA_permeability(k)
               phi(i,j,k) = phi(i,j,k) - dM /dV
               phiAbsorbedS =phiAbsorbedS + dM/dV
            ENDIF
         END DO
      END DO
   END DO
ENDIF

!----- Printing out the outputs for Monitoring Negative-phi issue
IF (Flag_Rectify_Neg_phi) THEN
   CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   CALL MPI_ALLREDUCE(Negative_phi_Counter, Negative_phi_Counter_Global, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, mpierr)
   CALL MPI_ALLREDUCE(Negative_phi_Total,   Negative_phi_Total_Global,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
   CALL MPI_ALLREDUCE(Negative_phi_Worst,   Negative_phi_Worst_Global,   1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)
   IF (myid .EQ. master) THEN
      IF (Negative_phi_Counter_Global .GE. 1) THEN
          write(2118,*) iter, Negative_phi_Counter_Global, Negative_phi_Total_Global/S_intrinsic, Negative_phi_Worst_Global/S_intrinsic, Negative_phi_Total_Global/(S_intrinsic*Negative_phi_Counter_Global)
          CALL FLUSH(2118)
      END IF
   END IF
ENDIF

!----- Monitoring the Over Saturation problem
!CALL MPI_ALLREDUCE(Over_Sat_Counter, Over_Sat_Counter_Global, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, mpierr)
!CALL MPI_ALLREDUCE(Over_Sat_Total,   Over_Sat_Total_Global,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
!CALL MPI_ALLREDUCE(Largest_phi,      Largest_phi_Global,      1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)

!IF (myid .EQ. master) THEN
!   IF (Over_Sat_Counter_Global .GE. 1) THEN 
!       write(2119,*) iter, Over_Sat_Counter_Global, Largest_phi_Global/Cs_mol, Over_Sat_Total_Global/(Over_Sat_Counter_Global*Cs_mol)
!       CALL FLUSH(2119)
!   END IF
!END IF

!----- Add the amount of scalar absorbed through the outer surfaces
phiAbsorbed = phiAbsorbedS 			
!===================================================================================================
END SUBROUTINE Scalar
!===================================================================================================





!===================================================================================================
END MODULE PassiveScalar
!===================================================================================================
