!===================================================================================================
MODULE Output						! Contains subroutines for printing to file
!===================================================================================================
USE SetPrecision
USE Setup
USE LBM
USE PassiveScalar
USE MPI			

IMPLICIT NONE 

CONTAINS



!===================================================================================================
SUBROUTINE Output_Setup								! sets up the output
!===================================================================================================
IMPLICIT NONE

ALLOCATE(filenum(0:nt))			! maximum number of output files (should only output ~numOuts times)
filenum = 0_lng				! initialize to 0
fileCount = 0_lng			! initialize to 0

ALLOCATE(parfilenum(0:nt))		! maximum number of output files (should only output ~numOuts times)
ALLOCATE(numparticleSubfile(0:nt))	! maximum number of output files (should only output ~numOuts times)
parfilenum = 0_lng			! initialize to 0
parfileCount = 0_lng			! initialize to 0
numparticleSubfile = 0_lng

ALLOCATE(radius(0:nz+1,0:500))		! 500 is an arbitrarily large number of output iterations...
radcount = 0_lng			! initialize the output count
!===================================================================================================
END SUBROUTINE Output_Setup
!===================================================================================================








!===================================================================================================
SUBROUTINE OpenOutputFiles							! opens output files
!===================================================================================================
IMPLICIT NONE

IF (myid .EQ. master) THEN
   !----- Monitoring negative phi----
   OPEN(2118,FILE='Negative-phi.dat',POSITION='APPEND')
   WRITE(2118,'(A120)') 'VARIABLES = iter,  Number of Negative phi Nodes,  Total Sum of Negative phi,  Worst Negative phi,  Average of Negative phi'
   CALL FLUSH(2118)

   !----- Monitoring over saturation
   OPEN(2119,FILE='Over_Saturation.dat',POSITION='APPEND')
   WRITE(2119,'(A120)') 'VARIABLES = iter,  Number of OverSaturated Nodes,  Worst Oversaturation'
   CALL FLUSH(2119)

   !----- Monitoring computational costs
   OPEN(5,FILE='status.dat')										
   CALL FLUSH(5)													

  !----- Surface Area-----
  !OPEN(2474,FILE='SA.dat',POSITION='APPEND')
  !WRITE(2474,'(A36)') 'VARIABLES = "period", "SA"'
  !WRITE(2474,*) 'ZONE F=POINT'
  !CALL FLUSH(2474)

  !----- Walll Flux-----
  !OPEN(4748,FILE='wall_flux.dat')
  !WRITE(4748,*) 'VARIABLES = "Axial Distance", "Flux"'
  !CALL FLUSH(4748)

  !---- Volume -----
  OPEN(2460,FILE='volume.dat')
  WRITE(2460,*) 'VARIABLES = "period", "volume"'
  WRITE(2460,*) 'ZONE F=POINT'
  CALL FLUSH(2460)

  !----- Drug Conservation
  OPEN(2472,FILE='Drug-Conservation-'//sub//'.dat')
  WRITE(2472,'(A120)') '#VARIABLES =iter,time, Drug_Initial, Drug_Released_Total, Drug_Absorbed, Drug_Remained_in_Domain, Drug_Loss_Percent, Drug_Loss_Modified_Percent'
  WRITE(2472,*) '#ZONE F=POINT'
  CALL FLUSH(2472)

  !----- Diensity Correction to improve mass conservation
  OPEN(2473,FILE='Density_Correction.dat')
  WRITE(2473,*) 'VARIABLES: iter, Density Correction'
  CALL FLUSH(2473)
END IF

!----- Mass
OPEN(2458,FILE='mass-'//sub//'.dat')
WRITE(2458,'(A120)') '#VARIABLES = period, time, mass_theory, mass_actual, mass_err'
WRITE(2458,*) '#ZONE F=POINT'
CALL FLUSH(2458)

!---- Test Output
!OPEN(9,FILE='testoutput-'//sub//'.dat',POSITION='append')
!CALL FLUSH(9)

!===================================================================================================
END SUBROUTINE OpenOutputFiles
!===================================================================================================








!===================================================================================================
SUBROUTINE CloseOutputFiles	! opens output files
!===================================================================================================
IMPLICIT NONE

IF(myid .EQ. master) THEN
  CLOSE(2118) 			! Monitoring negative phi
  CLOSE(2119)                   ! Monitorin Over Saturation
  CLOSE(5)			! Status								
! CLOSE(2474)			! Surface Area
  CLOSE(2460)			! Volume
END IF
CLOSE(2458)			! Mass
CLOSE(2472)			! Scalar
CLOSE(2473)                     ! Density Corrections
!CLOSE(9) 			! Test Output
!===================================================================================================
END SUBROUTINE CloseOutputFiles
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintStatus					! Writes Computational Costs to a File
!===================================================================================================
IMPLICIT NONE

IF (myid .EQ. master) THEN
   IF (iter .EQ. iter0-1_lng) THEN
      rate = 100_lng							! Set the rate of counting
      CALL SYSTEM_CLOCK(start,rate)					! Keep Track of Elapsed Time 
      WRITE(5,*) 'Running Simulation...'				! Print Current Status
      WRITE(5,*) '[',nt,'Iterations]'
      CALL FLUSH(5)													
   END IF
 
   IF (MOD(iter,2) .EQ. 0 .AND. (iter .NE. 0)) THEN 			! Print Current Status Periodically
      CALL SYSTEM_CLOCK(current,rate)					
      WRITE(5,*) iter, ((current-start)/REAL(rate))/(iter-iter0)
      CALL FLUSH(5)													
   END IF

   IF (iter .EQ. nt) THEN
      CALL SYSTEM_CLOCK(final,rate)					! End the Timer
      WRITE(5,*) 'Simulation Complete.'
      WRITE(5,*) 'Total Wall-Time (min):', ((final-start)/REAL(rate))/60.0_dbl
      CALL FLUSH(5)													
   END IF
END IF
!===================================================================================================
END SUBROUTINE PrintStatus
!===================================================================================================









!===================================================================================================
SUBROUTINE PrintRestart							      ! prints restart files
!===================================================================================================
IMPLICIT NONE

TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next
INTEGER(lng) :: i,j,k,m			
CHARACTER(7):: iter_char                        ! iteration stored as a character

!----- Creating a file called iter0.dat with iteration number in it --------------------------------
IF (myid .EQ. master) THEN
   OPEN(550,FILE='Restart-iter.dat')
   WRITE(550,*) iter
   CLOSE(550)
END IF

!----- Creating a file called restart-iter-sub.dat with all macroscopic fields and distribution functions -----
WRITE(iter_char(1:7),'(I7.7)') iter-1_lng
OPEN(500,FILE='Restart-Out-'//iter_char//'-'//sub//'.dat')			
DO k=0,nzSub+1
   DO j=0,nySub+1
      DO i=0,nxSub+1
         WRITE(500,'(I1)') node(i,j,k)
         WRITE(500,'(F8.6)') u(i,j,k)
         WRITE(500,'(F8.6)') v(i,j,k)
         WRITE(500,'(F8.6)') w(i,j,k)
         WRITE(500,'(F8.6)') rho(i,j,k)
         WRITE(500,'(F8.6)') phi(i,j,k)
         DO m=0,NumDistDirs
            WRITE(500,'(F10.8)') f(m,i,j,k)
         END DO
      END DO
   END DO
END DO
WRITE(500,*) Drug_Initial
WRITE(500,*) Drug_Released_Total
WRITE(500,*) Drug_Absorbed
WRITE(500,*) Drug_Remained_in_Domain
CLOSE(500)


!----- Creating a file called particle-restart-iter.dat with all the particle data in it ----------------
IF (myid .eq. master) THEN
   OPEN(156,FILE='Restart-Particles-'//iter_char//'.dat')
   write(156,*) np
   current => ParListHead%next							 
   DO WHILE (ASSOCIATED(current))
      next => current%next 	
      WRITE(156,*)   	current%pardata%parid,	        &
		     	current%pardata%xp, 	  	&
			current%pardata%yp,  	  	&
			current%pardata%zp, 	  	&
                        current%pardata%up, 	  	&
		       	current%pardata%vp, 	  	&
			current%pardata%wp, 	  	&
			current%pardata%rp, 	  	&
                        current%pardata%xpold,          &
                        current%pardata%ypold,          &
                        current%pardata%zpold,          &
                        current%pardata%upold,          &
                        current%pardata%vpold,          &
                        current%pardata%wpold,          &
                        current%pardata%rpold,          &
                        current%pardata%par_conc,       &
                        current%pardata%gamma_cont,	&
                        current%pardata%sh,	        &
                        current%pardata%S,    		&
                        current%pardata%Sst,   		&
                        current%pardata%Veff,           &
			current%pardata%Nbj,		&		
			current%pardata%bulk_conc,	&
			current%pardata%delNBbyCV, 	&
		 	current%pardata%cur_part,	&
		 	current%pardata%new_part  		
1001 format (I5,24F18.10,2I5)
     current => next  
  ENDDO
  CLOSE(160)
END IF


!===================================================================================================
END SUBROUTINE PrintRestart
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintFields			       ! print velocity, density, scalar to output files
!===================================================================================================
IMPLICIT NONE

INTEGER(lng):: i,j,k,ii,jj,kk,n			! index variables (local and global)
CHARACTER(7):: iter_char			! iteration stored as a character
REAL(lng)   :: pressure				

IF ((MOD(iter,(((nt+1_lng)-iter0)/numOuts)) .EQ. 0) .OR. &
   (iter .EQ. iter0-1_lng) .OR. (iter .EQ. iter0)  .OR. &
   (iter .EQ. phiStart) .OR. (iter .EQ. nt)) THEN
   !----- scale the iteration by 1/10 such that the numbers used in the output file aren't too large
   WRITE(iter_char(1:7),'(I7.7)') iter

   !----- store the current iteration in "filenum"
   filenum(fileCount) = iter
   fileCount = fileCount + 1_lng

   !----- open the proper output file
   OPEN(60,FILE='out-'//iter_char//'-'//sub//'.dat')
   WRITE(60,*) 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
   WRITE(60,'(A10,E15.5,A5,I4,A5,I4,A5,I4,A8)') 'ZONE T="',iter/(nt/nPers),'" I=',nxSub,' J=',nySub,' K=',nzSub,'F=POINT'
   DO k=1,nzSub
      DO j=1,nySub
         DO i=1,nxSub
            !----- convert local i,j,k, to global ii,jj,kk
            ii = ((iMin - 1_lng) + i)
            jj = ((jMin - 1_lng) + j)
            kk = ((kMin - 1_lng) + k)
            IF (phi(i,j,k) .LT. 1.0e-18) THEN
               phi(i,j,k)=0.0_lng
            END IF   
            pressure= (rho(i,j,k)-denL)*dcf*pcf
            WRITE(60,'(I3,2I4,3F11.7,E13.4,E12.4,I2)') ii,jj,kk, u(i,j,k)*vcf, v(i,j,k)*vcf,         &
                                                       w(i,j,k)*vcf, pressure, phi(i,j,k), node(i,j,k)
         END DO
      END DO
   END DO

   CLOSE(60)
 
   IF (myid .EQ. master) THEN 			 ! Store radius at this iteration     
      DO k=0,nz+1
         radius(k,radcount) = rDom(k)
      END DO
      radcount = radcount + 1_lng 
   END IF
END IF
!===================================================================================================
END SUBROUTINE PrintFields
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintParticles	! print particle position, velocity, radius, and concentrationto output files
!===================================================================================================
IMPLICIT NONE
INTEGER(lng) 	:: numParticlesSub
INTEGER(lng)	:: i,j,k,ii,jj,kk,n,xaxis,yaxis		! index variables (local and global)
CHARACTER(7)	:: iter_char				! iteration stored as a character
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

xaxis=ANINT(0.5_dbl*(nx+1))
yaxis=ANINT(0.5_dbl*(ny+1))

IF ((MOD(iter,(((nt+1_lng)-iter0)/numOuts)) .EQ. 0) &
   .OR. (iter .EQ. iter0-1_lng) .OR. (iter .EQ. iter0) &
   .OR. (iter .EQ. phiStart) .OR. (iter .EQ. nt)) THEN
   !------ scale the iteration by 1/10 such that the numbers used in the output file aren't too large
   WRITE(iter_char(1:7),'(I7.7)') iter

   !------ store the current iteration in "parfilenum"
   parfilenum(parfileCount) = iter
   numParticlesSub  = 0_lng

   !------ open the proper output file
   OPEN(160,FILE='pardat-'//iter_char//'-'//sub//'.csv')
   WRITE(160,*) '"CPU","x","y","z","u","v","w","ParID","Sh","rp","Cb/Cs","delNBbyCV","Sst","S","Veff","Nbj"'
   current => ParListHead%next							 

   DO WHILE (ASSOCIATED(current))
      numParticlesSub = numParticlesSub + 1_lng
      next => current%next 							! copy pointer of next node
      WRITE(160,1001)   current%pardata%cur_part  	,',', 	&
		 	current%pardata%xp 	  	,',',	&
			current%pardata%yp  	  	,',',	&
			current%pardata%zp 	  	,',',	&
                        current%pardata%up*vcf 	  	,',',	&
		       	current%pardata%vp*vcf 	  	,',',	&
			current%pardata%wp*vcf 	  	,',',	&
                        current%pardata%parid 	  	,',',	&
			current%pardata%sh 	  	,',',	&
			current%pardata%rp/xcf 	  	,',',	&
			current%pardata%bulk_conc/Cs_mol,',', 	&
			current%pardata%delNBbyCV 	,',', 	&
			current%pardata%Sst 	  	,',',	&
			current%pardata%S 	  	,',',	&
			current%pardata%Veff 	  	,',',	&
			current%pardata%Nbj
1001 format (I4,a2,F9.4,a2,F9.4,a2,F9.4,a2,F10.6,a2,F10.6,a2,F10.6,a2,I5,a2,F12.8,a2,F15.10,a2,F15.7,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2)
     current => next   						! point to next node in the list
  ENDDO

  CLOSE(160)
  numparticleSubfile(parfileCount) = numParticlesSub	
  parfileCount = parfileCount + 1_lng
ENDIF


!===================================================================================================
END SUBROUTINE PrintParticles	
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintTime	! print time information for scalability information
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: dest, src, tag			! send/recv variables: destination, source, message tag
INTEGER(lng) :: stat(MPI_STATUS_SIZE)		! status object: rank and tag of the sending processing unit/message
INTEGER(lng) :: mpierr				! MPI standard error variable

IF(iter .EQ. 0) THEN				! start the "timer"
  tStart = MPI_WTIME()				! start time [intrinsic: MPI_WTIME() tracks wall time]
ELSE IF(iter .EQ. nt+1) THEN
  tEnd	= MPI_WTIME()				! stop the "timer" 
  tTotal	= (tEnd - tStart)		! calculate the time (for each processing unit)

  !----- Send times from each processing unit to the master processing unit
  dest	= master				! send to master
  tag	= 50 					! message tag (arbitrary)

  IF(myid .NE. master) THEN			! send time to master
    CALL MPI_SEND(tTotal,1,MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,mpierr)
  ELSE
    !----- Set up time information output file
    OPEN(871,FILE='time.dat')  
    WRITE(871,*) 'TIMING INFORMATION'
    WRITE(871,*) '---------------------------------'
    WRITE(871,'(1X,"Subdomain Layout: ",I2," x ",I2," x ",I2)') NumSubsX, NumSubsY, NumSubsZ
    WRITE(871,'(1X,"Number of Time Steps: ",I8)') nt 
    WRITE(871,*)
    WRITE(871,'(1X,"Processing unit ",I2," took ",1F15.5," seconds.")') myid, tTotal

    tSum = tTotal										! initialize the sum to the time for the master processing unit
    DO src = 1,(numprocs-1)
       CALL MPI_RECV(tRecv,1,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,stat,mpierr)		! recieve time from each processing unit
       WRITE(871,'(1X,"Processing unit ",I2," took ",1F15.5," seconds.")') src, tRecv
       tSum = tSum + tRecv									! add the time from the sourceth processing unit to the sum
    END DO

    WRITE(871,*)
    WRITE(871,'(1X,"The sum of the wall times for each proccessing unit is: ",1F15.5," seconds.")') tSum
    WRITE(871,*)
    WRITE(871,'(1X,"The average wall time per processing unit is: ",5X,1F15.5," seconds.")') tSum/numprocs
    CLOSE(871)
  END IF
  CALL PrintFields							! output the velocity, density, and scalar fields
ELSE
  OPEN(1000,FILE="error.dat")
  WRITE(1000,*) 'PrintTime has been called, but iter is .NE. to iter0 or nt+1'
  WRITE(1000,*) 'iter=', iter
  CLOSE(1000)
  STOP
END IF

!===================================================================================================
END SUBROUTINE PrintTime
!===================================================================================================







!===================================================================================================
SUBROUTINE PrintMass					! checks the total mass in the system 
!===================================================================================================
IMPLICIT NONE

INTEGER(lng):: i,j,k				! index variables
REAL(dbl)   :: node_volume 			! Volume of each lattice cell
REAL(dbl)   :: volume		 		! domain volume 
REAL(dbl)   :: mass_theory			! mass in the system based on uniform density 
REAL(dbl)   :: mass_actual			! mass in the system based on local density 
REAL(dbl)   :: mass_err				! mass error due to denisty variations

node_volume= xcf*ycf*zcf			! calculate the node volume

mass_actual= 0.0_dbl
volume    = 0.0_dbl

!----- calculate the mass in the system based on the density and the number of fluid nodes
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         IF (node(i,j,k) .EQ. FLUID) THEN
            volume	= volume      + node_volume
            mass_actual = mass_actual + node_volume *rho(i,j,k)*dcf
         END IF 
      END DO
   END DO
END DO

mass_theory = volume *den 

mass_err= 100*(mass_theory-mass_actual)/mass_theory

WRITE(2458,'(I8,6E21.12)') iter, iter*tcf, mass_theory, mass_actual, mass_err 
CALL FLUSH(2458)  

!===================================================================================================
END SUBROUTINE PrintMass
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintVolume						!volume as a function of time
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k						! index variables
REAL(dbl) :: volume						! analytical volume

IF (myid .EQ. master) THEN
   volume = 0.0_dbl 						! initialize the volume
   DO k=1,nz							! cacluate volume in the system
      volume = volume + PI*rDom(k)*rDom(k)*zcf
   END DO
   WRITE(2460,'(I8,E15.7)') iter, volume
   CALL FLUSH(2460)  
END IF
!===================================================================================================
END SUBROUTINE PrintVolume
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintDrugConservation		! prints the total amount of scalar absorbed through the walls 
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k,mpierr				! index variables
INTEGER(lng) :: numFluids				! number of fluid nodes in the domain
REAL(dbl)    :: phiDomain, phiIC, Drug_Initial		! current amount of scalar in the domain
REAL(dbl)    :: phiAverage				! average scalar in the domain
REAL(dbl)    :: zcf3					! node volume in physical units
REAL(dbl)    :: phiTotal_Global, phiAbsorbed_Global, phiDomain_Global

TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

!----- Calculate the amount of scalar in the domain
numFluids = 0_lng
phiDomain = 0.0_dbl
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         IF (node(i,j,k) .EQ. FLUID) THEN
            phiDomain = phiDomain + phi(i,j,k)
            numFluids = numFluids + 1_lng
         END IF
      END DO
   END DO
END DO

!------ average scalar in the domain
IF (numFluids .GT. 1e-8) THEN
   phiAverage = phiDomain/numFluids		
ELSE
   phiAverage = 0.0_dbl
END IF

!------ node volume in physical units (cm3) which is consistent with phi units (mole/cm3) 
zcf3 =  1000000.0_dbl * zcf*zcf*zcf 

!------ Computing the total drug released from particles      
IF (myid .EQ. master) THEN
   IF ((ParticleTrack .EQ. ParticleOn) .AND. (iter .GE. phiStart)) THEN
      current => ParListHead%next
      DO WHILE (ASSOCIATED(current))
         next => current%next
         Drug_Released_Total = Drug_Released_Total + current%pardata%delNBbyCV * zcf3
         current => next
      ENDDO
   END IF
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
CALL MPI_ALLREDUCE(phiTotal   , phiTotal_Global   ,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(phiAbsorbed, phiAbsorbed_Global,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(phiDomain  , phiDomain_Global  ,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

Drug_Initial		= phiTotal_Global    * zcf3
Drug_Absorbed 		= phiAbsorbed_Global * zcf3 + Drug_Absorbed_restart
Drug_Remained_in_Domain = phiDomain_Global   * zcf3
Drug_Loss 		= (Drug_Released_Total+ Drug_Initial                            ) - (Drug_Absorbed + Drug_Remained_in_Domain)  
Drug_Loss_Modified 	= (Drug_Released_Total+ Drug_Initial- Negative_phi_Total_Global) - (Drug_Absorbed + Drug_Remained_in_Domain)

IF (Drug_Released_Total .LT. 1e-20) THEN
   Drug_Released_Total =1e-20
END IF

Drug_Loss_Percent = (Drug_Loss / (Drug_Released_Total+Drug_Initial)) * 100.0_lng
Drug_Loss_Modified_Percent = (Drug_Loss_Modified / (Drug_Released_Total+Drug_Initial)) * 100.0_lng  

IF (abs(Drug_Absorbed) .lt. 1.0e-40) THEN
   Drug_Absorbed = 0.0_lng
ENDIF

IF (myid .EQ. master) THEN
   WRITE(2472,'(I7, F9.3, 6E21.13)') iter, iter*tcf, Drug_Initial, Drug_Released_Total, Drug_Absorbed, Drug_Remained_in_Domain, Drug_Loss_Percent, Drug_Loss_Modified_Percent 
   CALL FLUSH(2472)
END IF
!===================================================================================================
END SUBROUTINE PrintDrugConservation 
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintParams	              ! prints the total amount of scalar absorbed through the walls 
!===================================================================================================
IMPLICIT NONE

IF (myid .EQ. 0) THEN
   OPEN(11,FILE='parameters.dat') 		! write input data to file
   WRITE(11,*) 'nx=',nx	 			! number of nodes in the x-direction
   WRITE(11,*) 'ny=',ny				! number of nodes in the y-direction
   WRITE(11,*) 'nz=',nz				! number of nodes in the z-direction
   WRITE(11,*)
   WRITE(11,*) 'NumSubsX=',NumSubsX		! number of subdomains in the X direction
   WRITE(11,*) 'NumSubsY=',NumSubsY		! number of subdomains in the Y direction
   WRITE(11,*) 'NumSubsY=',NumSubsZ		! number of subdomains in the Z direction
   WRITE(11,*)
   WRITE(11,*) 'NumSubsTotal=',NumSubsTotal	! number of total subdomains
   WRITE(11,*)
   WRITE(11,*) 'L=',L				! length
   WRITE(11,*) 'D=',D				! diameter
   WRITE(11,*)
   WRITE(11,*) 'epsOVERa1=',epsOVERa1		! peristaltic occlusion ratio (distance of occlusion/mean half-width)
   WRITE(11,*) 's1=',s1				! peristaltic wave speed
   WRITE(11,*) 'numw1=',numw1			! number of peristaltic waves
   WRITE(11,*) 'wc1=',wc1			! peristaltic weighting coefficient
   WRITE(11,*)
   WRITE(11,*) 'epsOVERa2=',epsOVERa2		! segmental occlusion ratio (distance of occlusion/mean half-width)
   WRITE(11,*) 'Ts=',Ts				! segmental contraction period
   WRITE(11,*) 'numw2=',numw2			! number of segmental waves
   WRITE(11,*) 'wc2=',wc2			! segmental weighting coefficient
   WRITE(11,*)
   WRITE(11,*) 'Tmix=',Tmix			! period of mixed mode simulation
   WRITE(11,*)
   WRITE(11,*) 'numVilliZ=',numVilliZ		! number of rows of villi in the Z-direction
   WRITE(11,*) 'numVilliTheta=',numVilliTheta	! number of villi in each row (theta-direction)
   WRITE(11,*) 
   WRITE(11,*) 'Lv=',Lv				! length of the villi (micrometers)
   WRITE(11,*) 'Rv=',Rv				! radius of the villi (micrometers)
   WRITE(11,*) 
   WRITE(11,*) 'freqRatioT=',freqRatioT		! villous frequency to macroscopic contraction frequency (azimuthal, theta direction)
   WRITE(11,*) 'freqRatioZ=',freqRatioZ		! villous frequency to macroscopic contraction frequency (axial, z direction)
   WRITE(11,*) 
   WRITE(11,*) 'randORord=',randORord		! flag to determine if the villous motion is random or ordered
   WRITE(11,*) 'villiAngle=',villiAngle		! maximum angle of active villous travel (degrees) 
   WRITE(11,*) 
   WRITE(11,*) 'den=',den			! density
   WRITE(11,*) 'nu=',nu				! kinematic viscosity
   WRITE(11,*) 'tau=',tau			! relaxation parameter
   WRITE(11,*)
   WRITE(11,*) 'Dm=',Dm*Dmcf			! diffusivity
   WRITE(11,*) 'Sc=',Sc				! Schmidt number
   WRITE(11,*) 'sclrIC=',sclrIC			! initial/maintained scalar distribution (1=BLOB,2=LINE,3=INLET,4=UNIFORM)
   WRITE(11,*) 'phiIC=',phiIC			! maximum scalar concentration
   WRITE(11,*) 'phiPer=',phiPer			! period at which to start the scalar
   WRITE(11,*) 'phiStart=',phiStart		! iteration at which to start the scalar
   WRITE(11,*)
   WRITE(11,*) 'nPers=',nPers			! total number of periods to run
   WRITE(11,*) 'numOuts=',numOuts		! number of output files (roughly)
   WRITE(11,*) 'restart=',restart		! use restart file? (0 if no, 1 if yes)
   WRITE(11,*)
   WRITE(11,*) 'xcf=', xcf			! x distance conversion factor
   WRITE(11,*) 'ycf=', ycf			! y distance conversion factor
   WRITE(11,*) 'zcf=', zcf			! z distance conversion factor
   WRITE(11,*)
   WRITE(11,*) 'tcf=', tcf			! time conversion factor
   WRITE(11,*) 'dcf=', dcf			! density conversion factor
   WRITE(11,*) 'vcf=', vcf			! velocity conversion factor
   WRITE(11,*) 'pcf=', pcf			! pressure conversion factor
   WRITE(11,*) 'Dmcf=',Dmcf			! diffusivity conversion factor
   WRITE(11,*)
   WRITE(11,*) 'nt=', nt			! number of time steps
   CLOSE(11)
END IF
!===================================================================================================
END SUBROUTINE PrintParams
!===================================================================================================








!===================================================================================================
SUBROUTINE MergeOutput	! combines the subdomain output files into one output file for the entire computational domain 
!===================================================================================================
IMPLICIT NONE

CALL MergeScalar
CALL MergeFields
CALL MergeMass
!===================================================================================================
END SUBROUTINE MergeOutput
!===================================================================================================








!===================================================================================================
SUBROUTINE MergeFields		! combines the subdomain output into an output file for the entire computational domain 
!===================================================================================================
IMPLICIT NONE

REAL(dbl), ALLOCATABLE	:: FieldData(:,:,:,:)			! u,v,w,density, and scalar for each node in the computational domain
!REAL(dbl), ALLOCATABLE	:: fluxField(:,:,:)			! scalar flux
!REAL(dbl), ALLOCATABLE	:: psi(:,:)				! stream function
INTEGER(lng) :: SubLimits(NumSubsTotal,3)			! array containing nxSub, nySub and nzSub for each subdomain
INTEGER(lng) :: nxSend(3),nxRecv(3)				! nx-,ny-,nzSub from each subdomain
INTEGER(lng) :: ii,jj,kk					! neighboring indices for writing
INTEGER(lng) :: i,j,k,n,nn,nnn					! loop variables
INTEGER(lng) :: dest, src, tag					! send/recv variables: destination, source, message tag
INTEGER(lng) :: stat(MPI_STATUS_SIZE)				! status object: rank and tag of the sending processing unit/message
INTEGER(lng) :: mpierr						! MPI standard error variable
INTEGER(lng) :: numLines					! number of lines to read
INTEGER(lng) :: combine1,combine2				! clock variables
CHARACTER(7) :: iter_char					! iteration stored as a character 
CHARACTER(5) :: nthSub						! current subdomain stored as a character

!---- send subdomain information to the master
tag = 60_lng										! starting message tag (arbitrary)

IF (myid .NE. master) THEN
   dest	= master									! send to master
   !----- fill out nxSend array
   nxSend(1) = nxSub
   nxSend(2) = nySub
   nxSend(3) = nzSub
   CALL MPI_SEND(nxSend(1:3),3,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,mpierr)		! send nx-,ny-,nzSub
ELSE

   ALLOCATE(FieldData(nx,ny,nz,6))

   !----- initialize FieldData to 0
   FieldData = 0.0_dbl

   !----- print combining status...
   CALL SYSTEM_CLOCK(combine1,rate)							! Restart the Timer
   OPEN(5,FILE='status.dat',POSITION='APPEND')										
   WRITE(5,*)
   WRITE(5,*)
   WRITE(5,*) 'Combining field output files and deleting originials...'
   WRITE(5,*)     
   CALL FLUSH(5)

   !----- fill out SubLimits for the master processor (1st subdomain)
   SubLimits(1,1) = nxSub
   SubLimits(1,2) = nySub
   SubLimits(1,3) = nzSub

   DO src = 1,(numprocs-1)
     CALL MPI_RECV(nxRecv(1:3),3,MPI_INTEGER,src,tag,MPI_COMM_WORLD,stat,mpierr)	! receive nx-,ny-,nzSub from each subdomain
     !----- fill out SubLimits for the master processor (1st subdomain)
     SubLimits(src+1,1) = nxRecv(1)
     SubLimits(src+1,2) = nxRecv(2)
     SubLimits(src+1,3) = nxRecv(3)  
   END DO

   !----- combine the output files from each subdomain
   DO n = 0,(fileCount-1)
      !----- print combining status...
      WRITE(5,*)
      WRITE(5,*) 'combining field output file',n+1,'of',fileCount
      WRITE(5,*) 'reading/deleting...'
      CALL FLUSH(5)
      DO nn = 1,NumSubsTotal
         WRITE(nthSub(1:5),'(I5.5)') nn							! write subdomain number to 'nthSub' for output file exentsions
         !----- open the nnth output file for the nth subdomain 
         WRITE(iter_char(1:7),'(I7.7)') filenum(n)					! write the file number (iteration) to a charater
         OPEN(60,FILE='out-'//iter_char//'-'//nthSub//'.dat')				! open file

         !----- read the output file
         numLines = SubLimits(nn,1)*SubLimits(nn,2)*SubLimits(nn,3)			! determine number of lines to read
         READ(60,*)									! first line is variable info
         READ(60,*)									! second line is zone info
         DO nnn = 1,numLines
            READ(60,*) i,j,k,							&	! i,j,k node location
            FieldData(i,j,k,1),FieldData(i,j,k,2),FieldData(i,j,k,3),		&	! u,v,w @ i,j,k
            FieldData(i,j,k,4),							&	! rho(i,j,k)
            FieldData(i,j,k,5),							&	! phi(i,j,k)
            FieldData(i,j,k,6)								! node(i,j,k)
         END DO	
         CLOSE(60,STATUS='DELETE')							! close and delete current output file (subdomain)
      END DO

      !----- print combining status...
      WRITE(5,*) 'writing...'
      CALL FLUSH(5)

      !----- open and write to new combined file
      OPEN(685,FILE='out-'//iter_char//'.dat')
      WRITE(685,*) 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
      WRITE(685,'(A10,E15.5,A5,I4,A5,I4,A5,I4,A8)') 'ZONE T="',filenum(n)/(nt/nPers),'" I=',nx,' J=',ny,' K=',nz,'F=POINT'
      DO k=1,nz
         DO j=1,ny
            DO i=1,nx
               WRITE (685,'(3I6,5E15.5,I6)') i,j,k,				&	! x,y,z node location
                     FieldData(i,j,k,1),FieldData(i,j,k,2),FieldData(i,j,k,3),	&	! u,v,w @ i,j,k
                     FieldData(i,j,k,4),					&	! rho(i,j,k)
                     FieldData(i,j,k,5),					&	! phi(i,j,k)
                     INT(FieldData(i,j,k,6))						! node(i,j,k)									
            END DO
         END DO
      END DO
      CLOSE(685)									! close current output file (combined)
   END DO

   !----- End timer and print the amount of time it took for the combining
   CALL SYSTEM_CLOCK(combine2,rate)							! End the Timer
   WRITE(5,*)
   WRITE(5,*)
   WRITE(5,*) 'Total Time to Combine Files (min.):', ((combine2-combine1)/REAL(rate))/60.0_dbl
   WRITE(5,*)
   WRITE(5,*)
   CLOSE(5)
   DEALLOCATE(FieldData)
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)							! synchronize all processing units before next loop [Intrinsic]

!===================================================================================================
END SUBROUTINE MergeFields
!===================================================================================================










!===================================================================================================
SUBROUTINE MergeMass					! combines the subdomain output into an output file for the entire computational domain 
!===================================================================================================
IMPLICIT NONE

REAL(dbl), ALLOCATABLE	:: MassData(:,:,:)		! mass data from each subdomain, stored to be rearranged for combined output
REAL(dbl)    :: mass1,mass2				! mass sums
INTEGER(lng) :: i,n,nn					! iteration, loop variables
INTEGER(lng) :: numLines				! number of lines to read
INTEGER(lng) :: combine1,combine2			! clock variables
INTEGER(lng) :: mpierr					! MPI standard error variable
CHARACTER(5) :: nthSub					! current subdomain stored as a character
REAL(dbl), ALLOCATABLE	:: timeread(:)			! time

IF(myid .EQ. master) THEN
  numLines = (nt-iter0) + 1_lng
  ALLOCATE(MassData(numLines,4,NumSubsTotal))
  ALLOCATE(timeread(numLines))

  ! initialize MassData to 0
  MassData = 0.0_dbl

  ! print combining status...
  CALL SYSTEM_CLOCK(combine1,rate)					! Restart the Timer
  OPEN(5,FILE='status.dat',POSITION='APPEND')	
  WRITE(5,*)
  WRITE(5,*)
  WRITE(5,*) 'Combining mass output files and deleting originials...'
  WRITE(5,*)     
  CALL FLUSH(5)

  DO n = 1,NumSubsTotal
    ! print combining status...
    WRITE(5,*)
    WRITE(5,*) 'combining mass output file',n,'of',NumSubsTotal
    WRITE(5,*) 'reading/deleting...'
    CALL FLUSH(5)
    WRITE(nthSub(1:5),'(I5.5)') n					! write subdomain number to 'nthSub' for output file exentsions

    ! open the output file from the nth subdomain 
    OPEN(2458,FILE='mass-'//nthSub//'.dat')				! open file

    ! read the output file
    READ(2458,*)							! first line is variable info
    READ(2458,*)							! second line is zone info
    DO nn = 1,numLines
      READ(2458,*) i,timeread(nn),MassData(nn,1,n),	&		
                     MassData(nn,2,n),MassData(nn,3,n),MassData(nn,4,n)
    END DO
    CLOSE(2458,STATUS='DELETE')						! close and delete current output file (subdomain)

  END DO

  ! print combining status...
  WRITE(5,*) 'writing...'
  CALL FLUSH(5)

  ! open and write to new combined file
  OPEN(2459,FILE='mass.dat')
  WRITE(2459,'(A81)') '#VARIABLES = "period","time" "mass_actual", "mass_theoretical","fmovingsum  ","fmovingrhosum  "'
  WRITE(2459,*) '#ZONE F=POINT'
  DO nn=1,numLines

    ! initialize the summations
    mass1 = 0.0_dbl
    mass2 = 0.0_dbl 

    DO n=1,NumSubsTotal
      mass1 = mass1 + MassData(nn,1,n)
      mass2 = mass2 + MassData(nn,2,n)
    END DO

    i = (iter0-1_lng) + nn
    WRITE(2459,'(3E15.5)') REAL(i/(nt/nPers)), mass1, mass2	
       
  END DO
  CLOSE(2459)																													! close current output file (combined)

  ! End timer and print the amount of time it took for the combining
  CALL SYSTEM_CLOCK(combine2,rate)																						! End the Timer
  WRITE(5,*)
  WRITE(5,*)
  WRITE(5,*) 'Time to Combine Files (min.):', ((combine2-combine1)/REAL(rate))/60.0_dbl
  CLOSE(5)

  DEALLOCATE(MassData)
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)																				! synchronize all processing units before next loop [Intrinsic]

!===================================================================================================
END SUBROUTINE MergeMass
!===================================================================================================








!===================================================================================================
SUBROUTINE MergeScalar	! combines the subdomain output into an output file for the entire computational domain 
!===================================================================================================
IMPLICIT NONE

REAL(dbl), ALLOCATABLE    :: ScalarData(:,:,:)			! scalar data from each subdomain, stored to be rearranged for combined output
REAL(dbl)		:: phiAbsTotal(iter0:nt)		! storage of total aborbed scalar for calulation of absorption rate
REAL(dbl)		:: phiAbsTotalS(iter0:nt)		! storage of total aborbed scalar for calulation of absorption rate (outer surface)
REAL(dbl)		:: phiAbsTotalV(iter0:nt)		! storage of total aborbed scalar for calulation of absorption rate (villi)
REAL(dbl)		:: phi1,phi2,phi3,phi4,phi5, phi6,phi7	! scalar sums
REAL(dbl)		:: SAtime,SA(iter0:nt)			! time from surface area file, surface area
REAL(dbl)		:: phiAverage(iter0:nt)			! period,surface area, average flux, bulk scalar concentration, diffusion resistance, USL thickness (2 values of phi*)
INTEGER(lng)	:: i,n,nn					! loop variables
INTEGER(lng)	:: numLines					! number of lines to read
INTEGER(lng)	:: combine1,combine2				! clock variables
INTEGER(lng)	:: mpierr					! MPI standard error variable
CHARACTER(5)	:: nthSub					! current subdomain stored as a character

IF(myid .EQ. master) THEN

  numLines = (nt-iter0) + 1_lng

  ALLOCATE(ScalarData(numLines,7,NumSubsTotal))

  ! initialize ScalarData and phiAbsTotal to 0
  ScalarData = 0.0_dbl
  phiAbsTotal = 0.0_dbl

  ! print combining status...
  CALL SYSTEM_CLOCK(combine1,rate)							! Restart the Timer
  OPEN(5,FILE='status.dat',POSITION='APPEND')						
  WRITE(5,*)
  WRITE(5,*)
  WRITE(5,*) 'Combining scalar output files and deleting originials...'
  WRITE(5,*)     
  CALL FLUSH(5)

  DO n = 1,NumSubsTotal

    ! print combining status...
    WRITE(5,*)
    WRITE(5,*) 'combining scalar output file',n,'of',NumSubsTotal
    WRITE(5,*) 'reading/deleting...'
    CALL FLUSH(5)

    WRITE(nthSub(1:5),'(I5.5)') n							! write subdomain number to 'nthSub' for output file exentsions

    ! open the output file from the nth subdomain 
    OPEN(2472,FILE='scalar-'//nthSub//'.dat')						! open file

    ! read the output file
    READ(2472,*)									! first line is variable info
    READ(2472,*)									! second line is zone info
    DO nn=1,numLines

      READ(2472,*) i,ScalarData(nn,1,n),	&
                     ScalarData(nn,6,n),	&
                     ScalarData(nn,7,n),	&
                     ScalarData(nn,2,n),	&
                     ScalarData(nn,3,n),	&
                     ScalarData(nn,4,n),	&
                     ScalarData(nn,5,n)

!      WRITE(6678,*) nn, i

    END DO
    CLOSE(2472,STATUS='DELETE')								! close and delete current output file (subdomain)

  END DO

  ! print combining status...
  WRITE(5,*) 'writing...'
  CALL FLUSH(5)

  ! open and write to files
  OPEN(2473,FILE='scalar.dat')
  WRITE(2473,'(A100)') 'VARIABLES = "period", "phiA", "phiAS", "phiAV", "phiT-phiD", "phiD", "phA+phiD", "phiAverage"'
  WRITE(2473,*) 'ZONE F=POINT'

!  OPEN(2474,FILE='SA.dat')
!  READ(2474,*)					! first line is variable info
!  READ(2474,*)					! second line is zone info

  DO nn=1,numLines

    ! initialize the summations
    phi1 = 0.0_dbl
    phi2 = 0.0_dbl 
    phi3 = 0.0_dbl
    phi4 = 0.0_dbl
    phi5 = 0.0_dbl
    phi6 = 0.0_dbl
    phi7 = 0.0_dbl

    DO n=1,NumSubsTotal
      phi1 = phi1 + ScalarData(nn,1,n)
      phi2 = phi2 + ScalarData(nn,2,n)
      phi3 = phi3 + ScalarData(nn,3,n)
      phi4 = phi4 + ScalarData(nn,4,n)
      phi5 = phi5 + ScalarData(nn,5,n)
      phi6 = phi6 + ScalarData(nn,6,n)
      phi7 = phi7 + ScalarData(nn,7,n)
    END DO

    i = (iter0-1_lng) + nn

    WRITE(2473,'(8E25.15)') REAL(i/(nt/nPers)), phi1, phi6, phi7, phi2, phi3, phi4, phi5/NumSubsTotal		! write to combined output file

    phiAbsTotal(i) = phi1    											! store phiAbsorbed for calculation of absorption rate (total)
    phiAbsTotalS(i) = phi6    											! store phiAbsorbed for calculation of absorption rate (outer surface)
    phiAbsTotalV(i) = phi7    											! store phiAbsorbed for calculation of absorption rate (villi)
    phiAverage(i) = phi5/NumSubsTotal										! store phiAverage for calculation of USL 												 		

!    READ(2474,*) SAtime, SA(i)											! read in surface area from file for calculation of scalar flux/USL

  END DO

  CLOSE(2473)													! close output file (combined)
! CLOSE(2474)													! close the surface area file

  ! End timer and print the amount of time it took for the combining
  CALL SYSTEM_CLOCK(combine2,rate)										! End the Timer
  WRITE(5,*)
  WRITE(5,*)
  WRITE(5,*) 'Time to Combine Files (min.):', ((combine2-combine1)/REAL(rate))/60.0_dbl
  CLOSE(5)

  IF(nt .GT. phiStart) THEN
    CALL PrintAbsRate(phiAbsTotal,phiAbsTotalS,phiAbsTotalV,phiAverage,SA)					! calculate and output the absorption rate
  END IF

  DEALLOCATE(ScalarData)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)										! synchronize all processing units before next loop [Intrinsic]

!===================================================================================================
END SUBROUTINE MergeScalar
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintAbsRate(phiAbsTotal,phiAbsTotalS,phiAbsTotalV,phiAverage,SA)	! calculate and output the absorption rate
!===================================================================================================
IMPLICIT NONE

REAL(dbl), INTENT(IN)	:: phiAbsTotal(iter0:nt)	! total aborbed scalar
REAL(dbl), INTENT(IN)	:: phiAbsTotalS(iter0:nt)	! total aborbed scalar (outer surface)
REAL(dbl), INTENT(IN)	:: phiAbsTotalV(iter0:nt)	! total aborbed scalar (villi)
REAL(dbl), INTENT(IN)	:: phiAverage(iter0:nt)		! average scalar in the domain
REAL(dbl), INTENT(IN)	:: SA(iter0:nt)			! surface area
REAL(dbl):: SAS, SAV					! surface area (surface), surface area (villi)
REAL(dbl):: AbsRate, AbsRateS, AbsRateV			! absorption rate at the nth time step
REAL(dbl):: Js,JsS,JsV					! average flux, average flux (outer surface), average flux (villi)
REAL(dbl):: phiBulk,Rw1,Rw2,USL1,USL2			! average flux,bulk scalar concentration,diffusion resistance, USL thicknesses (2 values of phi*)
INTEGER(lng):: n					! loop variable

phiBulk = 1.0_dbl 					! Define the nominal bulk concentration

SAV = numVilliActual*((2.0_dbl*PI*Rv*(Lv-Rv)) + (2.0_dbl*PI*Rv*Rv))	! surface area of the villi

!----- set up output file
OPEN(2479,FILE='phiRate.dat')
WRITE(2479,'(A100)') 'VARIABLES = "period", "AbsRate", "AbsRateS", "AbsRateV","flux", "fluxS", "fluxV" "usl1", "usl2"'
WRITE(2479,*) 'ZONE F=POINT'
IF(restart) THEN

  !----- iter0 (1st order forward differencing)
  AbsRate = (phiAbsTotal(iter0+1) - phiAbsTotal(iter0))/tcf
  AbsRateS = (phiAbsTotalS(iter0+1) - phiAbsTotalS(iter0))/tcf
  AbsRateV = (phiAbsTotalV(iter0+1) - phiAbsTotalV(iter0))/tcf
  Js = AbsRate/SA(iter0)
  JsS = AbsRateS/(SA(iter0)-SAV)
  JsV = AbsRateV/SAV

  !----- Calculate the resistance to diffusion
  IF(Js .GT. 1e-18) THEN
    Rw1 = phiBulk/Js	
    Rw2 = phiAverage(iter0)/Js	
  ELSE
    Rw1 = 0.0_dbl										! set to 0.0 while the scalar is diffusing to the surface (Js=0)
    Rw2 = 0.0_dbl	
  END IF

  USL1 = Rw1*Dm*Dmcf										! calculate the effective UWL thickness
  USL2 = Rw2*Dm*Dmcf
  
  WRITE(2479,'(9E25.15)') iter0/(nt/nPers), AbsRate, AbsRateS, AbsRateV, Js, JsS, JsV, USL1, USL2

  !----- iter0+1 to nt-1 (2nd order central differencing) 
  DO n=iter0+1,nt-1

    AbsRate = (phiAbsTotal(n+1) - phiAbsTotal(n-1))/(2.0_dbl*tcf)				
    AbsRateS = (phiAbsTotalS(n+1) - phiAbsTotalS(n-1))/(2.0_dbl*tcf)				
    AbsRateV = (phiAbsTotalV(n+1) - phiAbsTotalV(n-1))/(2.0_dbl*tcf)				
    Js = AbsRate/SA(n)
    JsS = AbsRateS/(SA(n)-SAV)
    JsV = AbsRateV/SAV

    !----- Calculate the resistance to diffusion
    IF(Js .GT. 1e-18) THEN
      Rw1 = phiBulk/Js	
      Rw2 = phiAverage(n)/Js	
    ELSE
      Rw1 = 0.0_dbl										! set to 0.0 while the scalar is diffusing to the surface (Js=0)
      Rw2 = 0.0_dbl	
    END IF

    USL1 = Rw1*Dm*Dmcf										! calculate the effective UWL thickness
    USL2 = Rw2*Dm*Dmcf
  
    WRITE(2479,'(9E25.15)') n/(nt/nPers), AbsRate, AbsRateS, AbsRateV, Js, JsS, JsV, USL1, USL2

  END DO

ELSE

  ! phiStart (1st order forward differencing)
  AbsRate = (phiAbsTotal(phiStart+1) - phiAbsTotal(phiStart))/tcf
  AbsRateS = (phiAbsTotalS(phiStart+1) - phiAbsTotalS(phiStart))/tcf
  AbsRateV = (phiAbsTotalV(phiStart+1) - phiAbsTotalV(phiStart))/tcf
  Js = AbsRate/SA(phiStart)
  JsS = AbsRateS/(SA(phiStart)-SAV)
  JsV = AbsRateV/SAV

  ! Calculate the resistance to diffusion
  IF(Js .GT. 1e-18) THEN
    Rw1 = phiBulk/Js	
    Rw2 = phiAverage(phiStart)/Js	
  ELSE
    Rw1 = 0.0_dbl										! set to 0.0 while the scalar is diffusing to the surface (Js=0)
    Rw2 = 0.0_dbl	
  END IF

  USL1 = Rw1*Dm*Dmcf										! calculate the effective UWL thickness
  USL2 = Rw2*Dm*Dmcf

  WRITE(2479,'(9E25.15)') phiStart/(nt/nPers), AbsRate, AbsRateS, AbsRateV, Js, JsS, JsV, USL1, USL2

  ! phiStart+1 to nt-1 (2nd order central differencing) 
  DO n=phiStart+1,nt-1

    AbsRate = (phiAbsTotal(n+1) - phiAbsTotal(n-1))/(2.0_dbl*tcf)				
    AbsRateS = (phiAbsTotalS(n+1) - phiAbsTotalS(n-1))/(2.0_dbl*tcf)	
    AbsRateV = (phiAbsTotalV(n+1) - phiAbsTotalV(n-1))/(2.0_dbl*tcf)	
    Js = AbsRate/SA(n)
    JsS = AbsRateS/(SA(n)-SAV)
    JsV = AbsRateV/SAV

    ! Calculate the resistance to diffusion
    IF(Js .GT. 1e-18) THEN
      Rw1 = phiBulk/Js	
      Rw2 = phiAverage(n)/Js	
    ELSE
      Rw1 = 0.0_dbl										! set to 0.0 while the scalar is diffusing to the surface (Js=0)
      Rw2 = 0.0_dbl	
    END IF

    USL1 = Rw1*Dm*Dmcf										! calculate the effective UWL thickness
    USL2 = Rw2*Dm*Dmcf
  
    WRITE(2479,'(9E25.15)') n/(nt/nPers), AbsRate, AbsRateS, AbsRateV, Js, JsS, JsV, USL1, USL2

  END DO

END IF

!---- nt (1st order backward differencing)
AbsRate = (phiAbsTotal(nt) - phiAbsTotal(nt-1))/tcf
AbsRateS = (phiAbsTotalS(nt) - phiAbsTotalS(nt-1))/tcf
AbsRateV = (phiAbsTotalV(nt) - phiAbsTotalV(nt-1))/tcf
Js = AbsRate/SA(nt)
JsS = AbsRateS/(SA(nt)-SAV)
JsV = AbsRateV/SAV

!---- Calculate the resistance to diffusion
IF(Js .GT. 1e-18) THEN
  Rw1 = phiBulk/Js	
  Rw2 = phiAverage(nt)/Js	
ELSE
  Rw1 = 0.0_dbl											! set to 0.0 while the scalar is diffusing to the surface (Js=0)
  Rw2 = 0.0_dbl	
END IF

USL1 = Rw1*Dm*Dmcf										! calculate the effective UWL thickness
USL2 = Rw2*Dm*Dmcf

WRITE(2479,'(9E25.15)') nt/(nt/nPers), AbsRate, AbsRateS, AbsRateV, Js, JsS, JsV, USL1, USL2

CLOSE(2479)

!===================================================================================================
END SUBROUTINE PrintAbsRate
!===================================================================================================










!===================================================================================================
END MODULE Output
!===================================================================================================
