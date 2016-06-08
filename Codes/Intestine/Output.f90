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

   !----- Status-----
   OPEN(5,FILE='status.dat')										
   CALL FLUSH(5)													

  !----- Surface Area-----
  OPEN(2474,FILE='SA.dat',POSITION='APPEND')
  WRITE(2474,'(A36)') 'VARIABLES = "period", "SA"'
  WRITE(2474,*) 'ZONE F=POINT'
  CALL FLUSH(2474)

  !----- Walll Flux-----
  OPEN(4748,FILE='wall_flux.dat')
  WRITE(4748,*) 'VARIABLES = "Axial Distance", "Flux"'
  CALL FLUSH(4748)

  !---- Volume -----
  OPEN(2460,FILE='volume.dat')
  WRITE(2460,*) 'VARIABLES = "period", "volume"'
  WRITE(2460,*) 'ZONE F=POINT'
  CALL FLUSH(2460)
END IF

!----- Mass
OPEN(2458,FILE='mass-'//sub//'.dat')
WRITE(2458,'(A120)') '#VARIABLES = period, time, mass_theory, mass_actual, mass_err'
WRITE(2458,*) '#ZONE F=POINT'
CALL FLUSH(2458)

!----- Drug Conservation
OPEN(2472,FILE='Drug-Conservation-'//sub//'.dat')
WRITE(2472,'(A120)') '#VARIABLES =iter,time, Drug_Initial, Drug_Released_Total, Drug_Absorbed, Drug_Remained_in_Domain, Drug_Loss_Percent, Drug_Loss_Modified_Percent'
WRITE(2472,*) '#ZONE F=POINT'
CALL FLUSH(2472)

!---- Test Output
OPEN(9,FILE='testoutput-'//sub//'.dat',POSITION='append')
CALL FLUSH(9)

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
  CLOSE(2474)			! Surface Area
  CLOSE(2460)			! Volume
END IF
CLOSE(2458)			! Mass
CLOSE(2472)			! Scalar
CLOSE(9) 			! Test Output
!===================================================================================================
END SUBROUTINE CloseOutputFiles
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintStatus	! Writes the current Status to a File
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
SUBROUTINE PrintPeriodicRestart		! prints restart file periodically to guard against a crash
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m		! index variables

!----- Write restart file (restart.XX) and corresponding starting iteration file (iter0.dat)
IF (MOD(iter,100) .EQ. 0) THEN
   !----- 1st Attempt
   OPEN(50,FILE='restart.'//sub)
   DO k=0,nzSub+1
      DO j=0,nySub+1
         DO i=0,nxSub+1
            WRITE(50,'(20E15.5)') u(i,j,k),v(i,j,k),w(i,j,k),rho(i,j,k),phi(i,j,k),(f(m,i,j,k),m=0,14)
         END DO
      END DO
   END DO

   CLOSE(50)

   IF (myid .EQ. master) THEN
      OPEN(55,FILE='iter0.dat')
      WRITE(55,*) iter+1
      CLOSE(55)
   END IF

   !----- 2nd Attempt (Backup in case of computer error during simulation execution)
   OPEN(51,FILE='restart-bak.'//sub)
   DO k=0,nzSub+1
      DO j=0,nySub+1
         DO i=0,nxSub+1
            WRITE(51,'(20E15.5)') u(i,j,k),v(i,j,k),w(i,j,k),rho(i,j,k),phi(i,j,k),(f(m,i,j,k),m=0,14)
         END DO
      END DO
   END DO

   CLOSE(51)

   IF (myid .EQ. master) THEN
      OPEN(56,FILE='iter0-bak.dat')
      WRITE(56,*) iter+1
      CLOSE(56)
   END IF
END IF

!===================================================================================================
END SUBROUTINE PrintPeriodicRestart
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintFinalRestart		! prints restart file periodically to guard against a crash
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m			! index variables

!----- Write restart file (restart.XX) and corresponding starting iteration file (iter0.dat)
OPEN(500,FILE='restart.'//sub)

DO k=0,nzSub+1
   DO j=0,nySub+1
      DO i=0,nxSub+1
         WRITE(500,*) node(i,j,k)
         WRITE(500,*) u(i,j,k)
         WRITE(500,*) v(i,j,k)
         WRITE(500,*) w(i,j,k)
         WRITE(500,*) rho(i,j,k)
         WRITE(500,*) phi(i,j,k)
         DO m=0,NumDistDirs
            WRITE(500,*) f(m,i,j,k)
         END DO
      END DO
   END DO
END DO

WRITE(500,*) phiAbsorbed
WRITE(500,*) phiAbsorbedS
WRITE(500,*) phiAbsorbedV
WRITE(500,*) phiInOut

CLOSE(500)

IF (myid .EQ. master) THEN
   OPEN(550,FILE='iter0.dat')
   WRITE(550,*) iter
   CLOSE(550)
END IF
!===================================================================================================
END SUBROUTINE PrintFinalRestart
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
   current => ParListHead%next					 ! Using linked lists

   DO WHILE (ASSOCIATED(current))
      numParticlesSub = numParticlesSub + 1_lng
      next => current%next 					! copy pointer of next node
      WRITE(160,1001)   current%pardata%cur_part  ,',', &
		 	current%pardata%xp 	  ,',',	&
			current%pardata%yp  	  ,',',	&
			current%pardata%zp 	  ,',',	&
                        current%pardata%up*vcf 	  ,',',	&
		       	current%pardata%vp*vcf 	  ,',',	&
			current%pardata%wp*vcf 	  ,',',	&
                        current%pardata%parid 	  ,',',	&
			current%pardata%sh 	  ,',',	&
			current%pardata%rp/xcf 	  ,',',	&
			current%pardata%bulk_conc/Cs_mol  ,',', &
			current%pardata%delNBbyCV ,',', &
			current%pardata%Sst 	  ,',',	&
			current%pardata%S 	  ,',',	&
			current%pardata%Veff 	  ,',',	&
			current%pardata%Nbj
1001 format (I4,a2,F9.4,a2,F9.4,a2,F9.4,a2,F10.6,a2,F10.6,a2,F10.6,a2,I5,a2,F12.8,a2,F15.10,a2,F15.7,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2)
     current => next   						! point to next node in the list
  ENDDO

  CLOSE(160)
  numparticleSubfile(parfileCount) = numParticlesSub	
  parfileCount = parfileCount + 1_lng
ENDIF


!===================================================================================================
END SUBROUTINE PrintParticles	! print particle position, velocity, radius, and concentrationto output files
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
SUBROUTINE CheckVariables	! checks to see where variables become "NaN" 
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k		! index variables
INTEGER(lng) :: uInt	   	! integer version of u
REAL(dbl) :: localu

DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         localu = u(i,j,k)
         uInt = INT(u(i,j,k))
         IF (uInt .GT. 1000_lng) THEN
            OPEN(1042,FILE='u-'//sub//'.dat')
            WRITE(1042,*) 'iter=', iter
            WRITE(1042,*) 'i=', i, 'j=', j, 'k=', k
            WRITE(1042,*) 'node=', node(i,j,k)
            WRITE(1042,*) 'u=', u(i,j,k), 'rho=', rho(i,j,k)
            WRITE(1042,*) 'node=', node(i,j,k)
            CLOSE(1042)
            STOP
         END IF 
      END DO
   END DO
END DO

!===================================================================================================
END SUBROUTINE CheckVariables
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

INTEGER(lng) :: i,j,k					! index variables
INTEGER(lng) :: numFluids				! number of fluid nodes in the domain
REAL(dbl)    :: phiDomain, phiIC, Drug_Initial		! current amount of scalar in the domain
REAL(dbl)    :: phiAverage				! average scalar in the domain
REAL(dbl)    :: zcf3					! node volume in physical units
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

CALL ScalarInOut    	 				! Calculate the amount of scalar that entered/left through the inlet/outlet

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
IF (ParticleTrack.EQ.ParticleOn .AND. iter .GE. phiStart) THEN
   current => ParListHead%next
   DO WHILE (ASSOCIATED(current))
      next => current%next
      Drug_Released_Total = Drug_Released_Total + current%pardata%delNBbyCV * zcf3
      current => next
   ENDDO
END IF

Drug_Initial=  phiTotal *zcf3
Drug_Absorbed = phiAbsorbed * zcf3
Drug_Remained_in_Domain = phiDomain * zcf3
Drug_Loss = (Drug_Released_Total + Drug_Initial) - (Drug_Absorbed + Drug_Remained_in_Domain)  
Drug_Loss_Modified = (Drug_Released_Total+ Drug_Initial- Negative_phi_Total) - (Drug_Absorbed + Drug_Remained_in_Domain)

IF (Drug_Released_Total .LT. 1e-20) THEN
   Drug_Released_Total =1e-20
END IF

Drug_Loss_Percent = (Drug_Loss / (Drug_Released_Total+Drug_Initial)) * 100.0_lng
Drug_Loss_Modified_Percent = (Drug_Loss_Modified / (Drug_Released_Total+Drug_Initial)) * 100.0_lng  

IF (abs(Drug_Absorbed) .lt. 1.0e-40) THEN
   Drug_Absorbed = 0.0_lng
ENDIF

WRITE(2472,'(I7, F9.3, 6E21.13)') iter, iter*tcf, Drug_Initial, Drug_Released_Total, Drug_Absorbed, Drug_Remained_in_Domain, Drug_Loss_Percent, Drug_Loss_Modified_Percent 
CALL FLUSH(2472)
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
!CALL MergeMass
IF (ParticleTrack.EQ.ParticleOn .AND. iter .GE. phiStart) THEN 	! If particle tracking is 'on' then do the following
   CALL MergeParticleOutput					! Merge particle output
ENDIF
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
SUBROUTINE MergeParticleOutput				! combines the subdomain particle output into an output file for the entire computational domain 
!===================================================================================================
IMPLICIT NONE

REAL(dbl), ALLOCATABLE	:: ParticleData(:,:)		! data for each particle in the computational domain
INTEGER(lng) :: ParticleDistribution(0:nt,NumSubsTotal) ! array containing numparticlesSub for each sub domain
INTEGER(lng) :: nxSend(3),nxRecv(3),npRecv,npSend	! nx-,ny-,nzSub from each subdomain
INTEGER(lng) :: ii,jj,kk				! neighboring indices for writing
INTEGER(lng) :: i,j,k,n,nn,nnn				! loop variables
INTEGER(lng) :: dest, src, tag				! send/recv variables: destination, source, message tag
INTEGER(lng) :: stat(MPI_STATUS_SIZE)			! status object: rank and tag of the sending processing unit/message
INTEGER(lng) :: mpierr					! MPI standard error variable
INTEGER(lng) :: numLines,numLineStart,numLineEnd	! number of lines to read
INTEGER(lng) :: combine1,combine2			! clock variables
CHARACTER(7) :: iter_char				! iteration stored as a character 
CHARACTER(5) :: nthSub					! current subdomain stored as a character
CHARACTER(2) :: tmpchar

!----- send subdomain information to the master
tag = 63_lng											! starting message tag (arbitrary)

IF (myid .NE. master) THEN
   dest	= master										! send to master
   CALL MPI_SEND(numparticleSubfile(0:nt),nt+1,MPI_INTEGER,dest,tag,MPI_COMM_WORLD,mpierr)	! send num particles in each subdomain
ELSE
   ALLOCATE(ParticleData(numparticlesDomain,15))
   !----- initialize FieldData to 0
   ParticleData = 0.0_dbl

   !-----  print combining status...
   CALL SYSTEM_CLOCK(combine1,rate)								! Restart the Timer
   OPEN(5,FILE='status.dat',POSITION='APPEND')										
   WRITE(5,*)
   WRITE(5,*)
   WRITE(5,*) 'Combining particle output files and deleting originials...'
   WRITE(5,*)     
   CALL FLUSH(5)

   ParticleDistribution(0:nt,1) = numparticleSubfile(0:nt)

   DO src = 1,(numprocs-1)
      CALL MPI_RECV(ParticleDistribution(0:nt,src+1),nt+1,MPI_INTEGER,src,tag,MPI_COMM_WORLD,stat,mpierr) ! receive np from each subdomain for each output
   END DO

   !----- combine the output files from each subdomain
   DO n= 0,(parfileCount-1)
      !----- print combining status...
      WRITE(5,*)
      WRITE(5,*) 'combining particle output file',n+1,'of',parfileCount
      WRITE(5,*) 'reading/deleting...'
      CALL FLUSH(5)
      numLineStart=  0_lng
      numLineEnd  =  0_lng
      DO nn = 1,NumSubsTotal
  	 WRITE(nthSub(1:5),'(I5.5)') nn	! write subdomain number to 'nthSub' for output file exentsions

         !----- open the nnth output file for the nth subdomain 
	 WRITE(iter_char(1:7),'(I7.7)') parfilenum(n)			! write the file number (iteration) to a charater
	 OPEN(60,FILE='pardat-'//iter_char//'-'//nthSub//'.csv')	! open file

	 !----- read the output file
  	 numLines = ParticleDistribution(n,nn)				! determine number of lines to read
	 READ(60,*)							! first line is variable info
			
	 numLineStart	=  numLineEnd 
	 numLineEnd	=  numLineStart + numLines
	 DO nnn = numLineStart+1_lng,numLineEnd
       	    READ(60,*) 	ParticleData(nnn,1),	& 
	    		ParticleData(nnn,2),	&	
	    		ParticleData(nnn,3),	&
			ParticleData(nnn,4),	&
			ParticleData(nnn,5),	&
			ParticleData(nnn,6),	&
			ParticleData(nnn,7),	&
			ParticleData(nnn,8),	&
			ParticleData(nnn,9),	&
			ParticleData(nnn,10),	&
			ParticleData(nnn,11),	&
			ParticleData(nnn,12),	&	
			ParticleData(nnn,13),	&
			ParticleData(nnn,14),	&
			ParticleData(nnn,15)
	 END DO
	 CLOSE(60,STATUS='DELETE') ! close and delete current output file (subdomain)
       END DO

       !----- print combining status...
       WRITE(5,*) 'writing...'
       CALL FLUSH(5)

       !----- Write to combined output file	
       OPEN(685,FILE='pardat-all-'//iter_char//'.csv')
       WRITE(685,*) '"x","y","z","u","v","w","ParID","Sh","rp","bulk_conc","delNBbyCV","Sst","S","Veff","Nbj"'

       DO nnn=1,numLineEnd
          WRITE(685,1001) ParticleData(nnn,1)  		,',', &
		 	  ParticleData(nnn,2)		,',', &
			  ParticleData(nnn,3) 		,',', &
			  ParticleData(nnn,4)		,',', &
			  ParticleData(nnn,5)		,',', &
			  ParticleData(nnn,6) 		,',', &
			  INT(ParticleData(nnn,7),lng)	,',', &
			  ParticleData(nnn,8)		,',', &
			  ParticleData(nnn,9) 		,',', &
			  ParticleData(nnn,10)		,',', &
			  ParticleData(nnn,11)		,',', &
			  ParticleData(nnn,12) 		,',', &
			  ParticleData(nnn,13)		,',', &
			  ParticleData(nnn,14)		,',', &
			  ParticleData(nnn,15)

1001 format (a3, F12.5,a2,F12.5,a2,F12.5,a2,F12.5,a2,F12.5,a2,F15.5,a2,1I5,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2,F15.10,a2)
	
	END DO
	CLOSE(685)						! close current output file (combined)

    END DO
    !----- End timer and print the amount of time it took for the combining
    CALL SYSTEM_CLOCK(combine2,rate)				! End the Timer
    WRITE(5,*)
    WRITE(5,*)
    WRITE(5,*) 'Total Time to Combine Particle Files (min.):', ((combine2-combine1)/REAL(rate))/60.0_dbl
    WRITE(5,*)
    WRITE(5,*)
    CLOSE(5)

    DEALLOCATE(ParticleData)
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)				! synchronize all processing units before next loop [Intrinsic]

!===================================================================================================
END SUBROUTINE MergeParticleOutput
!===================================================================================================








!===================================================================================================
SUBROUTINE StreamFunction(vz,nodeFlag,sf)		! integrates the velocity to get the streamfunction
!===================================================================================================
IMPLICIT NONE

REAL(dbl), INTENT(IN) 	:: vz(nx,ny,nz)			! w(:,:,:) (axial velocity)
INTEGER(lng), INTENT(IN):: nodeFlag(nx,ny,nz)		! node(:,:,:) (node flags)
REAL(dbl), INTENT(OUT) 	:: sf(nz,nx)			! streamfunction
REAL(dbl) 		:: vz2(nz,nx)			! plane velocity
INTEGER(lng) 	:: nodeFlag2(nz,nx)			! plane node flags
INTEGER(lng)	:: i,k					! index variables

! initialize the streamfunction array
sf = 0.0_dbl

! check to see if ny is odd or even, then integrate the velocity to get the stream function
IF(MOD(nx,2) .NE. 0) THEN ! (ny is odd - nodes exist on center plane)

  DO k=1,nz
    DO i=1,nx
      ! store the plane velocity and node flags
      vz2(k,i) = vz(i,1,k)
      nodeFlag2(k,i) = nodeFlag(i,1,k)
    END DO
    ! set the centerline streamfunction
    sf(k,1) = 0.0_dbl
  END DO

  ! integrate to obtain streamfunction
  DO k=1,nz
    DO i=2,nx
      IF(nodeFlag2(k,i) .EQ. FLUID) THEN
        sf(k,i)			= 0.5_dbl*(xx(i)*vz2(k,i) + xx(i-1)*vz2(k,i-1))*(xx(i)-xx(i-1)) + sf(k,i-1)
      ELSE
        EXIT
      END IF
    END DO
  END DO

ELSE ! (ny is even - no nodes on center plane)

  DO k=1,nz
    DO i=1,nx
      ! store the plane velocity and node flags
      vz2(k,i) = vz(i,1,k)
      nodeFlag2(k,i) = nodeFlag(i,1,k)
    END DO
    ! set the streamfunction at the nodes adjacent to the centerline
    sf(k,1) 	= 0.5_dbl*(xx(1)*vz2(k,1))*(0.5_dbl*xcf)
   END DO

  ! integrate to obtain streamfunction
  DO k=1,nz
    DO i=2,nx
      IF(nodeFlag2(k,i) .EQ. FLUID) THEN
        sf(k,i)			= 0.5_dbl*(xx(i)*vz2(k,i) + xx(i-1)*vz2(k,i-1))*(xx(i)-xx(i-1)) + sf(k,i-1)
      ELSE
        EXIT
      END IF
    END DO
  END DO

END IF

!===================================================================================================
END SUBROUTINE StreamFunction
!===================================================================================================








!===================================================================================================
SUBROUTINE ScalarFlux(n,phiFld,nodeFlag,fluxFld)		! calculates the scalar flux field and flux at the wall
!===================================================================================================
IMPLICIT NONE

REAL(dbl), INTENT(IN) 		:: phiFld(nx,ny,nz)		! phi(:,:,:)
INTEGER(lng), INTENT(IN) 	:: nodeFlag(nx,ny,nz)		! node(:,:,:) (node flags)
INTEGER(lng), INTENT(IN) 	:: n				! file number [filenum(n)=iter]
REAL(dbl), INTENT(OUT) 		:: fluxFld(nz,nx,3)		! flux field
REAL(dbl) 		:: phi2(nz,nx)				! plane scalar
REAL(dbl)		:: wflux(nz-1)				! wall flux
INTEGER(lng) 	:: nodeFlag2(nz,nx)				! plane node flags
INTEGER(lng)	:: i,k						! index variables

! initialize the flux array
fluxFld = 0.0_dbl

! store the scalar and node flags
DO k=1,nz
  DO i=1,nx

    phi2(k,i) = phiFld(i,1,k)
    nodeFlag2(k,i) = nodeFlag(i,1,k)

  END DO
END DO

! differentiate to get scalar flux in each direction
! inside nodes
DO k=2,nz-1
  DO i=2,nx-1

    IF(nodeFlag2(k,i) .EQ. FLUID) THEN
      fluxFld(k,i,1)	= (phi2(k+1,i) - phi2(k-1,i))/(2.0_dbl*zcf)
      fluxFld(k,i,2)	= (phi2(k,i+1) - phi2(k,i-1))/(2.0_dbl*xcf)
      fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)
    END IF

  END DO
END DO

! left boundary
k=1
DO i=2,nx-1
  IF(nodeFlag2(k,i) .EQ. FLUID) THEN
    fluxFld(k,i,1)	= (phi2(k+1,i) - phi2(k,i))/zcf
    fluxFld(k,i,2)	= (phi2(k,i+1) - phi2(k,i-1))/(2.0_dbl*xcf)
    fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)
  END IF
END DO

! right boundary
k=nz
DO i=2,nx-1
  IF(nodeFlag2(k,i) .EQ. FLUID) THEN
    fluxFld(k,i,1)	= (phi2(k,i) - phi2(k-1,i))/zcf
    fluxFld(k,i,2)	= (phi2(k,i+1) - phi2(k,i-1))/(2.0_dbl*xcf)
    fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)
  END IF
END DO

! top boundary
i=nx
DO k=2,nz-1
  IF(nodeFlag2(k,i) .EQ. FLUID) THEN
    fluxFld(k,i,1)	= (phi2(k+1,i) - phi2(k-1,i))/(2.0_dbl*zcf)
    fluxFld(k,i,2)	= (phi2(k,i) - phi2(k,i-1))/xcf
    fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)
  END IF
END DO

! bottom boundary
i=1
DO k=2,nz-1
  IF(nodeFlag2(k,i) .EQ. FLUID) THEN
    fluxFld(k,i,1)	= (phi2(k+1,i) - phi2(k-1,i))/(2.0_dbl*zcf)
    fluxFld(k,i,2)	= (phi2(k,i+1) - phi2(k,i))/xcf
    fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)
  END IF
END DO

! top left corner
k=1
i=nx
fluxFld(k,i,1)	= (phi2(k+1,i) - phi2(k,i))/zcf
fluxFld(k,i,2)	= (phi2(k,i) - phi2(k,i-1))/xcf
fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)

! top right corner
k=nz
i=nx
fluxFld(k,i,1)	= (phi2(k,i) - phi2(k-1,i))/zcf
fluxFld(k,i,2)	= (phi2(k,i) - phi2(k,i-1))/xcf
fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)

! bottom right corner
k=nz
i=1
fluxFld(k,i,1)	= (phi2(k,i) - phi2(k-1,i))/zcf
fluxFld(k,i,2)	= (phi2(k,i+1) - phi2(k,i))/xcf
fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)

! bottom left corner
k=1
i=1
fluxFld(k,i,1)	= (phi2(k+1,i) - phi2(k,i))/zcf
fluxFld(k,i,2)	= (phi2(k,i+1) - phi2(k,i))/xcf
fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)

! bottom left corner
k=1
i=1
fluxFld(k,i,1)	= (phi2(k+1,i) - phi2(k,i))/zcf
fluxFld(k,i,2)	= (phi2(k,i+1) - phi2(k,i))/xcf
fluxFld(k,i,3)	= SQRT(fluxFld(k,i,1)**2 + fluxFld(k,i,2)**2)

!-- PLANAR WALL FLUX --
IF(filenum(n) .GE. phiStart) THEN

  WRITE(4748,'(A8,E15.5,A4,I4,A8)') 'ZONE T="',filenum(n)/(nt/nPers),'" I=', nz-1,' F=POINT'
  CALL Compute_flux(n,phi2,nodeFlag2,wflux)
  DO k=1,nz-1
    WRITE(4748,'(2E15.5)') zz(k)+0.5_dbl*zcf, wflux(k)
  END DO
  CALL FLUSH(4748)

END IF

!===================================================================================================
END SUBROUTINE ScalarFlux
!===================================================================================================









!===================================================================================================
SUBROUTINE Compute_flux(n,phiPlane,nodePlane,wflux)! measures the flux of scalar through the walls (using Yanxing Wang's extrapolation method)
!===================================================================================================
IMPLICIT NONE

INTEGER(lng), INTENT(IN)	:: n				! file number [filenum(n)=iter]
REAL(dbl), INTENT(IN) 		:: phiPlane(nz,nx) 		! planar scalar field
INTEGER(lng), INTENT(IN)	:: nodePlane(nz,nx)		! planar node flag field
REAL(dbl), INTENT(OUT)		:: wflux(nz-1)			! wall flux
INTEGER(lng) :: kLLb,iLLb,kULb,iULb				! lattice coordinates of the nodes of the cell containing point b (left: lower, upper)
INTEGER(lng) :: kLRb,iLRb,kURb,iURb				! lattice coordinates of the nodes of the cell containing point b (right: lower, upper)
INTEGER(lng) :: kLLc,iLLc,kULc,iULc				! lattice coordinates of the nodes of the cell containing point c (left: lower, upper)
INTEGER(lng) :: kLRc,iLRc,kURc,iURc				! lattice coordinates of the nodes of the cell containing point c (right: lower, upper)
INTEGER(lng) :: numpointsb					! number of FLUID nodes on the cell containing point b
INTEGER(lng) :: numpointsc					! number of FLUID nodes of the cell containing point c
INTEGER(lng) :: k, kp1						! lattice indices
REAL(dbl) :: z0,za,zb,zc,z1					! points for calculating the scalar (z-coordinates)
REAL(dbl) :: x0,xa,xb,xc,x1					! points for calculating the scalar (x-coordinates)
REAL(dbl) :: kb,ib,kc,ic					! lattice coordinates of point b and point c
REAL(dbl) :: SinTheta, CosTheta					! sin and cos of the angle theta
REAL(dbl) :: phiLLb,phiULb,phiLRb,phiURb			! scalar the nodes of the cell containing point b (lower left, upper left, lower right, upper right)
REAL(dbl) :: phiLLc,phiULc,phiLRc,phiURc			! scalar the nodes of the cell containing point c (lower left, upper left, lower right, upper right)
REAL(dbl) :: phib, phic						! scalar at the points b and c
REAL(dbl) :: A, B, C, DD					! coefficients of the bilinear scalar equation
REAL(dbl) :: fluxZb, fluxXb, fluxNb				! fluxes at point b (z,x, and normal directions)
REAL(dbl) :: fluxZc, fluxXc, fluxNc				! fluxes at point c (z,x, and normal directions)
REAL(dbl) :: pointsb(4,3), pointsc(4,3)				! 2D array to store the i,j locations, and scalar value of the fluid points for 

DO k=1,nz-1
   kp1	= k + 1_lng						! k + 1
   !----- LINE CONNECTING THE POINTS
   !----- find the z,x coordinates of r(k) (point 0) and r(k+1) (point 1)
   z0 = zz(k)
   x0 = radius(k,n)
   z1 = zz(kp1)
   x1 = radius(kp1,n)

   !----- find the sin and cos of the angle of the line connecting point 0 with point 1
   SinTheta = (x1-x0)/SQRT((x1-x0)*(x1-x0) + zcf*zcf)
   CosTheta = zcf/SQRT((x1-x0)*(x1-x0) + zcf*zcf)

   !----- MIDPOINT
   !----- find the z,x coordinates at the midpoint between a line connecting point 0 with point 1
   za = 0.5_dbl*(z0+z1)
   xa = 0.5_dbl*(x0+x1)

   !----- FIRST POINT
   !----- find the coordinates of point b (1 grid spacing away from point a)
   zb = za + zcf*SinTheta			
   xb = xa + zcf*(-CosTheta)
 
  ! lattice units
  kb = (zb/zcf) + 0.5_dbl	
  ib = (xb/xcf) + 1.0_dbl

  ! find the coordinates of the points on the cell containing point b
  kLLb = FLOOR(kb)		      	! lower left lattice node
  kLRb = kLLb + 1_lng				! lower right lattice node

  ! periodic bcs
  IF(kLLb .EQ. 0) THEN
    kLLb = nz
  ELSE IF(kLRb .EQ. nz+1) THEN
    kLRb = 1_lng
  END IF

  kULb = kLLb				      	! upper left lattice node
  kURb = kLRb							! upper right lattice node

  iLLb = FLOOR(ib)
  iLRb = iLLb
  iULb = iLLb + 1_lng
  iURb = iULb

  ! find the scalar at the points on the cell containing point b
  phiLLb = phiPlane(kLLb,iLLb)	! lower left
  phiULb = phiPlane(kULb,iULb)	! upper left
  phiURb = phiPlane(kURb,iURb)	! upper right
  phiLRb = phiPlane(kLRb,iLRb)	! lower right

  ! find the number of points on the cell containing point b that are in the fluid domain (non-solid) and store the coordinates
  pointsb 		= 0.0_dbl	! initialize points
  numpointsb 	= 0_lng		! initialize numpointsb
 
  IF(nodePlane(kLLb,iLLb) .EQ. FLUID) THEN	! lower left
    numpointsb = numpointsb + 1_lng
    pointsb(numpointsb,1) = kLLb
    pointsb(numpointsb,2) = iLLb
    pointsb(numpointsb,3) = phiLLb
  END IF

  IF(nodePlane(kULb,iULb) .EQ. FLUID) THEN	! upper left
    numpointsb = numpointsb + 1_lng
    pointsb(numpointsb,1) = kULb
    pointsb(numpointsb,2) = iULb
    pointsb(numpointsb,3) = phiULb
  END IF

  IF(nodePlane(kURb,iURb) .EQ. FLUID) THEN	! upper right
    numpointsb = numpointsb + 1_lng
    pointsb(numpointsb,1) = kURb
    pointsb(numpointsb,2) = iURb
    pointsb(numpointsb,3) = phiURb
  END IF

  IF(nodePlane(kLRb,iLRb) .EQ. FLUID) THEN	! lower right
    numpointsb = numpointsb + 1_lng
    pointsb(numpointsb,1) = kLRb
    pointsb(numpointsb,2) = iLRb
    pointsb(numpointsb,3) = phiLRb
  END IF

  ! find the scalar and fluxes at point b
  IF(numpointsb .EQ. 4) THEN								! four-point bilinear interpolation

    A = phiLLb - phiLRb + phiURb - phiULb
    B = (phiURb - phiULb) - (iLLb + 1_lng)*A
    C = -(kLLb + 1_lng)*A - (phiLRb - phiURb)
    DD = phiULb - kLLb*(iLLb + 1_lng)*A - kLLb*B - (iLLb + 1_lng)*C

    phib 	= A*kb*ib + B*kb + C*ib + DD						! scalar
    fluxZb 	= (Dm*Dmcf*(A*ib + B))/zcf						! flux in the z-direction
    fluxXb 	= (Dm*Dmcf*(A*kb + C))/zcf						! flux in the x-direction
    fluxNb 	= fluxZb*SinTheta + fluxXb*(-CosTheta)					! flux in the n-direction

  ELSE IF(numpointsb .EQ. 3) THEN							! three-point bilinear interpolation

    IF(pointsb(1,2) .EQ. pointsb(2,2)) THEN						! case 1 (i1-i2 = 0)

      A =  (pointsb(1,3)-pointsb(2,3)) / (pointsb(1,1)-pointsb(2,1))			
      B = ((pointsb(2,3)-pointsb(3,3)) - (pointsb(2,1)-pointsb(3,1))) &
        /  (pointsb(2,2)-pointsb(3,2))	
      C =   pointsb(3,3)-pointsb(3,1)*A - pointsb(3,2)*B

    ELSE IF(pointsb(2,2) .EQ. pointsb(3,2)) THEN					! case 1 (i2-i3 = 0)

      A =  (pointsb(2,3)-pointsb(3,3)) / (pointsb(2,1)-pointsb(3,1))			
      B = ((pointsb(1,3)-pointsb(2,3)) - (pointsb(1,1)-pointsb(2,1))) &
        /  (pointsb(1,2)-pointsb(2,2))	
      C =   pointsb(3,3)-pointsb(3,1)*A - pointsb(3,2)*B

    ELSE										! case 3 ((i1-i2 != 0) and (i2-i3 != 0) - no singularities)

      A = ((pointsb(1,3)-pointsb(2,3)) / (pointsb(1,2)-pointsb(2,2))  &
        -  (pointsb(2,3)-pointsb(3,3)) / (pointsb(2,2)-pointsb(3,2))) &
        / ((pointsb(1,1)-pointsb(2,1)) / (pointsb(1,2)-pointsb(2,2))  &
        -  (pointsb(2,1)-pointsb(3,1)) / (pointsb(2,2)-pointsb(3,2))) 
      B =  (pointsb(1,3)-pointsb(2,3)) / (pointsb(1,2)-pointsb(2,2))  &
        - ((pointsb(1,1)-pointsb(2,1)) / (pointsb(1,2)-pointsb(2,2)))*A
      C =   pointsb(3,3)-pointsb(3,1)*A - pointsb(3,2)*B

  END IF

    phib 	= A*kb + B*ib + C							! scalar
    fluxZb 	= (Dm*Dmcf*A)/zcf							! flux in the z-direction
    fluxXb 	= (Dm*Dmcf*B)/zcf							! flux in the x-direction
    fluxNb 	= fluxZb*SinTheta + fluxXb*(-CosTheta)					! flux in the n-direction

  ELSE

    OPEN(1000,FILE="error.txt")
    WRITE(1000,*) "Error in MeasureFlux in Output.f90: numpointsb is not 3 or 4..."
	 WRITE(1000,*) "numpointsb=", numpointsb, "k:", k, "iteration=", filenum(n)
  	 CLOSE(1000)
  	 STOP

  END IF

  ! SECOND POINT
  ! find the coordinates of point c (2 grid spacings away from point a -- 1 away from point b)
  zc = zb + zcf*SinTheta			
  xc = xb + zcf*(-CosTheta)

  ! lattice units
  kc = (zc/zcf) + 0.5_dbl	
  ic = (xc/xcf) + 1.0_dbl

  ! find the coordinates of the points on the cell containing point c
  kLLc = FLOOR(kc)		      	! lower left lattice node
  kLRc = kLLc + 1_lng				! lower right lattice node

  ! periodic bcs
  IF(kLLc .EQ. 0) THEN
    kLLc = nz
  ELSE IF(kLRc .EQ. nz+1) THEN
    kLRc = 1_lng
  END IF

  kULc = kLLc				      	! upper left lattice node
  kURc = kLRc							! upper right lattice node

  iLLc = FLOOR(ic)
  iLRc = iLLc
  iULc = iLLc + 1_lng
  iURc = iULc

  ! find the scalar at the points on the cell containing point c
  phiLLc = phiPlane(kLLc,iLLc)	! lower left
  phiULc = phiPlane(kULc,iULc)	! upper left
  phiURc = phiPlane(kURc,iURc)	! upper right
  phiLRc = phiPlane(kLRc,iLRc)	! lower right

  ! find the number of points on the cell containing point c that are in the fluid domain (non-solid) and store the coordinates
  pointsc 		= 0.0_dbl	! initialize points
  numpointsc 	= 0_lng		! initialize numpointsc
 
  IF(nodePlane(kLLc,iLLc) .EQ. FLUID) THEN	! lower left
    numpointsc = numpointsc + 1_lng
    pointsc(numpointsc,1) = kLLc
    pointsc(numpointsc,2) = iLLc
    pointsc(numpointsc,3) = phiLLc
  END IF

  IF(nodePlane(kULc,iULc) .EQ. FLUID) THEN				! upper left
    numpointsc = numpointsc + 1_lng
    pointsc(numpointsc,1) = kULc
    pointsc(numpointsc,2) = iULc
    pointsc(numpointsc,3) = phiULc
  END IF

  IF(nodePlane(kURc,iURc) .EQ. FLUID) THEN				! upper right
    numpointsc = numpointsc + 1_lng
    pointsc(numpointsc,1) = kURc
    pointsc(numpointsc,2) = iURc
    pointsc(numpointsc,3) = phiURc
  END IF

  IF(nodePlane(kLRc,iLRc) .EQ. FLUID) THEN				! lower right
    numpointsc = numpointsc + 1_lng
    pointsc(numpointsc,1) = kLRc
    pointsc(numpointsc,2) = iLRc
    pointsc(numpointsc,3) = phiLRc
  END IF

  ! find the scalar and fluxes at point c
  IF(numpointsc .EQ. 4) THEN						! four-point bilinear interpolation

    A = phiLLc - phiLRc + phiURc - phiULc
    B = (phiURc - phiULc) - (iLLc + 1_lng)*A
    C = -(kLLc + 1_lng)*A - (phiLRc - phiURc)
    DD = phiULc - kLLc*(iLLc + 1_lng)*A - kLLc*B - (iLLc + 1_lng)*C

    phic 	= A*kc*ic + B*kc + C*ic + D				! scalar
    fluxZc 	= (Dm*Dmcf*(A*ic + B))/zcf				! flux in the z-direction
    fluxXc 	= (Dm*Dmcf*(A*kc + C))/zcf				! flux in the x-direction
    fluxNc 	= fluxZc*SinTheta + fluxXc*(-CosTheta)			! flux in the n-direction

  ELSE IF(numpointsc .EQ. 3) THEN					! three-point bilinear interpolation

    IF(pointsc(1,2) .EQ. pointsc(2,2)) THEN				! case 1 (i1-i2 = 0)

      A =  (pointsc(1,3)-pointsc(2,3)) / (pointsc(1,1)-pointsc(2,1))			
      B = ((pointsc(2,3)-pointsc(3,3)) - (pointsc(2,1)-pointsc(3,1))) &
        /  (pointsc(2,2)-pointsc(3,2))	
      C =   pointsc(3,3)-pointsc(3,1)*A - pointsc(3,2)*B

    ELSE IF(pointsc(2,2) .EQ. pointsc(3,2)) THEN			! case 1 (i2-i3 = 0)

      A =  (pointsc(2,3)-pointsc(3,3)) / (pointsc(2,1)-pointsc(3,1))			
      B = ((pointsc(1,3)-pointsc(2,3)) - (pointsc(1,1)-pointsc(2,1))) &
        /  (pointsc(1,2)-pointsc(2,2))	
      C =   pointsc(3,3)-pointsc(3,1)*A - pointsc(3,2)*B

    ELSE								! case 3 ((i1-i2 != 0) and (i2-i3 != 0) - no singularities)

      A = ((pointsc(1,3)-pointsc(2,3)) / (pointsc(1,2)-pointsc(2,2))  &
        -  (pointsc(2,3)-pointsc(3,3)) / (pointsc(2,2)-pointsc(3,2))) &
        / ((pointsc(1,1)-pointsc(2,1)) / (pointsc(1,2)-pointsc(2,2))  &
        -  (pointsc(2,1)-pointsc(3,1)) / (pointsc(2,2)-pointsc(3,2))) 
      B =  (pointsc(1,3)-pointsc(2,3)) / (pointsc(1,2)-pointsc(2,2))  &
        - ((pointsc(1,1)-pointsc(2,1)) / (pointsc(1,2)-pointsc(2,2)))*A
      C =   pointsc(3,3)-pointsc(3,1)*A - pointsc(3,2)*B

  END IF

    phic 	= A*kc + B*ic + C					! scalar
    fluxZc 	= (Dm*Dmcf*A)/zcf					! flux in the z-direction
    fluxXc 	= (Dm*Dmcf*B)/zcf					! flux in the x-direction
    fluxNc 	= fluxZc*SinTheta + fluxXc*(-CosTheta)			! flux in the n-direction

  ELSE

    OPEN(1000,FILE="error.txt")
    WRITE(1000,*) "Error in MeasureFlux in Output.f90: numpointsc is not 3 or 4..."
	 WRITE(1000,*) "numpointsc=", numpointsc, "k:", k, "iteration=", filenum(n)
  	 CLOSE(1000)
  	 STOP

  END IF

  ! COMPUTE FLUX
  wflux(k) = 2.0_dbl*(2.0_dbl*fluxNb - fluxNc)				! extrapolate the fluxes to the boundary (2x for top and bottom)

END DO

!===================================================================================================
END SUBROUTINE Compute_flux
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

  OPEN(2474,FILE='SA.dat')
  READ(2474,*)																																	! first line is variable info
  READ(2474,*)																																	! second line is zone info

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

    READ(2474,*) SAtime, SA(i)											! read in surface area from file for calculation of scalar flux/USL

  END DO

  CLOSE(2473)													! close output file (combined)
  CLOSE(2474)													! close the surface area file

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
SUBROUTINE FixMass			! enforces conservation of mass in the system 
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k			! index variables
INTEGER(lng) :: fluid_nodes		! number of fluid nodes
REAL(dbl) :: mass_actual		! mass in the system (per unit volume)
REAL(dbl) :: mass_theoretical		! mass in the system (per unit volume)
REAL(dbl) :: volume, node_volume	! total volume and volume of a sincle node (cell)
REAL(dbl) :: mass_error			! mass percent error (actual vs. theoretical)
REAL(dbl) :: mass_corrector		! mass correction factor

!----- calculate the node volume
!node_volume = xcf*ycf*zcf
node_volume = 1.0_dbl

!----- initialize the mass and node count to 0
mass_actual = 0.0_dbl
fluid_nodes = 0_lng

!----- cacluate the mass in the system based on the density and the number of fluid nodes
DO k=1,nzSub
   DO j=1,nySub
      DO i=1,nxSub
         IF (node(i,j,k) .EQ. FLUID) THEN
            mass_actual = mass_actual + (rho(i,j,k)*node_volume)
            fluid_nodes = fluid_nodes + 1_lng
         END IF 
      END DO
   END DO
END DO

!----- calcuate the theoretical amount of mass in the system
mass_theoretical = (fluid_nodes*node_volume)*denL

!----- calculate the deviation from the theoretical mass value
mass_error = mass_theoretical - mass_actual

!----- calculate the mass correction factor, mass to add to each distribution function
mass_corrector = (mass_error/fluid_nodes)/15.0_dbl

!----- correct the mass
f = f + mass_corrector

!===================================================================================================
END SUBROUTINE FixMass
!===================================================================================================








!===================================================================================================
END MODULE Output
!===================================================================================================
