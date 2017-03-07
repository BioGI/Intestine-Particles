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

ALLOCATE(filenum(0:nt))     ! maximum number of output files (should only output ~numOuts times)
filenum = 0_lng             ! initialize to 0
fileCount = 0_lng           ! initialize to 0

ALLOCATE(parfilenum(0:nt))		! maximum number of output files (should only output ~numOuts times)
ALLOCATE(numparticleSubfile(0:nt))	! maximum number of output files (should only output ~numOuts times)
parfilenum = 0_lng			! initialize to 0
parfileCount = 0_lng			! initialize to 0
numparticleSubfile = 0_lng

ALLOCATE(radius(0:nz+1,0:10000))		! 10000000 is  an arbitrarily large number of output iterations...
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

!  !----- Monitoring over saturation
!  OPEN(2119,FILE='Over_Saturation.dat',POSITION='APPEND')
!  WRITE(2119,'(A140)') 'VARIABLES = iter,  Number of oversaturated nodes, largest phi (C_max), Average of oversaturated nodes (C/Cs)'

   !----- Monitoring computational costs
   OPEN(5,FILE='Computational_Time.dat')										
   WRITE(5,'(A100)') '#VARIABLES = iter, Computational cost at each iteration, Average computational cost so far'

   !---- Volume -----
   OPEN(2460,FILE='volume.dat')
   WRITE(2460,*) 'VARIABLES = "period", "volume"'

   !----- Drug Conservation
   OPEN(2472,FILE='Drug-Conservation-'//sub//'.dat')
   WRITE(2472,'(A145)') '#VARIABLES =iter,time, Initial, Released_del_diff, Released_del_shear, Released_del_slip, Released_Total, Absorbed, Remained_in_Domain, Loss_Percent'
   
   !----- Diensity Correction to improve mass conservation
   OPEN(2473,FILE='Density_Correction.dat')
   WRITE(2473,*) 'VARIABLES: iter, Density Correction'
END IF

!----- Mass
OPEN(2458,FILE='mass-'//sub//'.dat')
WRITE(2458,'(A120)') '#VARIABLES = period, time, mass_theory, mass_actual, mass_err'

!===================================================================================================
END SUBROUTINE OpenOutputFiles
!===================================================================================================






!===================================================================================================
SUBROUTINE CloseOutputFiles	! opens output files
!===================================================================================================
IMPLICIT NONE

IF(myid .EQ. master) THEN
  CLOSE(2118) 			! Monitoring negative phi
! CLOSE(2119)                   ! Monitorin Over Saturation
  CLOSE(5)			! Status								
  CLOSE(2460)			! Volume
  CLOSE(2472)			! Scalar
  CLOSE(2473)                   ! Density Corrections
END IF
CLOSE(2458)			! Mass
!===================================================================================================
END SUBROUTINE CloseOutputFiles
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintComputationalTime(N_Clock)				! Writes Computational Costs to a File
!===================================================================================================
IMPLICIT NONE
INTEGER(lng) :: N_Clock
REAL(dbl)    :: Time(10000000)

IF (myid .EQ. master) THEN
   rate = 100_lng							! Set the rate of counting

   IF (iter .EQ. iter0-1_lng) THEN
      CALL SYSTEM_CLOCK(start,rate)					! Keep Track of Elapsed Time 
      Time(iter)= start	
   ELSE
      CALL SYSTEM_CLOCK(current,rate)					
      Time(iter)= current 
      WRITE(5,*) iter, N_Clock, (Time(iter)-Time(iter-1))/REAL(rate)
   END IF

   IF ((MOD(iter,50) .EQ. 0))  THEN
      CALL FLUSH(5)
   END IF
END IF
!===================================================================================================
END SUBROUTINE PrintComputationalTime
!===================================================================================================









!===================================================================================================
SUBROUTINE PrintRestart							      ! prints restart files
!===================================================================================================
IMPLICIT NONE

TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next
INTEGER(lng) :: i,j,k,m,Particle_Counter			
CHARACTER(7):: iter_char             			           ! iteration stored as a character

!----- Creating a file called iter0.dat with iteration number in it --------------------------------
IF (myid .EQ. master) THEN
   OPEN(550,FILE='Restart-iter.dat')
   WRITE(550,*) iter
   CLOSE(550)
END IF

!----- Creating a file called restart-iter-sub.dat with all macroscopic fields and distribution functions -----
WRITE(iter_char(1:7),'(I7.7)') iter
OPEN(500,FILE='Restart-Out-'//iter_char//'-'//sub//'.dat')			
DO k=0,nzSub+1
   DO j=0,nySub+1
      DO i=0,nxSub+1
         WRITE(500,'(I1)') node(i,j,k)
         IF (node(i,j,k) .EQ. FLUID) THEN
            WRITE(500,'(F8.6)') u(i,j,k)
            WRITE(500,'(F8.6)') v(i,j,k)
            WRITE(500,'(F8.6)') w(i,j,k)
            WRITE(500,'(F8.6)') rho(i,j,k)
            WRITE(500,'(F8.6)') phi(i,j,k)
            DO m=0,NumDistDirs
               WRITE(500,'(F10.8)') f(m,i,j,k)
            END DO
         ELSE 
            WRITE(500,'(I1)') 0
            WRITE(500,'(I1)') 0
            WRITE(500,'(I1)') 0
            WRITE(500,'(I1)') 0
            WRITE(500,'(I1)') 0
            DO m=0,NumDistDirs
               WRITE(500,'(I1)') 0
            END DO
         ENDIF
      END DO
   END DO
END DO
WRITE(500,*) Drug_Initial

IF (myid .EQ. master) THEN
   WRITE(500,*) Drug_Released_Total
ELSE
   WRITE(500,*) 0
END IF   
WRITE(500,*) Drug_Absorbed
WRITE(500,*) Drug_Remained_in_Domain
CLOSE(500)


!----- Creating a file called particle-restart-iter.dat with all the particle data in it ----------------
IF ((myid .EQ. master) .AND. (Flag_ParticleTrack) .AND. (iter .GE. iter_Start_phi)) THEN
   OPEN(156,FILE='Restart-Particles-'//iter_char//'.dat')

   !----- Counting hte number of the particles which are not compeletel dissolved yet ------------------
   Particle_Counter = 0
   current => ParListHead%next
   DO WHILE (ASSOCIATED(current))
      next => current%next
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                                           ! only write particle data when particle is not fully dissolved
         Particle_Counter = Particle_Counter + 1
      END IF
      current => next
   END DO

   write(156,*) Particle_Counter 
   current => ParListHead%next							 
   DO WHILE (ASSOCIATED(current))
      next => current%next 
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                                           ! only write particle data when particle is not fully dissolved
         WRITE(156,*)  	current%pardata%parid,	   &
                        current%pardata%xp, 	  	 &
                        current%pardata%yp,  	  	 &
                        current%pardata%zp, 	  	 &
                        current%pardata%up, 	  	 &
                        current%pardata%vp, 	     &
                        current%pardata%wp, 	  	 &
                        current%pardata%rp,        &
                        current%pardata%xpold,     &
                        current%pardata%ypold,     &
                        current%pardata%zpold,     &
                        current%pardata%upold,     & 
                        current%pardata%vpold,     &
                        current%pardata%wpold,     &
                        current%pardata%rpold,     &
                        current%pardata%par_conc,  &
                        current%pardata%gamma_cont,&
                        current%pardata%sh_conf,   &
                        current%pardata%sh_shear,  &
                        current%pardata%sh_slip,   &
                        current%pardata%S,    		 &
                        current%pardata%Sst,   		 &
                        current%pardata%Veff,      &
                        current%pardata%Nbj,		   &		
                        current%pardata%bulk_conc, &
                        current%pardata%delNBbyCV, &
                        current%pardata%cur_part,	 &
                        current%pardata%new_part  		
1001                    format (I5,24F18.10,2I5)
        END IF
     current => next  
  ENDDO
  CLOSE(156)
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

IF ((MOD(iter, Output_Intervals) .EQ. 0) 	   .OR. &
   (iter .EQ. iter0-1_lng) .OR. (iter .EQ. iter0)  .OR. &
   (iter .EQ. iter_Start_phi) .OR. (iter .EQ. nt)) THEN
   !----- scale the iteration by 1/10 such that the numbers used in the output file aren't too large
   WRITE(iter_char(1:7),'(I7.7)') iter

   !----- store the current iteration in "filenum"
   filenum(fileCount) = iter
   fileCount = fileCount + 1_lng

   !----- open the proper output file
   OPEN(60,FILE='out-'//iter_char//'-'//sub//'.dat')

   IF (iter .LT. iter_Freeze_LBM) THEN
      WRITE(60,*) 'VARIABLES = "x" "y" "z" "u(mm/s)" "v(mm/s)" "w(mm/s)" "P" "phi/Cs" "node"              &  
                       "dudx" "dvdx" "dwdx" "dudy" "dvdy" "dwdy" "dudz" "dvdz" "dwdz"                     &
                       "d2udx2" " d2vdx2" "d2wdx2" "d2udy2" "d2vdy2" "d2wdy2" "d2udz2" "d2vdz2" "d2wdz2"  &
                       "DUdt_x" "DUdt_y" "DUdt_z" "Laplacian_x" "Laplacian_y" "Laplacian_z"               &
                       "dA1dx" "dA1dy" "dA1dz" "dA2dx" "dA2dy" "dA2dz" "dA3dx" "dA3dy" "dA3dz"            &
                       "DLaplacianDt_x" "DLaplacianDt_y" "DLaplacianDt_z"'
      WRITE(60,'(A10,E15.5,A5,I4,A5,I4,A5,I4,A8)') 'ZONE T="',iter/(nt/nPers),'" I=',nxSub,' J=',nySub,' K=',nzSub,'F=POINT'
      DO k=1,nzSub
         DO j=1,nySub
            DO i=1,nxSub
               !----- convert local i,j,k, to global ii,jj,kk
               ii = ((iMin - 1_lng) + i)
               jj = ((jMin - 1_lng) + j)
               kk = ((kMin - 1_lng) + k)
               pressure= (rho(i,j,k)-denL)*dcf*pcf
               IF (node(i,j,k) .EQ. FLUID) THEN
                  WRITE(60,'(I3,2I4,3F7.2,E11.3,F9.5,I2,36E18.10)') ii, jj,kk,                                                                           &
                       1000.0_dbl*u(i,j,k)*vcf, 1000.0_dbl*v(i,j,k)*vcf, 1000.0_dbl*w(i,j,k)*vcf, pressure, phi(i,j,k)/S_intrinsic, node(i,j,k),         &
                       dudx(i,j,k),dvdx(i,j,k), dwdx(i,j,k),dudy(i,j,k),dvdy(i,j,k),dwdy(i,j,k),dudz(i,j,k),dvdz(i,j,k),dwdz(i,j,k),                     &
                       d2udx2(i,j,k), d2vdx2(i,j,k),d2wdx2(i,j,k),d2udy2(i,j,k),d2vdy2(i,j,k),d2wdy2(i,j,k),d2udz2(i,j,k),d2vdz2(i,j,k),d2wdz2(i,j,k),   &
                       DUdt_x(i,j,k),DUdt_y(i,j,k),DUdt_z(i,j,k),                                                                                        &
                       Laplacian_x(i,j,k),Laplacian_y(i,j,k),Laplacian_z(i,j,k),                                                                         &
                       dA1dx(i,j,k),dA1dy(i,j,k),dA1dz(i,j,k),dA2dx(i,j,k),dA2dy(i,j,k),dA2dz(i,j,k),dA3dx(i,j,k),dA3dy(i,j,k),dA3dz(i,j,k),             &
                       DLaplacianDt_x(i,j,k),DLaplacianDt_y(i,j,k),DLaplacianDt_z(i,j,k)
               ELSE
                  WRITE(60,'(I3,2I4,6I2,27I2)') ii,jj,kk,0,0,0,0,0,node(i,j,k),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

               END IF
            END DO
         END DO
      END DO
      CLOSE(60)
   ELSE IF (iter .GE. iter_Freeze_LBM) THEN
      WRITE(60,*) 'VARIABLES = "x" "y" "z" "phi/S_interinsic" "node"'
      WRITE(60,'(A10,E15.5,A5,I4,A5,I4,A5,I4,A8)') 'ZONE T="',iter/(nt/nPers),'" I=',nxSub,' J=',nySub,' K=',nzSub,'F=POINT'
      DO k=1,nzSub
         DO j=1,nySub
            DO i=1,nxSub
               !----- convert local i,j,k, to global ii,jj,kk
               ii = ((iMin - 1_lng) + i)
               jj = ((jMin - 1_lng) + j)
               kk = ((kMin - 1_lng) + k)
               IF (node(i,j,k) .EQ. FLUID) THEN
                  WRITE(60,'(I3,2I4,F9.5,I2)') ii,jj,kk, phi(i,j,k)/S_intrinsic, node(i,j,k)
               ELSE
                  WRITE(60,'(I3,2I4,2I2)')     ii,jj,kk,0,node(i,j,k)
               END IF
            END DO
         END DO
      END DO
      CLOSE(60)
   END IF 
 
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
SUBROUTINE PrintParticles						  ! prints all particle data
!===================================================================================================
IMPLICIT NONE
INTEGER(lng) 	:: numParticlesSub
INTEGER(lng)	:: i,j,k,ii,jj,kk,n		! index variables (local and global)
CHARACTER(7)	:: iter_char				! iteration stored as a character
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

IF (myid .EQ. master) THEN
   IF ((MOD(iter, Output_Intervals) .EQ. 0)           		& 
      .OR. (iter .EQ. iter0-1_lng) .OR. (iter .EQ. iter0) 	&
      .OR. (iter .EQ. iter_Start_phi) .OR. (iter .EQ. nt)) THEN
   
       WRITE(iter_char(1:7),'(I7.7)') iter

      !------ store the current iteration in "parfilenum"
      parfilenum(parfileCount) = iter
      numParticlesSub  = 0_lng

      !------ open the proper output file
      OPEN(160,FILE='pardat-'//iter_char//'-'//sub//'.csv')
      WRITE(160,*) '"x","y","z","u","v","w","ParID","Sh_conf","Sh_shear","Sh_slip","rp","Cb/Cs","delNBbyCV","Sst","S","C_surface","CPU"'
      current => ParListHead%next							 

      DO WHILE (ASSOCIATED(current))
         numParticlesSub = numParticlesSub + 1_lng
         next => current%next 	
         IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                                           ! only write particle data when particle is not fully dissolved
             WRITE(160,1001) current%pardata%xp          ,',',	&
                             current%pardata%yp          ,',',	&
                             current%pardata%zp          ,',',	&
                             current%pardata%up*vcf 	   ,',',	&
                             current%pardata%vp*vcf 	   ,',',	&
                             current%pardata%wp*vcf 	   ,',',	&
                             current%pardata%U_slp  	   ,',',	&
                             current%pardata%parid 	  	 ,',',	&
                             current%pardata%sh_conf 	   ,',',	&
                             current%pardata%sh_shear    ,',',	&
                             current%pardata%sh_slip 	   ,',',	&
                             current%pardata%rp          ,',',	&
                             current%pardata%bulk_conc/S_intrinsic,',', 	&
                             current%pardata%delNBbyCV   ,',', 	&
                             current%pardata%Sst 	     	 ,',',	&
                             current%pardata%S 	         ,',', 	&
                             current%pardata%par_conc    ,',', 	&
                             current%pardata%cur_part    
         END IF	
1001     format (6(F8.3,a2),I6,a2,3(F9.5,a2),3(F11.8,a2),3(F13.8,a2),I4)
         current => next
      ENDDO

      CLOSE(160)
      numparticleSubfile(parfileCount) = numParticlesSub	
      parfileCount = parfileCount + 1_lng
  ENDIF 
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
volume     = 0.0_dbl

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
IF ((MOD(iter,50) .EQ. 0))  THEN
   CALL FLUSH(2458)  
ENDIF
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
   IF ((MOD(iter,50) .EQ. 0))  THEN
      CALL FLUSH(2460)  
   ENDIF   
END IF
!===================================================================================================
END SUBROUTINE PrintVolume
!===================================================================================================








!===================================================================================================
SUBROUTINE PrintDrugConservation		! prints the total amount of scalar absorbed through the walls 
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: i,j,k,mpierr						! index variables
INTEGER(lng) :: numFluids						! number of fluid nodes in the domain
REAL(dbl)    :: phiDomain, phiIC, Drug_Initial,Drug_Negative		! current amount of scalar in the domain
REAL(dbl)    :: phiAverage						! average scalar in the domain
REAL(dbl)    :: zcf3							! node volume in physical units
REAL(dbl)    :: phiTotal_Global, phiAbsorbed_Global, phiDomain_Global
REAL(dbl)    :: Sherwood, Drug_Released_del_diff, Drug_Released_del_shear, Drug_Released_del_slip 
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
   IF ((Flag_ParticleTrack) .AND. (iter .GE. iter_Start_phi)) THEN
      Drug_Released_del_diff = 0.0_dbl
      Drug_Released_del_shear= 0.0_dbl
      Drug_Released_del_slip = 0.0_dbl

      current => ParListHead%next
      DO WHILE (ASSOCIATED(current))
         next => current%next
         IF (current%pardata%rp .GT. Min_R_Acceptable) THEN
            Sherwood               = 1.0_dbl + current%pardata%sh_shear + current%pardata%sh_slip
            Drug_Released_del_diff = Drug_Released_del_diff  + current%pardata%delNBbyCV*zcf3*(1.0_dbl                 /Sherwood)                                                   
            Drug_Released_del_shear= Drug_Released_del_shear + current%pardata%delNBbyCV*zcf3*(current%pardata%sh_shear/Sherwood)                                                  
            Drug_Released_del_slip = Drug_Released_del_slip  + current%pardata%delNBbyCV*zcf3*(current%pardata%sh_slip /Sherwood) 
         END IF
         current => next
      ENDDO
      Drug_Released_Total = Drug_Released_Total + Drug_Released_del_diff + Drug_Released_del_shear + Drug_Released_del_slip  
   END IF
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
CALL MPI_ALLREDUCE(phiTotal   , phiTotal_Global   ,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(phiAbsorbed, phiAbsorbed_Global,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(phiDomain  , phiDomain_Global  ,1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

Drug_Initial		= phiTotal_Global           * zcf3
Drug_Absorbed 		= phiAbsorbed_Global 	    * zcf3 + Drug_Absorbed_restart
Drug_Remained_in_Domain = phiDomain_Global   	    * zcf3
Drug_Negative		= Negative_phi_Total_Global * zcf3
Drug_Loss 		= (Drug_Released_Total+ Drug_Initial               ) - (Drug_Absorbed + Drug_Remained_in_Domain)  
Drug_Loss_Modified 	= (Drug_Released_Total+ Drug_Initial- Drug_Negative) - (Drug_Absorbed + Drug_Remained_in_Domain)

IF (Drug_Released_Total .LT. 1e-10) THEN
    Drug_Released_Total = 0.000000000_dbl
END IF
IF (Drug_Absorbed .lt. 1.0e-10) THEN
    Drug_Absorbed = 0.00000000_dbl
ENDIF
IF (Drug_Released_del_diff .lt. 1.0e-10) THEN
    Drug_Released_del_diff = 0.0000000000_dbl
ENDIF
IF (Drug_Released_del_shear .lt. 1.0e-10) THEN
    Drug_Released_del_shear = 0.0000000000_dbl
ENDIF
IF (Drug_Released_del_slip .lt. 1.0e-10) THEN
    Drug_Released_del_slip = 0.00000000_dbl
ENDIF

IF ((Drug_Released_Total+Drug_Initial).GT.1e-10) THEN
    Drug_Loss_Percent = (Drug_Loss / (Drug_Released_Total+Drug_Initial)) * 100.0_lng
    Drug_Loss_Modified_Percent = (Drug_Loss_Modified / (Drug_Released_Total+Drug_Initial)) * 100.0_lng  
ELSE
    Drug_Loss_Percent = 0.0000_dbl
    Drug_Loss_Modified_Percent = 0.0000_dbl
END IF


IF (myid .EQ. master) THEN
   WRITE(2472,'(I8, F13.4, 7E18.10, F11.6)') iter, iter*tcf, Drug_Initial, Drug_Released_del_diff, Drug_Released_del_shear, Drug_Released_del_slip, Drug_Released_Total, Drug_Absorbed, Drug_Remained_in_Domain, Drug_Loss_Percent  
   IF ((MOD(iter,50) .EQ. 0))  THEN
      CALL FLUSH(2472)
   ENDIF   
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
   WRITE(11,*) 'iter_Start_phi=',iter_Start_phi		! iteration at which to start the particle tracking & scalar calculations
   WRITE(11,*)
   WRITE(11,*) 'nPers=',nPers			! total number of periods to run
   WRITE(11,*) 'numOuts=',numOuts		! number of output files (roughly)
   WRITE(11,*) 'Flag_Restart=',Flag_Restart		! use restart file? (0 if no, 1 if yes)
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
SUBROUTINE MergeOutput          	  ! combines the subdomain output files into one output file
!===================================================================================================
IMPLICIT NONE

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
   IF ((MOD(iter,50) .EQ. 0))  THEN
      CALL FLUSH(5)
   ENDIF    

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
END MODULE Output
!===================================================================================================
