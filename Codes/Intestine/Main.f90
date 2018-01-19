!==================================================================================================
PROGRAM LBM3D	! 3D Parallelized LBM Simulation
                ! Gino Banco (2008-2010) - original LBM method parallelization for intestine
                ! Balaji Jayaraman (2014-2015) - Improved LBM method , particle tracking in parallel and drug release model
!==================================================================================================
USE SetPrecision 
USE Setup
USE Parallel    
USE LBM      
USE Geometry
USE PassiveScalar
USE IC
USE BClbm	
USE Output
USE ParticleTracking
USE ParticleDrug

IMPLICIT NONE

INTEGER(lng) :: mpierr						! MPI standard error variable
INTEGER(lng) :: i,j,k,ii,jj					! lattice indices
REAL(dbl)    :: phidomf,phidomfs				! current amount of scalar in the domain
INTEGER, allocatable :: seed(:)
INTEGER      :: seed_size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MPI Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALL MPI_INIT(mpierr)						! initialize parallelization [Intrinsic]
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,mpierr)		! get the size of the parallel "world" (number of processing units) [Intrinsic]
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,mpierr)			! assign each prossessing unit a number (myid) for identification [Intrinsic]

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALL RANDOM_SEED(size=seed_size)
ALLOCATE(seed(seed_size))
seed=10972
CALL RANDOM_SEED(put=seed)
DEALLOCATE(seed)
	
CALL Global_Setup				    ! set up the simulation {MODULE: Setup]
CALL MPI_Setup						  ! set up MPI component of the simulation [MODULE: Parallel]
CALL LBM_Setup							! set up LBM simulation [MODULE: LBM]
CALL Scalar_Setup						! set up the passive scalar component of the simluation [MODULE: Scalar]
CALL Output_Setup						! set up the output [MODULE: Output]
CAll AdvanceGeometry 
CALL ICs							      ! set initial conditions [MODULE: ICBC]
CALL OpenOutputFiles				! opens output files for writing [MODULE: Output.f90]
CALL PrintParams						! print simulation info [MODULE: Output]
CALL PrintFields						! output the velocity, density, and scalar fields [MODULE: Output]
!CALL PrintComputationalTime(0) 						! Start simulation timer, print status [MODULE: Output]

IF (Flag_ParticleTrack) THEN 					! If particle tracking is 'on' then do the following
   CALL IniParticles
   CALL Particle_Setup
ENDIF

IF (Pw .GT. 0.0) THEN
   CALL Permeability_Nodes
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)				! synchronize all processes before starting simulation [Intrinsic]

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DO iter = iter0-0_lng,nt
   IF (iter .LT. iter_Freeze_LBM) THEN
      CALL AdvanceGeometry						! advance the geometry to the next time step [MODULE: Geometry]
      CALL Collision							! collision step [MODULE: Algorithm]
      CALL MPI_Transfer						        ! transfer the data (distribution functions, density, scalar) [MODULE: Parallel]
      IF (domaintype .EQ. 0) THEN  					! only needed when planes of symmetry exist
         CALL SymmetryBC						! enforce symmetry boundary condition at the planes of symmetry [MODULE: ICBC]
      ENDIF
      CALL Stream							! perform the streaming operation (with Lallemand 2nd order BB) [MODULE: Algorithm]
      CALL Macro							! calcuate the macroscopic quantities [MODULE: Algorithm]
      IF ((Flag_Convection_Effects) .AND. ((MOD(iter,Output_Intervals).EQ.0).OR.(MOD(iter,Restart_Intervals).EQ.0))) THEN
         IF ((NumSubsX*NumSubsY*NumSubsZ).EQ.1)THEN
            CALL Compute_vel_derivatives
         ENDIF   
      ENDIF   
   ELSEIF ((iter .GE. iter_Freeze_LBM) .AND. (iter .EQ. iter0)) THEN 
      CALL AdvanceGeometry						! advance the geometry to the next time step [MODULE: Geometry]
      CALL Collision							! collision step [MODULE: Algorithm]
      CALL MPI_Transfer					        	! transfer the data (distribution functions, density, scalar) [MODULE: Parallel]
      IF (domaintype .EQ. 0) THEN  					! only needed when planes of symmetry exist
          CALL SymmetryBC						! enforce symmetry boundary condition at the planes of symmetry [MODULE: ICBC]
      ENDIF
   ENDIF
   CALL MPI_Transfer 
   CALL PrintComputationalTime(1) 						! Start simulation timer, print status [MODULE: Output]
   IF ((Flag_ParticleTrack) .AND. (iter .GE. iter_Start_phi)) THEN 	! If particle tracking is 'on' then do the following
      CALL Particle_Track
   ENDIF
   IF (iter .GE. iter_Start_phi) THEN
      CALL Scalar							! calcuate the evolution of scalar in the domain [MODULE: Algorithm]
   END IF
   CALL PrintComputationalTime(10)						! Start simulation timer, print status [MODULE: Output]
   !----- Outputs------------------------------------------------------
   CALL PrintFields					 	! output the velocity, pressure and scalar fields 
   IF ((Flag_ParticleTrack) .AND. (iter .GE. iter_Start_phi)) THEN 	! If particle tracking is 'on' then do the following
      CALL PrintParticles						! output the particle velocity, radius, position and con. [MODULE: Output]
   ENDIF
   CALL PrintDrugConservation						! print the total absorbed/entering/leaving scalar as a function of time [MODULE: Output]
   !CALL PrintMass							! print the total mass in the system (TEST)
   CALL PrintVolume						! print the volume in the system (TEST)
   CALL PrintComputationalTime(11)						! print current status [MODULE: Output]
   IF (MOD(iter,Restart_Intervals) .EQ. 0) THEN
      CALL PrintRestart
   END IF
   CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)				! synchronize all processing units before next loop [Intrinsic]
  
END DO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End Simulation Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CALL PrintRestart							! print a final set of restart files to continue if desired [MODULE: Output]
CALL DEAllocateArrays						! clean up the memory [MODULE: Setup]
CALL CloseOutputFiles						! closes output files [MODULE: Output.f90]
!CALL MergeOutput							! combine the subdomain output into an output file for the entire computational domain [MODULE: Output]
CALL MPI_TYPE_FREE(mpipartransfertype,mpierr)
CALL MPI_FINALIZE(mpierr)						! end the MPI simulation [Intrinsic]

!================================================ 
END PROGRAM LBM3D
!================================================
