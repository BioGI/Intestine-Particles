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
	USE ICBC	
	USE Output

	IMPLICIT NONE

	INTEGER(lng) :: mpierr						! MPI standard error variable
	INTEGER(lng) :: i,j,k,ii,jj					! lattice indices
	REAL(dbl) :: phidomf,phidomfs					! current amount of scalar in the domain
	INTEGER, allocatable :: seed(:)
	INTEGER :: seed_size

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
	
	CALL Global_Setup						! set up the simulation {MODULE: Setup]
	CALL MPI_Setup							! set up MPI component of the simulation [MODULE: Parallel]
        CALL LBM_Setup							! set up LBM simulation [MODULE: LBM]
	CALL Geometry_Setup						! set up the geometry of the physcial simulation [MODULE: Geometry]
	CALL Scalar_Setup						! set up the passive scalar component of the simluation [MODULE: Scalar]
	CALL Output_Setup						! set up the output [MODULE: Output]
	CALL ICs							! set initial conditions [MODULE: ICBC]
	CALL OpenOutputFiles						! opens output files for writing [MODULE: Output.f90]
	CALL PrintParams						! print simulation info [MODULE: Output]
	CALL PrintFields						! output the velocity, density, and scalar fields [MODULE: Output]
	CALL PrintStatus 						! Start simulation timer, print status [MODULE: Output]

	IF(ParticleTrack.EQ.ParticleOn) THEN 				! If particle tracking is 'on' then do the following
           CALL IniParticles
           CALL Particle_Setup
	ENDIF

	IF(restart) THEN						! calculate the villous locations/angles at iter0-1 [MODULE: Geometry]
           CALL AdvanceGeometry
	END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)				! synchronize all processes before starting simulation [Intrinsic]

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO iter = iter0-0_lng,nt
	CALL AdvanceGeometry						! advance the geometry to the next time step [MODULE: Geometry]
	CALL Collision							! collision step [MODULE: Algorithm]
	CALL MPI_Transfer						! transfer the data (distribution functions, density, scalar) [MODULE: Parallel]

	IF(domaintype .EQ. 0) THEN  					! only needed when planes of symmetry exist
          CALL SymmetryBC						! enforce symmetry boundary condition at the planes of symmetry [MODULE: ICBC]
	ENDIF

	IF(ParticleTrack.EQ.ParticleOn .AND. iter .GE. phiStart) THEN 	! If particle tracking is 'on' then do the following
	   CALL Calc_Global_Bulk_Scalar_Conc				! Estimate bluk	scalar concentration in each partition
	   CALL Collect_Distribute_Global_Bulk_Scalar_Conc		! Collect Cb_Global from different processors, average it and distribute it to all the processors.  
	   CALL Particle_Track
	   CALL Particle_MPI_Transfer
	ENDIF

	CALL Stream							! perform the streaming operation (with Lallemand 2nd order BB) [MODULE: Algorithm]
	CALL Macro							! calcuate the macroscopic quantities [MODULE: Algorithm]

	IF(iter .GE. phiStart) THEN
           CALL Scalar							! calcuate the evolution of scalar in the domain [MODULE: Algorithm]
	END IF

	IF (myid .EQ. 0) THEN
	   open(70,file='t-history.dat',position='append')
           write(70,*) iter, u(41,41,5)*vcf,v(41,41,5)*vcf,w(41,41,5)*vcf,sqrt(u(60,61,10)**2+v(60,61,10)**2)*vcf,u(70,71,10)*vcf,v(70,71,10)*vcf,w(70,71,10)*vcf,sqrt(u(70,71,10)**2+v(70,71,10)**2)*vcf
	   close(70)
	ENDIF

        CALL PrintFields						! output the velocity, density, and scalar fields [MODULE: Output]

        IF(ParticleTrack.EQ.ParticleOn .AND. iter .GE. phiStart) THEN 	! If particle tracking is 'on' then do the following
     	   CALL PrintParticles						! output the particle velocity, radius, position and con. [MODULE: Output]
        ENDIF

    	CALL PrintScalar						! print the total absorbed/entering/leaving scalar as a function of time [MODULE: Output]
     	CALL PrintMass							! print the total mass in the system (TEST)
     	CALL PrintVolume						! print the volume in the system (TEST)
	CALL PrintStatus 						! print current status [MODULE: Output]
	CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)				! synchronize all processing units before next loop [Intrinsic]

      END DO

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End Simulation Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL PrintFinalRestart						! print a final set of restart files to continue if desired [MODULE: Output]
      CALL DEAllocateArrays						! clean up the memory [MODULE: Setup]
      CALL CloseOutputFiles						! closes output files [MODULE: Output.f90]
      CALL MergeOutput							! combine the subdomain output into an output file for the entire computational domain [MODULE: Output]
      CALL MPI_TYPE_FREE(mpipartransfertype,mpierr)
      CALL MPI_FINALIZE(mpierr)						! end the MPI simulation [Intrinsic]

!================================================ 
END PROGRAM LBM3D
!================================================
