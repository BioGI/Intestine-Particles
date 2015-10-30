!==================================================================================================
MODULE Setup	! Defines global variables, reads input from file, allocates arrays
!==================================================================================================
USE SetPrecision

IMPLICIT NONE

!**************************** Global Simulation Quanitities ***************************************

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LBM Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

INTEGER(lng), PARAMETER :: NumDistDirs	= 14_lng				! number of distribution function directions minus one (ex. D3Q15 -> 14)
REAL(dbl),		ALLOCATABLE :: f(:,:,:,:)				! distribution function
REAL(dbl), 		ALLOCATABLE :: fplus(:,:,:,:)				! post-collision distribution function
REAL(dbl), 		ALLOCATABLE :: u(:,:,:),v(:,:,:),w(:,:,:)		! x,y, and z components of the fluid velocity vector
REAL(dbl), 		ALLOCATABLE :: dudx(:,:,:),dudy(:,:,:),dudz(:,:,:)	! x,y, and z components of the fluid velocity vector
REAL(dbl), 		ALLOCATABLE :: rho(:,:,:)				! density
INTEGER(lng), 	ALLOCATABLE :: node(:,:,:)    					! node flags (FLUID/SOLID)
REAL(dbl), 		ALLOCATABLE :: ex(:),ey(:),ez(:)			! LBM discretized velocity vectors
INTEGER(lng), 	ALLOCATABLE :: bb(:), sym(:,:)					! bounceback and symmetry directions
REAL(dbl), 		ALLOCATABLE :: wt(:)    				! weighting coefficients for the equilibrium distribution functions
REAL(dbl)		:: den, denL						! density (physical, lattice units)
REAL(dbl)		:: nu, nuL						! kinematic viscosity (physical, lattice units)
REAL(dbl) 		:: cs							! speed of sound on the lattice
REAL(dbl)		:: tau							! relaxation parameters of coarse and fine blocks
REAL(dbl)		:: oneOVERtau						! reciprical of tau
INTEGER(lng)	:: nx,ny,nz							! number of global nodes in the x, y, and z directions respectively
INTEGER(lng)	:: nxSub,nySub,nzSub						! number of local nodes in the each direction
INTEGER(lng)	:: iter0,iter,nt						! initial time step, timestep index, total number of timesteps
INTEGER(lng)	:: domaintype							! a flag to denote domain type - 0 for 1/4th cylinder and 1 for full cylinder
INTEGER(lng), PARAMETER :: FLUID		= 0_lng				! fluid inbetween
INTEGER(lng), PARAMETER :: SOLID		= 1_lng				! solid interior (moving)
INTEGER(lng), PARAMETER :: SOLID2		= 2_lng				! solid exterior (stationary)
INTEGER(lng), PARAMETER :: qitermax = 15_lng 					! max number of q iterations
LOGICAL :: restart								! Restart Flag

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Scalar Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dbl), ALLOCATABLE :: phi(:,:,:)						! passive scalar
REAL(dbl), ALLOCATABLE :: delphi_particle(:,:,:)				! passive scalar contribution from particles
REAL(dbl), ALLOCATABLE :: tausgs_particle_x(:,:,:)				! passive scalar contribution from particles
REAL(dbl), ALLOCATABLE :: tausgs_particle_y(:,:,:)				! passive scalar contribution from particles
REAL(dbl), ALLOCATABLE :: tausgs_particle_z(:,:,:)				! passive scalar contribution from particles
REAL(dbl), ALLOCATABLE :: phiTemp(:,:,:)					! temporary storage of passive scalar
REAL(dbl) :: Sc 								! Schmidt number
REAL(dbl) :: Dm,Dmcf								! binary molecular diffusivity (passive scalar in fluid), diffusivity conversion factor
REAL(dbl) :: Delta								! scalar parameter
REAL(dbl) :: phiIC, phiWall							! values of scalar: initial, wall, contribution from boundary
REAL(dbl) :: phiAbsorbed							! total amount of scalar absorbed up to current time
REAL(dbl) :: phiAbsorbedS							! total amount of scalar absorbed up to current time - through the macroscopic surface
REAL(dbl) :: phiAbsorbedV							! total amount of scalar absorbed up to current time - through the villi
REAL(dbl) :: phiInOut								! total amount of scalar leaving/entering the domain
REAL(dbl) :: phiTotal								! total intial amount of scalar in the domain
REAL(dbl) :: sigma								! standard deviation for scalar distributions
REAL(dbl) :: phiPer								! period at which to start the scalar
INTEGER(lng) :: phiStart							! iteration to start scalar calculation scalar
INTEGER(lng) :: sclrIC								! initial condition to use (BLOB, LINE, or INLET)
REAL(dbl), PARAMETER :: ee = 2.71828182846					! e^1
INTEGER(lng), PARAMETER :: BLOB=1						! scalar initial condition: circular gaussian distribution of scalar at the center of the domain
INTEGER(lng), PARAMETER :: LINE=2						! scalar initial condition: gaussian distribution of scalar in the x,y-directions along the centerline
INTEGER(lng), PARAMETER :: INLET=3						! scalar initial condition: gaussian distribution of scalar in the z-direction along the inlet
INTEGER(lng), PARAMETER :: UNIFORM=4						! scalar initial condition: uniform scalar in the entire domain (phi=phiIC)
REAL(dbl) :: phiInNodes,phiOutNodes						! total amount of scalar leaving/entering the domain

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parallel (MPI) Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Number of communication directions (3D LBM)
INTEGER(lng), PARAMETER :: NumCommDirs		= 26_lng			! number of communication directions (6 faces + 12 sides + 8 corners = 26)
INTEGER(lng), PARAMETER :: MaxDistFns		= 5_lng				! maximum number of transfered distribution functions (actual number: 5 for faces, 2 for sides, 1 for corners)
INTEGER(lng), PARAMETER :: NumFs_face		= 5_lng				! number of distribution functions transferred between neighboring faces
INTEGER(lng), PARAMETER :: NumFs_side		= 2_lng				! number of distribution functions transferred between neighboring side
INTEGER(lng), PARAMETER :: NumFs_corner	= 1_lng					! number of distribution functions transferred between neighboring corner

! MPI Arrays (arranged by descending size for storage efficiency)
INTEGER(lng), ALLOCATABLE :: f_Comps(:,:)					! specifies the components of the distribution functions to transfer in each MPI communication direction
INTEGER(lng), ALLOCATABLE :: Corner_SendIndex(:,:)				! i, j, and k indices for each corner
INTEGER(lng), ALLOCATABLE :: Corner_RecvIndex(:,:)				! i, j, and k indices for each corner (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: Z_SendIndex(:,:)					! i and j indices for each Z side 
INTEGER(lng), ALLOCATABLE :: Z_RecvIndex(:,:)					! i and j indices for each Z side (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: X_SendIndex(:,:)					! j and k indices for each X side 
INTEGER(lng), ALLOCATABLE :: X_RecvIndex(:,:)					! j and k indices for each X side (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: Y_SendIndex(:,:)					! i and k indices for each Y side 
INTEGER(lng), ALLOCATABLE :: Y_RecvIndex(:,:)					! i and k indices for each Y side (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: YZ_SendIndex(:)					! i index for each YZ face 
INTEGER(lng), ALLOCATABLE :: YZ_RecvIndex(:)					! i index for each YZ face (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: ZX_SendIndex(:)					! j index for each ZX face 
INTEGER(lng), ALLOCATABLE :: ZX_RecvIndex(:)					! j index for each ZX face (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: XY_SendIndex(:)					! k index for each XY face 
INTEGER(lng), ALLOCATABLE :: XY_RecvIndex(:)					! k index for each XY face (phantom node for recieving data)
INTEGER(lng), ALLOCATABLE :: SubID(:)						! id number of neighboring subdomains (same as rank of processing unit working on domain)
INTEGER(lng), ALLOCATABLE :: OppCommDir(:) 					! opposite MPI communication directions (like bounceback) 
INTEGER(lng), ALLOCATABLE :: CommDataStart_f(:)					! array of starting indices in the send arrays for the distribution functions from each communication direction 
INTEGER(lng), ALLOCATABLE :: CommDataStart_rho(:)				! array of starting indices in the send arrays for the density from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_phi(:)				! array of starting indices in the send arrays for the scalar from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_u(:)					! array of starting indices in the send arrays for the scalar from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_v(:)					! array of starting indices in the send arrays for the scalar from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_w(:)					! array of starting indices in the send arrays for the scalar from each communication direction
INTEGER(lng), ALLOCATABLE :: CommDataStart_node(:)				! array of starting indices in the send arrays for the node index from each communication direction
INTEGER(lng), ALLOCATABLE :: fSize(:)						! array of the number of elements sent for each communication direction (distribution functions)
INTEGER(lng), ALLOCATABLE :: dsSize(:)						! array of the number of elements sent for each communication direction (density and scalar)
INTEGER(lng), ALLOCATABLE :: uvwSize(:)						! array of the number of elements sent for each communication direction (density and scalar)
INTEGER(lng), ALLOCATABLE :: nodeSize(:)					! array of the number of elements sent for each communication direction (density and scalar)
INTEGER(lng), ALLOCATABLE :: msgSize(:)						! array of the number of elements sent for each communication direction (density and scalar)
INTEGER(lng), ALLOCATABLE :: req(:)						! array of MPI send/receive requests 
INTEGER(lng), ALLOCATABLE :: waitStat(:,:)					! array of MPI_WAITALL status objects

REAL(dbl), ALLOCATABLE :: msgSend(:)						! array of ALL of the sent information (total)
REAL(dbl), ALLOCATABLE :: msgRecv(:)						! array of ALL of the received information (total)

! MPI Variables
INTEGER(lng), PARAMETER :: master = 0_lng					! rank (id) of the master processing unit
INTEGER(lng) 	:: numprocs, myid, mySub 					! number of processing units, rank of current processing unit, subdomain of current processing unit

REAL(dbl) 	:: CommTime_f0, CommTime_fEnd, CommTime_f			! communication time - distribution functions: start time, end time, current time
REAL(dbl) 	:: CommTime_ds0, CommTime_dsEnd, CommTime_ds			! communication time - distribution functions: start time, end time, current time

! Number of Subdomains in each direction
INTEGER(lng) :: NumSubsX							! number of subdomains in the X direction
INTEGER(lng) :: NumSubsY							! number of subdomains in the Y direction
INTEGER(lng) :: NumSubsZ							! number of subdomains in the Z direction
INTEGER(lng) :: NumSubsTotal							! total number of subdomains

! Starting/Ending indices for each subdomain
INTEGER(lng) :: iMin								! starting local i index
INTEGER(lng) :: iMax								! ending local i index
INTEGER(lng) :: jMin								! starting local j index
INTEGER(lng) :: jMax								! ending local j index
INTEGER(lng) :: kMin								! starting local k index
INTEGER(lng) :: kMax								! ending local k index

! extension for output files
CHARACTER(5) :: sub	! rank + 1 for output file extensions

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Geometry Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REAL(dbl), ALLOCATABLE		:: villiLoc(:,:)				! location of each villous
REAL(dbl), ALLOCATABLE 		:: x(:),y(:),z(:)				! physical coordinate arrays
REAL(dbl), ALLOCATABLE		:: xx(:),yy(:),zz(:)				! x,y,z arrays (global)
REAL(dbl), ALLOCATABLE 		:: ub(:),vb(:),wb(:)				! x,y, and z components of the solid boundary velocity vector
!REAL(dbl), ALLOCATABLE 		:: rDom0(:),rDom(:),r(:)		! initial, and current radius at each z-location (global), radius at each location (local)
REAL(dbl), ALLOCATABLE 		:: rDom0In(:),rDomIn(:),rIn(:)			! initial, and current radius at each z-location (global), radius at each location (local)
REAL(dbl), ALLOCATABLE 		:: rDom0Out(:),rDomOut(:),rOut(:)		! initial, and current radius at each z-location (global), radius at each location (local)
!REAL(dbl), ALLOCATABLE		:: velDom(:),vel(:)				! global and local wall velocities 
REAL(dbl), ALLOCATABLE		:: velDomIn(:),velIn(:)				! global and local wall velocities 
REAL(dbl), ALLOCATABLE		:: velDomOut(:),velOut(:)			! global and local wall velocities 
REAL(dbl), ALLOCATABLE		:: rnd(:)					! array of random numbers for random villi phase angles
INTEGER(lng), ALLOCATABLE	:: villiGroup(:)				! array of which groups the villi are in 
REAL(dbl)		:: xcf, ycf, zcf					! length conversion factors
REAL(dbl)		:: dcf, vcf, pcf					! density, velocity, pressure conversion factors
REAL(dbl)		:: tcf							! time conversion factor
REAL(dbl)		:: nPers						! number of time periods simulated
REAL(dbl)		:: Lv							! length of the villi
REAL(dbl)		:: Rv							! radius of the villi
REAL(dbl)		:: villiAngle						! maximum angle of active villous travel
INTEGER(lng)	:: iLv								! length of the villi in lattice units
INTEGER(lng)	:: Ci,Cj,Ck							! center node location (global)

INTEGER(lng)	:: randORord							! flag to determine if the villous motion is random or ordered
INTEGER(lng), PARAMETER :: RANDOM=1						! random flag for randORord
INTEGER(lng), PARAMETER :: ORDERED=2						! ordered flag for randORord

REAL(dbl), PARAMETER :: PI = 3.1415926535897932384626433832			! Pi
REAL(dbl) :: D, L								! diameter, length of the intestinal segment
REAL(dbl) :: a1, a2								! half height of the passages
REAL(dbl) :: eps1, eps2								! occlusional distances
REAL(dbl) :: amp1, amp2								! amplitude of the waves
REAL(dbl) :: epsOVERa1, epsOVERa2						! occlusional distance to half-height ratios
REAL(dbl) :: aOVERlam1,	aOVERlam2						! half height to wavelength ratios
REAL(dbl) :: lambda1, lambda2							! wavelengths
REAL(dbl) :: kw1								! wave number (peristalsis)
REAL(dbl) :: s1, s2	            						! mode velocities
REAL(dbl) :: Ts, Tp								! segmental period, peristaltic period
REAL(dbl) :: wc1, wc2								! weighting coefficients for the different modes
REAL(dbl) :: Re1, Re2								! weighting coefficients for the different modes
REAL(dbl) :: shift2								! amplitude of the segmental contraction
REAL(dbl) :: Tmix								! calculated period (mix)
REAL(dbl) :: period								! period of current simulation
REAL(dbl) :: freqRatioT								! villous frequency to macroscopic contraction frequency (theta, azimuthal)
REAL(dbl) :: freqRatioZ								! villous frequency to macroscopic contraction frequency (z, axial)
REAL(dbl) :: vFreqT								! villous frequency in the theta direction (azimuthal)
REAL(dbl) :: vFreqZ								! villous frequency in the z direction (axial)
REAL(dbl) :: activeVflagT							! flag to turn on or off active villous motion in the theta direction (azimuthal)
REAL(dbl) :: activeVflagZ							! flag to turn on or off active villous motion in the z direction (axial)
INTEGER(lng) :: nlambda2							! wavelength - number of nodes (segmental)
INTEGER(lng) :: numw1, numw2							! number of waves
INTEGER(lng) :: nzSL, nzSR							! left and right indicies of scalar domain
INTEGER(lng) :: segment								! length of the segments in the portions of the segmental contractions
INTEGER(lng) :: seg1L, seg1R							! left/right point of slope segement 1
INTEGER(lng) :: seg2L, seg2R							! left/right point of slope segement 2
INTEGER(lng) :: numVilliZ, numVilliTheta					! number of villi rows in the axial direction, number of villi per row
INTEGER(lng) :: numVilli, numVilliActual					! number of total villi, actual number of total villi
INTEGER(lng) :: numVilliGroups							! number of groups of villi

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Output Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

REAL(dbl), ALLOCATABLE		:: radius(:,:)					! radius stored during output iterations
INTEGER(lng), ALLOCATABLE 	:: filenum(:)					! array of output file numbers
INTEGER(lng)    :: numOuts							! number of output files
INTEGER(lng)	:: fileCount							! current output file number (out of total number of output files)
INTEGER(lng)	:: outFlag							! specifies whether to output in readable format (1), binaries (2), or both (3)
INTEGER(lng)    :: radcount							! counts the number of output iterations for storing the radius

! System Clock Variables (for PrintStatus)
INTEGER(lng)	:: start, current, final, rate					! timing varibles for SYSTEM_CLOCK()
REAL(dbl)	:: tStart,tEnd,tTotal,tRecv,tSum,tAvg				! timing variables for parallel scalability [MPI_WTIME()]

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Particle Tracking Variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

INTEGER(lng)	:: ParticleTrack						! a flag to denote if particle track is on (1) or off (0) 
INTEGER(lng), PARAMETER :: ParticleOn=1						! flag to signify Particle Tracking is on
INTEGER(lng), PARAMETER :: ParticleOff=0					! flag for signify if particle tracking is off
INTEGER(lng)    :: np								! number of particles

!REAL(dbl), PARAMETER :: molarvol=2.65e-10,diffm=6.7e-10			! drug properties
!REAL(dbl), PARAMETER :: molarvol=92.73_dbl,diffm=4.235e-6_dbl,R0=0.0026_dbl,Cs_mol=3.14854e-6!1.2e-6		! drug properties
REAL(dbl), PARAMETER :: molarvol=92.73_dbl,diffm=8.47e-7_dbl,R0=0.0026_dbl,Cs_mol=3.14854e-6!1.2e-6		! drug properties
REAL(dbl):: Cb_global		! Global bulk scalar Concentration
INTEGER(lng):: Cb_numFluids	! Number of fluid nodes in the process for Global bulk scalar Concentration
INTEGER(lng):: num_particles	! Total number of particles in domain

!REAL(dbl), ALLOCATABLE		:: xp(:),yp(:),zp(:)				! particle physical location coordinates
!INTEGER(lng), ALLOCATABLE 	:: ipar(:),jpar(:),kpar(:)			! particle computational nodal coordinates
!REAL(dbl), ALLOCATABLE		:: up(:),vp(:),wp(:)				! particle velocities
!REAL(dbl), ALLOCATABLE		:: rp(:),delNBbyCV(:),par_conc(:),bulk_conc(:),rpold(:)		! particle radius and drug release rate
!REAL(dbl), ALLOCATABLE		:: sh(:),gamma_cont(:)					! particle sherwood number
! Particle tracking parallel variables

INTEGER(lng), ALLOCATABLE :: iMaxDomain(:),iMinDomain(:) ! List of starting/enning i indices for each subdomain
INTEGER(lng), ALLOCATABLE :: jMaxDomain(:),jMinDomain(:) ! List of starting/enning j indices for each subdomain
INTEGER(lng), ALLOCATABLE :: kMaxDomain(:),kMinDomain(:) ! List of starting/enning k indices for each subdomain
REAL(dbl), ALLOCATABLE  :: partransfersend(:,:),partransferrecv(:,:)
INTEGER(lng),ALLOCATABLE :: parreqid(:),parwtstat(:,:)		! number of send/recv requests
INTEGER(lng),ALLOCATABLE :: probestat(:)	! MPI status object
INTEGER(lng),ALLOCATABLE :: numpartransfer(:)	! Particles to be transferred in each direction
INTEGER(lng) :: NumCommDirsPar = 26_lng
INTEGER(lng) :: NumParVar = 16_lng

!TYPE ParRecord
!	TYPE(ParRecord), POINTER :: prev => NULL()! pointer to prev record
!	TYPE(ParRecord), POINTER :: next => NULL()	! pointer to next record
!	INTEGER(lng)	:: parid ! particle id in the overall list - a tag that can be used to track the particle
!	INTEGER(lng)	:: cur_part	! current sub-domain id / partition number
!	INTEGER(lng)	:: new_part	! current sub-domain id / partition number
!	REAL(dbl)	:: xp	! particle x-position
!	REAL(dbl)	:: yp ! particle y-position
!	REAL(dbl)	:: zp	! particle z-position
!	REAL(dbl)	:: up	! particle u-velocity
!	REAL(dbl)	:: vp   ! particle v-velocity
!	REAL(dbl)	:: wp	! particle w-velocity
!	REAL(dbl)	:: rp	! particle radius
!	REAL(dbl)	:: delNBbyCV ! particle drug release concentration 
!	REAL(dbl)	:: par_conc ! particle concentration
!	REAL(dbl)	:: bulk_conc ! bulk concentration at particle location
!	REAL(dbl)	:: xpold	! particle x-position
!	REAL(dbl)	:: ypold 	! particle y-position
!	REAL(dbl)	:: zpold	! particle z-position
!	REAL(dbl)	:: upold	! particle u-velocity
!	REAL(dbl)	:: vpold   	! particle v-velocity
!	REAL(dbl)	:: wpold	! particle w-velocity
!	REAL(dbl)	:: rpold	! old particle radius
!	REAL(dbl)	:: sh 	! Sherwood number
!	REAL(dbl)	:: gamma_cont	! gamma - container effect
!END TYPE ParRecord

TYPE ParRecordTransfer
	SEQUENCE
	INTEGER(lng)	:: parid ! particle id in the overall list - a tag that can be used to track the particle
	INTEGER(lng)	:: cur_part	! current sub-domain id / partition number
	INTEGER(lng)	:: new_part	! current sub-domain id / partition number
	REAL(dbl)	:: xp	! particle x-position
	REAL(dbl)	:: yp ! particle y-position
	REAL(dbl)	:: zp	! particle z-position
	REAL(dbl)	:: up	! particle u-velocity
	REAL(dbl)	:: vp   ! particle v-velocity
	REAL(dbl)	:: wp	! particle w-velocity
	REAL(dbl)	:: rp	! particle radius
	REAL(dbl)	:: delNBbyCV ! particle drug release concentration 
	REAL(dbl)	:: par_conc ! particle concentration
	REAL(dbl)	:: bulk_conc ! bulk concentration at particle location
	REAL(dbl)	:: xpold	! particle x-position
	REAL(dbl)	:: ypold 	! particle y-position
	REAL(dbl)	:: zpold	! particle z-position
	REAL(dbl)	:: upold	! particle u-velocity
	REAL(dbl)	:: vpold   	! particle v-velocity
	REAL(dbl)	:: wpold	! particle w-velocity
	REAL(dbl)	:: rpold	! old particle radius
	REAL(dbl)	:: sh 	! Sherwood number
	REAL(dbl)	:: gamma_cont	! gamma - container effect
	REAL(dbl)	:: S 	! Shear rate at particle location
	REAL(dbl)	:: Sst 	! Shear peclet number
	REAL(dbl)	:: Veff ! effective particle container volume
	REAL(dbl)	:: Nbj  ! number of moles associated with the particlnumber of moles associated with the particle
END TYPE ParRecordTransfer

TYPE ParRecord
	TYPE(ParRecord), POINTER :: prev => NULL()! pointer to prev record
	TYPE(ParRecord), POINTER :: next => NULL()	! pointer to next record
	INTEGER(lng)	:: parid ! particle id in the overall list - a tag that can be used to track the particle
	TYPE(ParRecordTransfer) :: pardata
END TYPE ParRecord

!TYPE(ParRecord),ALLOCATABLE	:: ParList(:)
TYPE(ParRecordTransfer),ALLOCATABLE	:: ParSendArray(:,:),ParRecvArray(:,:)
TYPE(ParRecord), POINTER	:: ParListHead,ParListEnd
TYPE(ParRecordTransfer) :: ParInit
LOGICAL :: ParticleTransfer
INTEGER :: mpipartransfertype
INTEGER :: numparticlesSub
INTEGER(lng), PARAMETER :: der_type_count = 26_lng,numparticlesDomain = 1000_lng
INTEGER :: mpidblextent,mpiintextent
INTEGER(lng), DIMENSION(der_type_count) :: der_block_len,der_block_types,der_block_offsets
REAL(dbl) :: fmovingsum,fmovingrhosum
INTEGER(lng), ALLOCATABLE 	:: parfilenum(:),numparticleSubfile(:) ! array of particle output file numbers and number of particles in each of these files
INTEGER(lng)	:: parfileCount				! current output file number (out of total number of output files)


!************************************************

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE Global_Setup		! sets up simulation
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

CALL ReadInput					! read input from file
!CALL SubDomainSetup			! set up the MPI subdomains
CALL SubDomainSetupNew			! set up the MPI subdomains
CALL AllocateArrays			! allocate global variable arrays

!------------------------------------------------
END SUBROUTINE Global_Setup
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ReadInput			! read the input file
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Read input from input file
OPEN(10,FILE='input.txt')
READ(10,*) domaintype	 				! a flag to denote domain type - 0 for 1/4th cylinder and 1 for full cylinder
READ(10,*) nx	 				! number of nodes in the x-direction
READ(10,*) ny					! number of nodes in the y-direction
READ(10,*) nz					! number of nodes in the z-direction

READ(10,*) NumSubsX			! number of subdomains in the X direction
READ(10,*) NumSubsY			! number of subdomains in the Y direction
READ(10,*) NumSubsZ			! number of subdomains in the Z direction

READ(10,*) L					! length
READ(10,*) D					! diameter

READ(10,*) epsOVERa1			! peristaltic occlusion ratio (distance of occlusion/mean half-width)
READ(10,*) s1					! peristaltic wave speed
READ(10,*) numw1				! number of peristaltic waves
READ(10,*) wc1					! peristaltic weighting coefficient

READ(10,*) epsOVERa2			! segmental occlusion ratio (distance of occlusion/mean half-width)
READ(10,*) Ts					! segmental contraction period
READ(10,*) numw2				! number of segmental waves
READ(10,*) wc2					! segmental weighting coefficient

READ(10,*) Tmix				! period of mixed mode simulation

READ(10,*) numVilliZ			! number of rows of villi in the Z-direction
READ(10,*) numVilliTheta	! number of villi in each row (theta-direction)
READ(10,*) numVilliGroups	! number of villi in each row (theta-direction)

READ(10,*) Lv					! length of the villi (micrometers)
READ(10,*) Rv					! radius of the villi (micrometers)

READ(10,*) freqRatioT		! villous frequency to macroscopic contraction frequency (azimuthal, theta direction)
READ(10,*) freqRatioZ		! villous frequency to macroscopic contraction frequency (axial, z direction)

READ(10,*) randORord			! flag to determine if the villous motion is random or ordered
READ(10,*) villiAngle		! maximum angle of active villous travel (degrees) 

READ(10,*) den					! density
READ(10,*) nu					! kinematic viscosity
READ(10,*) tau					! relaxation parameter

READ(10,*) Sc					! Schmidt number
READ(10,*) sclrIC				! initial/maintained scalar distribution (1=BLOB,2=LINE,3=INLET,4=UNIFORM)
READ(10,*) phiPer				! period at which to start the scalar
READ(10,*) phiIC				! maximum scalar concentration

READ(10,*) nPers				! total number of periods to run
READ(10,*) numOuts			! number of output files (roughly)
READ(10,*) restart			! use restart file? (0 if no, 1 if yes)
READ(10,*) ParticleTrack		! A flag to indicate if particle is on or off (0 if off, 1 if on)
CLOSE(10)

!tau=1.0_dbl

! Check to make sure the number of processors (numprocs) and the number of subdomains are equal
NumSubsTotal = NumSubsX*NumSubsY*NumSubsZ
IF(NumSubsTotal .NE. numprocs) THEN
  OPEN(1000,FILE="error.dat")
  WRITE(1000,*) 'NumSubsTotal is .NE. to numprocs'
  WRITE(1000,*) 'NumSubsTotal', NumSubsTotal
  WRITE(1000,*) 'NumSubsX', NumSubsX
  WRITE(1000,*) 'NumSubsY', NumSubsY
  WRITE(1000,*) 'NumSubsZ', NumSubsZ
  WRITE(1000,*) 'numprocs', numprocs
  WRITE(1000,*) 'check input file and PBS queueing file...'
  CLOSE(1000)
  STOP
END IF

! Check to make sure the number of villi in the axial direction is a multiple of the number of villi groups
IF(MOD(numVilliZ,numVilliGroups) .NE. 0) THEN
  OPEN(1000,FILE="error.dat")
  WRITE(1000,*) 'numVilliZ is not a multiple of numVilliGroups'
  WRITE(1000,*) 'check input file...'
  CLOSE(1000)
  STOP
END IF

!------------------------------------------------
END SUBROUTINE ReadInput
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE SubDomainSetupNew	! generates the information (ID number, starting/ending indices) of each neighboring subdomain (using the subroutine SetSubIDBC in this module)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: CDx(NumCommDirs), CDy(NumCommDirs), CDz(NumCommDirs)		! communication direction vectors in the x, y, and z directions respectively
INTEGER(lng) :: thisSub																	! ID of the current subdomain
INTEGER(lng) :: iComm,iSub,jSub,kSub,iiSub,jjSub,kkSub						! index variables
INTEGER(lng) :: quotientX, quotientY, quotientZ					 				! variables for determining the local subdomain bounds

ALLOCATE(SubID(NumCommDirs))															! id number of neighboring subdomains (same as rank of processing unit working on domain)

ALLOCATE(iMaxDomain(NumSubsTotal))
ALLOCATE(iMinDomain(NumSubsTotal))
ALLOCATE(jMaxDomain(NumSubsTotal))
ALLOCATE(jMinDomain(NumSubsTotal))
ALLOCATE(kMaxDomain(NumSubsTotal))
ALLOCATE(kMinDomain(NumSubsTotal))

! fill out communication direction vectors
CDx(1) =   1_lng
CDy(1) =   0_lng
CDz(1) =   0_lng

CDx(2) =  -1_lng
CDy(2) =   0_lng
CDz(2) =   0_lng

CDx(3) =   0_lng
CDy(3) =   1_lng
CDz(3) =   0_lng

CDx(4) =   0_lng
CDy(4) =  -1_lng
CDz(4) =   0_lng

CDx(5) =   0_lng
CDy(5) =   0_lng
CDz(5) =   1_lng

CDx(6) =   0_lng
CDy(6) =   0_lng
CDz(6) =  -1_lng

CDx(7) =   1_lng
CDy(7) =   1_lng
CDz(7) =   0_lng

CDx(8) =  -1_lng
CDy(8) =  -1_lng
CDz(8) =   0_lng

CDx(9) =   1_lng
CDy(9) =  -1_lng
CDz(9) =   0_lng

CDx(10) = -1_lng
CDy(10) =  1_lng
CDz(10) =  0_lng

CDx(11) =  0_lng
CDy(11) =  1_lng
CDz(11) =  1_lng

CDx(12) =  0_lng
CDy(12) = -1_lng
CDz(12) = -1_lng

CDx(13) =  0_lng
CDy(13) =  1_lng
CDz(13) = -1_lng

CDx(14) =  0_lng
CDy(14) = -1_lng
CDz(14) =  1_lng

CDx(15) =  1_lng
CDy(15) =  0_lng
CDz(15) =  1_lng

CDx(16) = -1_lng
CDy(16) =  0_lng
CDz(16) = -1_lng

CDx(17) = -1_lng
CDy(17) =  0_lng
CDz(17) =  1_lng

CDx(18) =  1_lng
CDy(18) =  0_lng
CDz(18) = -1_lng

CDx(19) =  1_lng
CDy(19) =  1_lng
CDz(19) =  1_lng

CDx(20) = -1_lng
CDy(20) = -1_lng
CDz(20) = -1_lng

CDx(21) =  1_lng
CDy(21) =  1_lng
CDz(21) = -1_lng

CDx(22) = -1_lng
CDy(22) = -1_lng
CDz(22) =  1_lng

CDx(23) = -1_lng
CDy(23) =  1_lng
CDz(23) =  1_lng

CDx(24) =  1_lng
CDy(24) = -1_lng
CDz(24) = -1_lng

CDx(25) =  1_lng
CDy(25) = -1_lng
CDz(25) =  1_lng

CDx(26) = -1_lng
CDy(26) =  1_lng
CDz(26) = -1_lng

! Number of the current subdomain
mySub = myid + 1_lng							! subdomain number
!WRITE(sub(1:2),'(I2.2)') mySub			! write subdomain number to 'sub' for output file exentsions
WRITE(sub(1:5),'(I5.5)') mySub			! write subdomain number to 'sub' for output file exentsions

! Define the local computational domain bounds (iMin:iMax,jMin:jMax,kMin:kMax)
quotientX	= CEILING(REAL(nx)/NumSubsX)						! divide the number of nodes by the number of subdomains (round up)
quotientY	= CEILING(REAL(ny)/NumSubsY)						! divide the number of nodes by the number of subdomains (round up)
quotientZ	= CEILING(REAL(nz)/NumSubsZ)						! divide the number of nodes by the number of subdomains (round up)

iMin = MOD(myid,NumSubsX)*quotientX + 1_lng					! starting local i index 
iMax = iMin + (quotientX - 1_lng)								! ending local i index

jMin = MOD((myid/NumSubsX),NumSubsY)*quotientY + 1_lng	! starting local j index
jMax = jMin + (quotientY - 1_lng)								! ending local j index

kMin = (myid/(NumSubsX*NumSubsY))*quotientZ + 1_lng		! starting local k index 
kMax = kMin + (quotientZ - 1_lng)								! ending local k index

! Check the bounds
IF(iMax .GT. nx) THEN
  iMax = nx																! if iMax is greater than nx, correct it
END IF

IF(jMax .GT. ny) THEN
  jMax = ny																! if jMax is greater than ny, correct it
END IF

IF(kMax .GT. nz) THEN
  kMax = nz																! if kMax is greater than nz, correct it
END IF


! Loop through the subdomains
DO kSub=1,NumSubsZ
  DO jSub=1,NumSubsY
    DO iSub=1,NumSubsX

      thisSub = iSub + (jSub-1)*NumSubsX + (kSub-1)*NumSubsX*NumSubsY	! get the ID of the current Subdomain

      IF(mySub .EQ. thisSub) THEN													! fill out the SubID array of the current subdomain is the 
     
        ! Loop through the communication directions for the current subdomain
        DO iComm=1,NumCommDirs

          iiSub = iSub + CDx(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
          jjSub = jSub + CDy(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
          kkSub = kSub + CDz(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
      
          !CALL SetSubID(iComm,iiSub,jjSub,kkSub)								! identify the neighboring subdomains (SubID)
          CALL SetSubIDNew(iComm,iiSub,jjSub,kkSub)								! identify the neighboring subdomains (SubID)

        END DO

      END IF

! Fill up arrays containging iMax, iMin, jMax,jMin,kMax, kMin for all subdomains


	iMinDomain(thisSub) = MOD((thisSub-1_lng),NumSubsX)*quotientX + 1_lng			! starting local i index 
	iMaxDomain(thisSub) = iMinDomain(thisSub) + (quotientX - 1_lng)				! ending local i index
	
	jMinDomain(thisSub) = MOD(((thisSub-1_lng)/NumSubsX),NumSubsY)*quotientY + 1_lng	! starting local j index
	jMaxDomain(thisSub) = jMinDomain(thisSub) + (quotientY - 1_lng)				! ending local j index
	
	kMinDomain(thisSub) = ((thisSub-1_lng)/(NumSubsX*NumSubsY))*quotientZ + 1_lng		! starting local k index 
	kMaxDomain(thisSub) = kMinDomain(thisSub) + (quotientZ - 1_lng)				! ending local k index
	
	! Check the bounds
	IF(iMaxDomain(thisSub) .GT. nx) THEN
	  iMaxDomain(thisSub) = nx																! if iMax is greater than nx, correct it
	END IF
	
	IF(jMaxDomain(thisSub) .GT. ny) THEN
	  jMaxDomain(thisSub) = ny																! if jMax is greater than ny, correct it
	END IF
	
	IF(kMaxDomain(thisSub) .GT. nz) THEN
	  kMaxDomain(thisSub) = nz																! if kMax is greater than nz, correct it
	END IF

    END DO
  END DO
END DO

! Commented out by Balaji 02/25/2015
!! Define the local computational domain bounds (iMin:iMax,jMin:jMax,kMin:kMax)
!quotientX	= CEILING(REAL(nx)/NumSubsX)						! divide the number of nodes by the number of subdomains (round up)
!quotientY	= CEILING(REAL(ny)/NumSubsY)						! divide the number of nodes by the number of subdomains (round up)
!quotientZ	= CEILING(REAL(nz)/NumSubsZ)						! divide the number of nodes by the number of subdomains (round up)
!
!iMin = MOD(myid,NumSubsX)*quotientX + 1_lng					! starting local i index 
!iMax = iMin + (quotientX - 1_lng)								! ending local i index
!
!jMin = MOD((myid/NumSubsX),NumSubsY)*quotientY + 1_lng	! starting local j index
!jMax = jMin + (quotientY - 1_lng)								! ending local j index
!
!kMin = (myid/(NumSubsX*NumSubsY))*quotientZ + 1_lng		! starting local k index 
!kMax = kMin + (quotientZ - 1_lng)								! ending local k index
!
!! Check the bounds
!IF(iMax .GT. nx) THEN
!  iMax = nx																! if iMax is greater than nx, correct it
!END IF
!
!IF(jMax .GT. ny) THEN
!  jMax = ny																! if jMax is greater than ny, correct it
!END IF
!
!IF(kMax .GT. nz) THEN
!  kMax = nz																! if kMax is greater than nz, correct it
!END IF
!

! Determine the number of nodes in each direction
nxSub = (iMax - iMin) + 1_lng
nySub = (jMax - jMin) + 1_lng
nzSub = (kMax - kMin) + 1_lng

! Write the local bounds to a file [TEST]
!OPEN(171,FILE='localBounds-'//sub//'.dat')
!WRITE(171,*) 'iMin =', iMin, 'iMax=', iMax
!WRITE(171,*) 'jMin =', jMin, 'jMax=', jMax 
!WRITE(171,*) 'kMin =', kMin, 'kMax=', kMax
!WRITE(171,*) 
!WRITE(171,*) 'nx =', nx, 'ny=', ny, 'nz=', nz
!WRITE(171,*) 'nxSub =', nxSub, 'nySub=', nySub, 'nzSub=', nzSub
!WRITE(171,*) 
!WRITE(171,*) 'quotientX =', quotientX
!WRITE(171,*) 'quotientY =', quotientY
!WRITE(171,*) 'quotientZ =', quotientZ
!CLOSE(171)

! Write the subID to a file [TEST]
!OPEN(172,FILE='sub-'//sub//'.dat')
!DO iComm=1,NumCommDirs
!  WRITE(172,*)'iComm=',iComm,'SubID(iComm)=', SubID(iComm)
!END DO
!CLOSE(172)
!STOP

!------------------------------------------------
END SUBROUTINE SubDomainSetupNew
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE SubDomainSetup	! generates the information (ID number, starting/ending indices) of each neighboring subdomain (using the subroutine SetSubIDBC in this module)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: CDx(NumCommDirs), CDy(NumCommDirs), CDz(NumCommDirs)		! communication direction vectors in the x, y, and z directions respectively
INTEGER(lng) :: thisSub																	! ID of the current subdomain
INTEGER(lng) :: iComm,iSub,jSub,kSub,iiSub,jjSub,kkSub						! index variables
INTEGER(lng) :: quotientX, quotientY, quotientZ					 				! variables for determining the local subdomain bounds

ALLOCATE(SubID(NumCommDirs))															! id number of neighboring subdomains (same as rank of processing unit working on domain)

! fill out communication direction vectors
CDx(1) =   1_lng
CDy(1) =   0_lng
CDz(1) =   0_lng

CDx(2) =  -1_lng
CDy(2) =   0_lng
CDz(2) =   0_lng

CDx(3) =   0_lng
CDy(3) =   1_lng
CDz(3) =   0_lng

CDx(4) =   0_lng
CDy(4) =  -1_lng
CDz(4) =   0_lng

CDx(5) =   0_lng
CDy(5) =   0_lng
CDz(5) =   1_lng

CDx(6) =   0_lng
CDy(6) =   0_lng
CDz(6) =  -1_lng

CDx(7) =   1_lng
CDy(7) =   1_lng
CDz(7) =   0_lng

CDx(8) =  -1_lng
CDy(8) =  -1_lng
CDz(8) =   0_lng

CDx(9) =   1_lng
CDy(9) =  -1_lng
CDz(9) =   0_lng

CDx(10) = -1_lng
CDy(10) =  1_lng
CDz(10) =  0_lng

CDx(11) =  0_lng
CDy(11) =  1_lng
CDz(11) =  1_lng

CDx(12) =  0_lng
CDy(12) = -1_lng
CDz(12) = -1_lng

CDx(13) =  0_lng
CDy(13) =  1_lng
CDz(13) = -1_lng

CDx(14) =  0_lng
CDy(14) = -1_lng
CDz(14) =  1_lng

CDx(15) =  1_lng
CDy(15) =  0_lng
CDz(15) =  1_lng

CDx(16) = -1_lng
CDy(16) =  0_lng
CDz(16) = -1_lng

CDx(17) = -1_lng
CDy(17) =  0_lng
CDz(17) =  1_lng

CDx(18) =  1_lng
CDy(18) =  0_lng
CDz(18) = -1_lng

CDx(19) =  1_lng
CDy(19) =  1_lng
CDz(19) =  1_lng

CDx(20) = -1_lng
CDy(20) = -1_lng
CDz(20) = -1_lng

CDx(21) =  1_lng
CDy(21) =  1_lng
CDz(21) = -1_lng

CDx(22) = -1_lng
CDy(22) = -1_lng
CDz(22) =  1_lng

CDx(23) = -1_lng
CDy(23) =  1_lng
CDz(23) =  1_lng

CDx(24) =  1_lng
CDy(24) = -1_lng
CDz(24) = -1_lng

CDx(25) =  1_lng
CDy(25) = -1_lng
CDz(25) =  1_lng

CDx(26) = -1_lng
CDy(26) =  1_lng
CDz(26) = -1_lng

! Number of the current subdomain
mySub = myid + 1_lng							! subdomain number
!WRITE(sub(1:2),'(I2.2)') mySub			! write subdomain number to 'sub' for output file exentsions
WRITE(sub(1:5),'(I5.5)') mySub			! write subdomain number to 'sub' for output file exentsions

! Loop through the subdomains
DO kSub=1,NumSubsZ
  DO jSub=1,NumSubsY
    DO iSub=1,NumSubsX

      thisSub = iSub + (jSub-1)*NumSubsX + (kSub-1)*NumSubsX*NumSubsY	! get the ID of the current Subdomain

      IF(mySub .EQ. thisSub) THEN													! fill out the SubID array of the current subdomain is the 
     
        ! Loop through the communication directions for the current subdomain
        DO iComm=1,NumCommDirs

          iiSub = iSub + CDx(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
          jjSub = jSub + CDy(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
          kkSub = kSub + CDz(iComm)													! subdomain index of neighboring subdomain in the iCommth communication direction
      
          CALL SetSubID(iComm,iiSub,jjSub,kkSub)								! identify the neighboring subdomains (SubID)

        END DO

      END IF

    END DO
  END DO
END DO

! Define the local computational domain bounds (iMin:iMax,jMin:jMax,kMin:kMax)
quotientX	= CEILING(REAL(nx)/NumSubsX)						! divide the number of nodes by the number of subdomains (round up)
quotientY	= CEILING(REAL(ny)/NumSubsY)						! divide the number of nodes by the number of subdomains (round up)
quotientZ	= CEILING(REAL(nz)/NumSubsZ)						! divide the number of nodes by the number of subdomains (round up)

iMin = MOD(myid,NumSubsX)*quotientX + 1_lng					! starting local i index 
iMax = iMin + (quotientX - 1_lng)								! ending local i index

jMin = MOD((myid/NumSubsX),NumSubsY)*quotientY + 1_lng	! starting local j index
jMax = jMin + (quotientY - 1_lng)								! ending local j index

kMin = (myid/(NumSubsX*NumSubsY))*quotientZ + 1_lng		! starting local k index 
kMax = kMin + (quotientZ - 1_lng)								! ending local k index

! Check the bounds
IF(iMax .GT. nx) THEN
  iMax = nx																! if iMax is greater than nx, correct it
END IF

IF(jMax .GT. ny) THEN
  jMax = ny																! if jMax is greater than ny, correct it
END IF

IF(kMax .GT. nz) THEN
  kMax = nz																! if kMax is greater than nz, correct it
END IF

! Determine the number of nodes in each direction
nxSub = (iMax - iMin) + 1_lng
nySub = (jMax - jMin) + 1_lng
nzSub = (kMax - kMin) + 1_lng

! Write the local bounds to a file [TEST]
!OPEN(171,FILE='localBounds-'//sub//'.dat')
!WRITE(171,*) 'iMin =', iMin, 'iMax=', iMax
!WRITE(171,*) 'jMin =', jMin, 'jMax=', jMax 
!WRITE(171,*) 'kMin =', kMin, 'kMax=', kMax
!WRITE(171,*) 
!WRITE(171,*) 'nx =', nx, 'ny=', ny, 'nz=', nz
!WRITE(171,*) 'nxSub =', nxSub, 'nySub=', nySub, 'nzSub=', nzSub
!WRITE(171,*) 
!WRITE(171,*) 'quotientX =', quotientX
!WRITE(171,*) 'quotientY =', quotientY
!WRITE(171,*) 'quotientZ =', quotientZ
!CLOSE(171)

! Write the subID to a file [TEST]
!OPEN(172,FILE='sub-'//sub//'.dat')
!DO iComm=1,NumCommDirs
!  WRITE(172,*)'iComm=',iComm,'SubID(iComm)=', SubID(iComm)
!END DO
!CLOSE(172)
!STOP

!------------------------------------------------
END SUBROUTINE SubDomainSetup
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE SetSubIDNew(iComm,iiSub,jjSub,kkSub)									! sets SubID based on neighboring subdomains
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng), INTENT(IN) :: iComm, iiSub, jjSub, kkSub 					! index variables
INTEGER(lng) :: nSub, kkSub2,jjSub2														! neighboring subdomain ID, kkSub (reset for periodicity)

IF(((iiSub .LT. 1) .OR. (iiSub .GT. NumSubsX))		&
!   .OR. ((kkSub .LT. 1) .OR. (kkSub .GT. NumSubsZ))	& 					! comment out for periodic BCs in the k-direction
!   .OR. ((jjSub .LT. 1) .OR. (jjSub .GT. NumSubsY))	&
	)	THEN 

  SubID(iComm) = 0_lng	! no neighbor

ELSE IF((kkSub .LT. 1)) THEN
	kkSub2 = NumSubsZ	! reset kkSub for periodicity in the z-direction
	IF((jjSub .LT. 1)) THEN
		  jjSub2 = NumSubsY																	! reset kkSub for periodicity in the z-direction
		  nSub = iiSub + (jjSub2-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		  SubID(iComm) = nSub			! set SubID(iComm) to neighboring sudomain ID
	
	ELSE IF((jjSub .GT. NumSubsY)) THEN
		jjSub2 = 1_lng																		! reset kkSub for periodicity in the z-direction
		nSub = iiSub + (jjSub2-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		SubID(iComm) = nSub			! set SubID(iComm) to neighboring sudomain ID
	ELSE
		nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		SubID(iComm) = nSub		! set SubID(iComm) to neighboring sudomain ID
	END IF
ELSE IF((kkSub .GT. NumSubsZ)) THEN

	kkSub2 = 1_lng	! reset kkSub for periodicity in the z-direction
	IF((jjSub .LT. 1)) THEN
		  jjSub2 = NumSubsY																	! reset kkSub for periodicity in the z-direction
		  nSub = iiSub + (jjSub2-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		  SubID(iComm) = nSub			! set SubID(iComm) to neighboring sudomain ID
	
	ELSE IF((jjSub .GT. NumSubsY)) THEN
		jjSub2 = 1_lng																		! reset kkSub for periodicity in the z-direction
		nSub = iiSub + (jjSub2-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		SubID(iComm) = nSub			! set SubID(iComm) to neighboring sudomain ID
	ELSE
		nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		SubID(iComm) = nSub		! set SubID(iComm) to neighboring sudomain ID
	END IF
ELSE
  
	IF((jjSub .LT. 1)) THEN
		  jjSub2 = NumSubsY																	! reset kkSub for periodicity in the z-direction
		  nSub = iiSub + (jjSub2-1)*NumSubsX + (kkSub-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		  SubID(iComm) = nSub			! set SubID(iComm) to neighboring sudomain ID
	
	ELSE IF((jjSub .GT. NumSubsY)) THEN
		jjSub2 = 1_lng																		! reset kkSub for periodicity in the z-direction
		nSub = iiSub + (jjSub2-1)*NumSubsX + (kkSub-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		SubID(iComm) = nSub			! set SubID(iComm) to neighboring sudomain ID
	ELSE
		nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
		SubID(iComm) = nSub		! set SubID(iComm) to neighboring sudomain ID
	END IF


END IF

!write(*,*) SubID(iComm),NumSubsX,NumSubsY,NumSubsZ,iiSub,jjSub,kkSub,SubID(icomm)

!------------------------------------------------
END SUBROUTINE SetSubIDNew
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE SetSubID(iComm,iiSub,jjSub,kkSub)									! sets SubID based on neighboring subdomains
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng), INTENT(IN) :: iComm, iiSub, jjSub, kkSub 					! index variables
INTEGER(lng) :: nSub, kkSub2														! neighboring subdomain ID, kkSub (reset for periodicity)

IF(((jjSub .LT. 1) .OR. (jjSub .GT. NumSubsY))	.OR.	&
!   ((kkSub .LT. 1) .OR. (kkSub .GT. NumSubsZ))	.OR. 	& 					! comment out for periodic BCs in the k-direction
   ((iiSub .LT. 1) .OR. (iiSub .GT. NumSubsX)))	THEN 

  SubID(iComm) = 0_lng																! no neighbor

ELSE IF((kkSub .LT. 1)) THEN

  kkSub2 = NumSubsZ																	! reset kkSub for periodicity in the z-direction
  nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
  SubID(iComm) = nSub																! set SubID(iComm) to neighboring sudomain ID

ELSE IF((kkSub .GT. NumSubsZ)) THEN

  kkSub2 = 1_lng																		! reset kkSub for periodicity in the z-direction
  nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub2-1)*NumSubsX*NumSubsY	! neighboring subdomain ID
  SubID(iComm) = nSub																! set SubID(iComm) to neighboring sudomain ID
    
ELSE
  
  nSub = iiSub + (jjSub-1)*NumSubsX + (kkSub-1)*NumSubsX*NumSubsY		! neighboring subdomain ID
  SubID(iComm) = nSub																! set SubID(iComm) to neighboring sudomain ID

END IF

!------------------------------------------------
END SUBROUTINE SetSubID
!------------------------------------------------

!-------------------------------------------------------------------------------------------------
!!!!!!! SUBROUTINES TO HANDLE POINTERS AND LINKED LISTS FOR PARTICLES
!-------------------------------------------------------------------------------------------------


! Initialize a head node SELF and optionally store the provided DATA.
!------------------------------------------------
SUBROUTINE list_init(self)
!------------------------------------------------
  TYPE(ParRecord), POINTER :: self

  ALLOCATE(self) ! Note: When self is allocated, all the 
  
!  ALLOCATE(self%next)
!  ALLOCATE(self%prev)
!  ALLOCATE(self%xp)
!  ALLOCATE(self%yp)
!  ALLOCATE(self%zp)
!  ALLOCATE(self%up)
!  ALLOCATE(self%vp)
!  ALLOCATE(self%wp)
!  ALLOCATE(self%rp)
!  ALLOCATE(self%delNBbyCV)
!  ALLOCATE(self%par_conc)
!  ALLOCATE(self%bulk_conc)
!  ALLOCATE(self%rpold)
!  ALLOCATE(self%sh)
!  ALLOCATE(self%gamma_cont)
!  ALLOCATE(self%cur_part)

  NULLIFY(self%next)
  NULLIFY(self%prev)

!------------------------------------------------
END SUBROUTINE list_init
!------------------------------------------------

! Insert a list node after SELF (an arbitrary node)
! NOTE: Remember to assign data to these new nodes
!------------------------------------------------
SUBROUTINE list_insert(self)
!------------------------------------------------
  TYPE(ParRecord), POINTER :: self
  TYPE(ParRecord), POINTER :: new

! Allocate and initialize the new pointer
  ALLOCATE(new)
  NULLIFY(new%next)
  NULLIFY(new%prev)

  IF (ASSOCIATED(self%next)) THEN
	new%next => self%next
	self%next%prev => new
  ELSE 
	NULLIFY(new%next)
  END IF
  new%prev => self
  self%next => new

!------------------------------------------------
END SUBROUTINE list_insert
!------------------------------------------------

! Delete a list node pointed by SELF (an arbitrary node)
!------------------------------------------------
SUBROUTINE list_delete(self)
!------------------------------------------------
  TYPE(ParRecord), POINTER :: self
  TYPE(ParRecord), POINTER :: next

  !ALLOCATE(next)
  self%prev%next => self%next
  IF (ASSOCIATED(self%next)) THEN
	self%next%prev => self%prev
  ENDIF

  !self%prev => NULL()
  !self%next => NULL()
  NULLIFY(self%prev)
  NULLIFY(self%next)
  DEALLOCATE(self)

!------------------------------------------------
END SUBROUTINE list_delete
!------------------------------------------------


! Free the entire list and all data, beginning at SELF
!------------------------------------------------
SUBROUTINE list_free(self)
!------------------------------------------------
  TYPE(ParRecord), POINTER :: self
  TYPE(ParRecord), POINTER :: current
  TYPE(ParRecord), POINTER :: next
  current => self
  DO WHILE (ASSOCIATED(current))
     next => current%next ! copy pointer of next node
     DEALLOCATE(current)
     NULLIFY(current)
     ! point to next node in the list
     current => next
     !write(*,*) i
  END DO
!------------------------------------------------
end subroutine list_free
!------------------------------------------------



!--------------------------------------------------------------------------------------------------
SUBROUTINE AllocateArrays	! allocates array space
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Distribution Functions
ALLOCATE(f(0:NumDistDirs,0:nxSub+1,0:nySub+1,0:nzSub+1),			&
         fplus(0:NumDistDirs,0:nxSub+1,0:nySub+1,0:nzSub+1))
! Velocity, Density
ALLOCATE(u(0:nxSub+1,0:nySub+1,0:nzSub+1),							&
         v(0:nxSub+1,0:nySub+1,0:nzSub+1),							&
         w(0:nxSub+1,0:nySub+1,0:nzSub+1))
ALLOCATE(rho(0:nxSub+1,0:nySub+1,0:nzSub+1))

! Scalar
ALLOCATE(phi(0:nxSub+1,0:nySub+1,0:nzSub+1), 						&
         phiTemp(0:nxSub+1,0:nySub+1,0:nzSub+1))
ALLOCATE(delphi_particle(0:nxSub+1,0:nySub+1,0:nzSub+1))
ALLOCATE(tausgs_particle_x(0:nxSub+1,0:nySub+1,0:nzSub+1))
ALLOCATE(tausgs_particle_y(0:nxSub+1,0:nySub+1,0:nzSub+1))
ALLOCATE(tausgs_particle_z(0:nxSub+1,0:nySub+1,0:nzSub+1))

! Node Flags
ALLOCATE(node(0:nxSub+1,0:nySub+1,0:nzSub+1))

! LBM Miscellaneous
ALLOCATE(ex(0:NumDistDirs),ey(0:NumDistDirs),ez(0:NumDistDirs))
ALLOCATE(bb(0:NumDistDirs),sym(0:NumDistDirs,0:NumDistDirs))
ALLOCATE(wt(0:NumDistDirs))

! MPI Communication Arrays
ALLOCATE(f_Comps(NumCommDirs,MaxDistFns))					! specifies the components of the distribution functions to transfer in each MPI communication direction
ALLOCATE(Corner_SendIndex(19:26,3))							! i, j, and k indices for each corner
ALLOCATE(Corner_RecvIndex(19:26,3))							! i, j, and k indices for each corner (phantom node for recieving data)
ALLOCATE(Z_SendIndex(7:10,2))									! i and j indices for each Z side 
ALLOCATE(Z_RecvIndex(7:10,2))									! i and j indices for each Z side (phantom node for recieving data)
ALLOCATE(X_SendIndex(11:14,2))								! j and k indices for each X side 
ALLOCATE(X_RecvIndex(11:14,2))								! j and k indices for each X side (phantom node for recieving data)
ALLOCATE(Y_SendIndex(15:18,2))								! i and k indices for each Y side 
ALLOCATE(Y_RecvIndex(15:18,2))								! i and k indices for each Y side (phantom node for recieving data)
ALLOCATE(YZ_SendIndex(1:2))									! i index for each YZ face 
ALLOCATE(YZ_RecvIndex(1:2))									! i index for each YZ face (phantom node for recieving data)
ALLOCATE(ZX_SendIndex(3:4))									! j index for each ZX face 
ALLOCATE(ZX_RecvIndex(3:4))									! j index for each ZX face (phantom node for recieving data)
ALLOCATE(XY_SendIndex(5:6))									! k index for each XY face 
ALLOCATE(XY_RecvIndex(5:6))									! k index for each XY face (phantom node for recieving data)
ALLOCATE(OppCommDir(NumCommDirs)) 							! opposite MPI communication directions (like bounceback) 
ALLOCATE(CommDataStart_f(NumCommDirs))						! array of starting indices in the send arrays for the distribution functions from each communication direction 
ALLOCATE(CommDataStart_rho(NumCommDirs))					! array of starting indices in the send arrays for the density from each communication direction
ALLOCATE(CommDataStart_phi(NumCommDirs))					! array of starting indices in the send arrays for the scalar from each communication direction
ALLOCATE(CommDataStart_u(NumCommDirs))					! array of starting indices in the send arrays for the scalar from each communication direction
ALLOCATE(CommDataStart_v(NumCommDirs))					! array of starting indices in the send arrays for the scalar from each communication direction
ALLOCATE(CommDataStart_w(NumCommDirs))					! array of starting indices in the send arrays for the scalar from each communication direction
ALLOCATE(CommDataStart_node(NumCommDirs))					! array of starting indices in the send arrays for the node index from each communication direction
ALLOCATE(fSize(NumCommDirs))									! array of the number of elements sent for each communication direction (distribution functions)
ALLOCATE(dsSize(NumCommDirs))									! array of the number of elements sent for each communication direction (density and scalar)
ALLOCATE(uvwSize(NumCommDirs))									! array of the number of elements sent for each communication direction (density and scalar)
ALLOCATE(nodeSize(NumCommDirs))									! array of the number of elements sent for each communication direction (density and scalar)
ALLOCATE(msgSize(NumCommDirs))								! array of the number of elements sent for each communication direction (total)
ALLOCATE(req(2*NumCommDirs))									! allocate the MPI send request array

! Geometry Arrays
!ALLOCATE(rDom0(0:nz+1),rDom(0:nz+1),r(0:nzSub+1))		! intial and current radius (global), current radius (local)
ALLOCATE(rDom0In(0:nz+1),rDomIn(0:nz+1),rIn(0:nzSub+1))		! intial and current radius (global), current radius (local)
ALLOCATE(rDom0Out(0:nz+1),rDomOut(0:nz+1),rOut(0:nzSub+1))	! intial and current radius (global), current radius (local)
!ALLOCATE(velDom(0:nz+1),vel(0:nzSub+1))					! global and local wall velocities
ALLOCATE(velDomIn(0:nz+1),velIn(0:nzSub+1))					! global and local wall velocities
ALLOCATE(velDomOut(0:nz+1),velOut(0:nzSub+1))					! global and local wall velocities
ALLOCATE(x(0:nxSub+1),y(0:nySub+1),z(0:nzSub+1))		! x, y, z, physical coordinate arrays (local)
ALLOCATE(xx(0:nx+1),yy(0:ny+1),zz(0:nz+1))				! x, y, z, physical coordinate arrays (global)

!------------------------------------------------
END SUBROUTINE AllocateArrays
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE DEAllocateArrays	! allocates array space
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Distribution Functions
DEALLOCATE(f,fplus)

! Velocity, Density
DEALLOCATE(u,v,w,rho)

! Scalar
DEALLOCATE(phi,phiTemp,delphi_particle,tausgs_particle_x,tausgs_particle_y,tausgs_particle_z)

! Node Flags
DEALLOCATE(node)

! LBM Miscellaneous
DEALLOCATE(ex,ey,ez)
DEALLOCATE(bb,sym)
DEALLOCATE(wt)

! MPI Communication Arrays
DEALLOCATE(f_Comps)	! specifies the components of the distribution functions to transfer in each MPI communication direction
DEALLOCATE(Corner_SendIndex)		! i, j, and k indices for each corner
DEALLOCATE(Corner_RecvIndex)		! i, j, and k indices for each corner (phantom node for recieving data)
DEALLOCATE(Z_SendIndex)				! i and j indices for each Z side 
DEALLOCATE(Z_RecvIndex)				! i and j indices for each Z side (phantom node for recieving data)
DEALLOCATE(X_SendIndex)				! j and k indices for each X side 
DEALLOCATE(X_RecvIndex)				! j and k indices for each X side (phantom node for recieving data)
DEALLOCATE(Y_SendIndex)				! i and k indices for each Y side 
DEALLOCATE(Y_RecvIndex)				! i and k indices for each Y side (phantom node for recieving data)
DEALLOCATE(YZ_SendIndex)			! i index for each YZ face 
DEALLOCATE(YZ_RecvIndex)			! i index for each YZ face (phantom node for recieving data)
DEALLOCATE(ZX_SendIndex)			! j index for each ZX face 
DEALLOCATE(ZX_RecvIndex)			! j index for each ZX face (phantom node for recieving data)
DEALLOCATE(XY_SendIndex)			! k index for each XY face 
DEALLOCATE(XY_RecvIndex)			! k index for each XY face (phantom node for recieving data)
DEALLOCATE(OppCommDir) 				! opposite MPI communication directions (like bounceback) 
DEALLOCATE(CommDataStart_f)		! array of starting indices in the send arrays for the distribution functions from each communication direction 
DEALLOCATE(CommDataStart_rho)		! array of starting indices in the send arrays for the density from each communication direction
DEALLOCATE(CommDataStart_phi)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(CommDataStart_u)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(CommDataStart_v)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(CommDataStart_w)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(CommDataStart_node)		! array of starting indices in the send arrays for the scalar from each communication direction
DEALLOCATE(fSize)						! array of the number of elements sent for each communication direction (distribution functions)
DEALLOCATE(dsSize)					! array of the number of elements sent for each communication direction (density and scalar)
DEALLOCATE(uvwSize)					! array of the number of elements sent for each communication direction (density and scalar)
DEALLOCATE(nodeSize)					! array of the number of elements sent for each communication direction (density and scalar)
DEALLOCATE(msgSize)					! array of the number of elements sent for each communication direction (density and scalar)
DEALLOCATE(req)						! array of MPI send/receive requests
DEALLOCATE(waitStat)					! array of MPI_WAITALL requests
DEALLOCATE(SubID)

! Geometry Arrays
!DEALLOCATE(rDom0,rDom,r)			! intial and current radius (global), current radius (local)
DEALLOCATE(rDom0In,rDomIn,rIn)			! intial and current radius (global), current radius (local)
DEALLOCATE(rDom0Out,rDomOut,rOut)			! intial and current radius (global), current radius (local)
!DEALLOCATE(velDom,vel)				! global and local wall velocities
DEALLOCATE(velDomIn,velIn)				! global and local wall velocities
DEALLOCATE(velDomOut,velOut)				! global and local wall velocities
DEALLOCATE(x,y,z)						! x, y, z, physical coordinate arrays (local)
DEALLOCATE(villiLoc)					! location of the villi
IF(randORord .EQ. RANDOM) THEN
  DEALLOCATE(rnd)						! array of random numbers for random villi phase angles
END IF

!Particle arrays
IF(ParticleTrack.EQ.ParticleOn) THEN
	!DEALLOCATE(xp,yp,zp,up,vp,wp,ipar,jpar,kpar,rp,delNBbyCV)
	!DEALLOCATE(par_conc,bulk_conc,sh,gamma_cont,rpold)
	!DEALLOCATE(ParList)
	CALL list_free(ParListHead)
END IF
! Particle MPI arrays
DEALLOCATE(iMinDomain)
DEALLOCATE(iMaxDomain)
DEALLOCATE(jMinDomain)
DEALLOCATE(jMaxDomain)
DEALLOCATE(kMinDomain)
DEALLOCATE(kMaxDomain)
DEALLOCATE(partransfersend)
DEALLOCATE(partransferrecv)
DEALLOCATE(numpartransfer)
DEALLOCATE(parreqid)
DEALLOCATE(parwtstat)
DEALLOCATE(probestat)
!DEALLOCATE(ParSendArray)
!DEALLOCATE(ParRecvArray)
!------------------------------------------------
END SUBROUTINE DEAllocateArrays
!------------------------------------------------

!================================================
END MODULE Setup
!================================================
