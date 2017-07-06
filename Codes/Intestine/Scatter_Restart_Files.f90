IMPLICIT NONE 

LOGICAL :: Flag_Convection_Effects     ! Flag for including convection effects in Sherwood number
CHARACTER(7) :: iter_char                                  ! iteration stored as a character
CHARACTER(5) :: sub
INTEGER :: iter0,i,j,k,m                                    ! index variables
INTEGER :: nx,ny,nz                                   ! number of global nodes in the x, y, and z directions respectively
INTEGER :: nxSub,nySub,nzSub					                 ! number of local nodes in the each direction
INTEGER :: NumSubsX,NumSubsY,NumSubsZ					    	 ! number of local nodes in the each direction
INTEGER :: iMin,iMax,jMin,jMax,kMin,kMax
INTEGER :: quotient_X,quotient_Y,quotient_Z
INTEGER :: CPU,N_CPU
INTEGER, PARAMETER :: NumDistDirs= 14                 ! number of distribution function directions minus one (ex. D3Q15 -> 14)
REAL*8  :: tmp
REAL*8  :: Drug_Initial
REAL*8  :: Drug_Released
REAL*8  :: Drug_Absorbed 
REAL*8  :: Drug_Remained 
REAL*8  :: delphi_particle 
REAL*8,	ALLOCATABLE :: f(:,:,:,:)                    ! distribution function
REAL*8,	ALLOCATABLE :: phi(:,:,:)                    ! distribution function
REAL*8,	ALLOCATABLE :: u(:,:,:),v(:,:,:),w(:,:,:)    ! x,y, and z components of the fluid velocity vector
REAL*8,	ALLOCATABLE :: DUdt_x(:,:,:),DUdt_y(:,:,:),DUdt_z(:,:,:) 
REAL*8,	ALLOCATABLE :: Laplacian_x(:,:,:),Laplacian_y(:,:,:),Laplacian_z(:,:,:)
REAL*8,	ALLOCATABLE :: DLaplacianDt_x(:,:,:),DLaplacianDt_y(:,:,:),DLaplacianDt_z(:,:,:)
REAL*8, ALLOCATABLE :: rho(:,:,:)                    ! density
INTEGER,ALLOCATABLE :: node(:,:,:)                   ! node flags (FLUID/SOLID)

OPEN(10,FILE='input.txt')
READ(10,*) tmp        ! a flag to denote domain type - 0 for 1/4th cylinder and 1 for full cylinder
READ(10,*) nx	 				! number of nodes in the x-direction
READ(10,*) ny					! number of nodes in the y-direction
READ(10,*) nz					! number of nodes in the z-direction
READ(10,*) NumSubsX		! number of subdomains in the X direction
READ(10,*) NumSubsY		! number of subdomains in the Y direction
READ(10,*) NumSubsZ		! number of subdomains in the Z direction
READ(10,*)  !Width      ! Width (only in case of Couette simulation)
READ(10,*)  !D					 ! diameter
READ(10,*)  !L					  ! length
READ(10,*)  !epsOVERa1		! peristaltic occlusion ratio (distance of occlusion/mean half-width)
READ(10,*)  !s1					  ! peristaltic wave speed
READ(10,*)  !s_movingF    ! Moving Frame of Reference speed
READ(10,*)  !numw1				! number of peristaltic waves
READ(10,*)  !wc1					! peristaltic weighting coefficient
READ(10,*)  !epsOVERa2		! segmental occlusion ratio (distance of occlusion/mean half-width)
READ(10,*)  !Ts					  ! segmental contraction period
READ(10,*)  !numw2				! number of segmental waves
READ(10,*)  !wc2					! segmental weighting coefficient
READ(10,*)  !Tmix				  ! period of mixed mode simulation
READ(10,*)  !den					! Liquid's density
READ(10,*)  !nu					  ! Liquid's kinematic viscosity
READ(10,*)  !S_intrinsic	   ! Drug solubility: intrinsic (mu mole/cm^3)
READ(10,*)  !S_bulk     	   ! Drug solubility: at bulk pH of 6.5  (mu mole/cm^3)
READ(10,*)  !diffm           ! Drug's diffusivity (cm2/s)
READ(10,*)  !molarvol        ! Drug's molar volume (cm^3/mole)
READ(10,*)  !den_P           ! Drug's density  (kg/m3)
READ(10,*)  !tau             ! relaxation parameter
READ(10,*)  !Sc              ! Schmidt number
READ(10,*)  !sclrIC          ! initial/maintained scalar distribution (1=BLOB,2=LINE,3=INLET,4=UNIFORM)
READ(10,*)  !iter_Start_phi  ! iteration at which to start particle tracking & scalar calculation  
READ(10,*)  !iter_Freeze_LBM ! iteration at wich steady state (for P & V) has reached so all LBM related functions can be turned OFF  
READ(10,*)  !phiIC           ! maximum scalar concentration
!----- Coefficients for the generalized scalar BC (coeffPhi*phiWall + coeffGrad*dPhiDn_wall = coeffConst). 'n' is normal vector from  wall into fluid.
READ(10,*)  !coeffPhi
READ(10,*)  !coeffGrad
READ(10,*)  !coeffConst
READ(10,*)  !nPers                       ! total number of periods to run
READ(10,*)  !Output_Intervals            ! number of iterations between writing the output files 
READ(10,*)  !Restart_Intervals           ! number of iterations between writing the restart files 
READ(10,*)  !Flag_Buffer                 ! Flag for Buffer Capacity: False-->0mM, TRUE-->10.5mM 
READ(10,*)  !Flag_Couette                ! Flag to run the Couette simulation
READ(10,*)  !Flag_Correcting_Mass        ! Flag for mass correction by bringing back rho to 1.0        
READ(10,*)  !Flag_BounceBack_2nd_Order   ! Flag for 2nd order LBM BC. If False --> 1st order LBM BC 
READ(10,*)  !Flag_ParticleTrack          ! Flag for tracking particles           
READ(10,*)  !Flag_Particle_Init_Sphere   ! Flag to initiate particles in a sphere (TRUE) or in the whole domain (False) 
READ(10,*)  !Flag_Shear_Effects          ! Flag for including shear effects in Sherwood number        
READ(10,*) Flag_Convection_Effects     ! Flag for including convection effects in Sherwood number        
READ(10,*)  !Flag_Confinement_Effects    ! Flag for including confinement effectgs in Sherwood number 
READ(10,*)  !Flag_Rectify_Neg_phi        ! Flag for rectifying negative phi (make it zero) or leave it as is
READ(10,*)  !Flag_Restart                ! Falg for using restart files instead of starting from zero
CLOSE(10)

ALLOCATE(node(0:nx+1,0:ny+1,0:nz+1))
ALLOCATE(u(0:nx+1,0:ny+1,0:nz+1),							&
         v(0:nx+1,0:ny+1,0:nz+1),							&
         w(0:nx+1,0:ny+1,0:nz+1))
ALLOCATE(DUdt_x(0:nx+1,0:ny+1,0:nz+1),				&
         DUdt_y(0:nx+1,0:ny+1,0:nz+1),				&
         DUdt_z(0:nx+1,0:ny+1,0:nz+1))				
ALLOCATE(Laplacian_x(0:nx+1,0:ny+1,0:nz+1),				&
         Laplacian_y(0:nx+1,0:ny+1,0:nz+1),				&
         Laplacian_z(0:nx+1,0:ny+1,0:nz+1))				
ALLOCATE(DLaplacianDt_x(0:nx+1,0:ny+1,0:nz+1),				&
         DLaplacianDt_y(0:nx+1,0:ny+1,0:nz+1),				&
         DLaplacianDt_z(0:nx+1,0:ny+1,0:nz+1))				

ALLOCATE(rho(0:nx+1,0:ny+1,0:nz+1))
ALLOCATE(phi(0:nx+1,0:ny+1,0:nz+1))
ALLOCATE(f(0:NumDistDirs,0:nx+1,0:ny+1,0:nz+1))	

OPEN(11,FILE='Restart-iter.dat')              ! open initial iteration file
READ(11,*) iter0                              ! read and set initial iteration
CLOSE(11)

WRITE(iter_char(1:7),'(I7.7)') iter0
OPEN(12,FILE='Restart-Out-'//iter_char//'-00001.dat')
DO k=0,nz+1
   DO j=0,ny+1
      DO i=0,nx+1
         READ(12,*) node(i,j,k)
         READ(12,*) u(i,j,k)
         READ(12,*) v(i,j,k)
         READ(12,*) w(i,j,k)
         READ(12,*) rho(i,j,k)
         READ(12,*) phi(i,j,k)
         DO m=0,NumDistDirs
            READ(12,*) f(m,i,j,k)
         END DO
         IF (Flag_Convection_Effects) THEN
            READ(12,*)  DUdt_x(i,j,k) 
            READ(12,*)  DUdt_y(i,j,k) 
            READ(12,*)  DUdt_z(i,j,k) 
            READ(12,*)  Laplacian_x(i,j,k)
            READ(12,*)  Laplacian_y(i,j,k)
            READ(12,*)  Laplacian_z(i,j,k)
            READ(12,*)  DLaplacianDt_x(i,j,k)
            READ(12,*)  DLaplacianDt_y(i,j,k)
            READ(12,*)  DLaplacianDt_z(i,j,k)
         ENDIF
     END DO
   END DO
END DO
READ(12,*) Drug_Initial
READ(12,*) Drug_Released
READ(12,*) Drug_Absorbed
READ(12,*) Drug_Remained
delphi_particle = 0.0                               	! Initialize the scalar contirbution from particles to 0.0. Once the particle
CLOSE(12)

! ----- Scattering the single Restart file

N_CPU = NumSubsX * NumSubsY * NumSubsZ 

quotient_X	= CEILING(REAL(nx)/NumSubsX)                ! divide the number of nodes by the number of subdomains (round up)
quotient_Y	= CEILING(REAL(ny)/NumSubsY)                ! divide the number of nodes by the number of subdomains (round up)
quotient_Z	= CEILING(REAL(nz)/NumSubsZ)                ! divide the number of nodes by the number of subdomains (round up)

DO CPU = 0, N_CPU-1
   iMin = MOD(CPU,NumSubsX)*quotient_X + 1               ! starting local i index 
   iMax = iMin + (quotient_X - 1)                      ! ending local i index
   jMin = MOD((CPU/NumSubsX),NumSubsY)*quotient_Y + 1    ! starting local j index
   jMax = jMin + (quotient_Y - 1)							         ! ending local j index
   kMin = (CPU/(NumSubsX*NumSubsY))*quotient_Z + 1	     ! starting local k index 
   kMax = kMin + (quotient_Z - 1)                      ! ending local k index

   !----- Check the bounds  ----------------------------------------------------------
   IF (iMax .GT. nx) THEN
      iMax = nx																! if iMax is greater than nx, correct it
   END IF
   IF (jMax .GT. ny) THEN
      jMax = ny																! if jMax is greater than ny, correct it
   END IF
   IF (kMax .GT. nz) THEN
      kMax = nz																! if kMax is greater than nz, correct it
   END IF
   WRITE(sub(1:5),'(I5.5)') CPU+1	
   OPEN (13,FILE='Restart-Out-'//iter_char//'-'//sub//'.dat')
   DO k=kMin-1,kMax+1
      DO j=jMin-1,jMax+1
         DO i=iMin-1,iMax+1
            WRITE(13,'(I1)') node(i,j,k)
            WRITE(13,'(F8.6)') u(i,j,k)
            WRITE(13,'(F8.6)') v(i,j,k)
            WRITE(13,'(F8.6)') w(i,j,k)
            WRITE(13,'(F8.6)') rho(i,j,k)
            WRITE(13,'(F8.6)') phi(i,j,k)
            DO m=0,NumDistDirs
               WRITE(13,'(F10.8)') f(m,i,j,k)
            END DO
            IF (Flag_Convection_Effects) THEN
               WRITE(13,'(E16.9)')  DUdt_x(i,j,k) 
               WRITE(13,'(E16.9)')  DUdt_y(i,j,k) 
               WRITE(13,'(E16.9)')  DUdt_z(i,j,k) 
               WRITE(13,'(E16.9)')  Laplacian_x(i,j,k)
               WRITE(13,'(E16.9)')  Laplacian_y(i,j,k)
               WRITE(13,'(E16.9)')  Laplacian_z(i,j,k)
               WRITE(13,'(E16.9)')  DLaplacianDt_x(i,j,k)
               WRITE(13,'(E16.9)')  DLaplacianDt_y(i,j,k)
               WRITE(13,'(E16.9)')  DLaplacianDt_z(i,j,k)
            ENDIF
         END DO
      END DO
   END DO
   WRITE(13,*) Drug_Initial
   IF (CPU .EQ. 0) THEN
      WRITE(13,*) Drug_Released
   ELSE
      WRITE(13,*) 0
   END IF   
   WRITE(13,*) Drug_Absorbed
   WRITE(13,*) Drug_Remained
   CLOSE(13)
END DO

end
