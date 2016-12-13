IMPLICIT NONE 

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
CLOSE(10)

ALLOCATE(node(0:nx+1,0:ny+1,0:nz+1))
ALLOCATE(u(0:nx+1,0:ny+1,0:nz+1),							&
         v(0:nx+1,0:ny+1,0:nz+1),							&
         w(0:nx+1,0:ny+1,0:nz+1))
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
   WRITE(sub(1:5),'(I5.5)') CPU	
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
         END DO
      END DO
   END DO
   WRITE(13,*) Drug_Initial
   IF (CPU .EQ. 1) THEN
      WRITE(13,*) Drug_Released
   ELSE
      WRITE(13,*) 0
   END IF   
   WRITE(13,*) Drug_Absorbed
   WRITE(13,*) Drug_Remained
   CLOSE(13)
END DO

end
