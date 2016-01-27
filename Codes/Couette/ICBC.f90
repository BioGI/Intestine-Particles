!==================================================================================================
MODULE ICBC		! Sets Initial and Boundary Conditions
					! Contains setup subroutines (ICs,BounceBack,IniParticles)
!==================================================================================================
USE SetPrecision
USE Setup  
USE MPI

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE ICs	! sets the initial conditions
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE  

! Define local variables
INTEGER(lng) :: i,j,k,m									! index variables
REAL(dbl) :: feq

INTEGER(lng) :: imintemp,imaxtemp	! index variables
REAL(dbl) :: alpha,xmin,xmax,xmid

IF(restart) THEN											! restart from  file 
  
  OPEN(50,FILE='restart.'//sub)						! open correct restart file
  
  DO k=0,nzSub+1_lng
    DO j=0,nySub+1_lng
      DO i=0,nxSub+1_lng

        READ(50,*) node(i,j,k)
        READ(50,*) u(i,j,k)
        READ(50,*) v(i,j,k)
        READ(50,*) w(i,j,k)
        READ(50,*) rho(i,j,k)
        READ(50,*) phi(i,j,k)

        DO m=0,NumDistDirs
          READ(50,*) f(m,i,j,k)
        END DO

      END DO
    END DO
  END DO

  READ(50,*) phiAbsorbed
  READ(50,*) phiAbsorbedS
  READ(50,*) phiAbsorbedV
  READ(50,*) phiInOut
 
  CLOSE(50)

  OPEN(55,FILE='iter0.dat')							! open initial iteration file
  READ(55,*) iter0										! read and set initial iteration
  CLOSE(55)

  iter = iter0-1_lng										! set the initial iteration to the last iteration from the previous run

  IF(randORord .EQ. RANDOM) THEN
    ALLOCATE(rnd(2_lng*numVilli))
    OPEN(1778,FILE='rnd.dat')							! read in the rnd array
    DO i=1,2_lng*numVilli
      READ(1778,*) rnd(i)
    END DO
    CLOSE(1778)
  END IF

ELSE															! clean start
  alpha = 0.4_dbl
  imintemp = -ANINT(alpha*(nx-1_lng))+ANINT(0.5_dbl*(nx+1))
  imaxtemp = +ANINT(alpha*(nx-1_lng))+ANINT(0.5_dbl*(nx+1))
  xmin = 0.5_dbl*(xx(imintemp)+xx(imintemp-1))
  xmax = 0.5_dbl*(xx(imaxtemp)+xx(imaxtemp+1))
  xmid = 0.5_dbl*(xmax+xmin)

  ! Initial conditions on velocity, density, and scalar
  DO k=0,nzSub+1_lng
    DO j=0,nySub+1_lng
      DO i=0,nxSub+1_lng
        u(i,j,k)   = 0.0_dbl!0.0_dbl!0.01_dbl							! x-velocity
        v(i,j,k)   = 0.0_dbl!0.01_dbl							! y-velocity
	IF (node(i,j,k).EQ.FLUID) THEN        
		w(i,j,k)   = (0.5_dbl*s1/vcf)*(xmid-xx(i+iMin-1))/(xmid-xmin)! z-velocity
		!w(i,j,k)   = (s1/vcf)*(xmax-xx(i+iMin-1))/(xmax-xmin)! z-velocity
		!w(i,j,k)   = (s1/vcf)! z-velocity
		!w(i,j,k)   = (s1/vcf)*(alpha*D-xx(i))/(2.0_dbl*alpha*D)!s1/vcf!0.0_dbl!0.0_dbl	! z-velocity
	ELSE
		w(i,j,k) = 0.0_dbl
	ENDIF
        rho(i,j,k) = denL
	!IF (node(i,j,k).EQ.FLUID) THEN
	!	write(*,*) i,j,k,xx(i),w(i,j,k)
	!ENDIF								! density
	! Balaji added
	! distribution functions (set to equilibrium)
	DO m=0,NumDistDirs
	  CALL Equilibrium_LOCAL(m,rho(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),feq)	! distribution functions
	  f(m,i,j,k) = feq
	END DO

      END DO
    END DO
  END DO

  ! Starting iteration
  iter0 = 1_lng
  !iter0 = 0_lng

  ! Initialize scalar values
  phiAbsorbed	= 0.0_dbl								! total amount of scalar absorbed
  phiAbsorbedS	= 0.0_dbl								! total amount of scalar absorbed through the macroscopic surface
  phiAbsorbedV	= 0.0_dbl								! total amount of scalar absorbed through the villi
  phiInOut	= 0.0_dbl								! total amount of scalar leaving the inlet/outlet
  delphi_particle = 0.0_dbl								! Initialize the scalar contirbution from particles to 0.0. Once the particle
											! data is read, we can interpolate to get the delphi_particle. In any event, thi
											! variable is designed to store temporary data. 

END IF

!------------------------------------------------
END SUBROUTINE ICs
!------------------------------------------------

!------------------------------------------------
SUBROUTINE IniParticles
!-----------------------------------------------
IMPLICIT NONE
INTEGER(lng)   :: i, parid,particle_partition,ipartition
REAL(dbl) :: xp,yp,zp,par_radius
TYPE(ParRecord), POINTER	:: CurPar
IF (restart) THEN
	! Read particle number and position along with it's radius,concentration.
	! Interpolate to calculate particle velocities.
	! Interpolate particle concentration to nodes into delphi_particle.

ELSE
	! Linked list approach
	OPEN(60,FILE='particle.dat')
	READ(60,*) np
	num_particles = np

	! Initialize Header Pointer
	
	!ALLOCATE(ParListHead)
	!ParListHead%next => null()!ParListHead
	!ParListHead%prev => null()!ParListHead
	CALL list_init(ParListHead)
	CurPar => ParListHead

    !IF (myid .EQ. master) THEN
	! Recursively allocate all the particle records and build the linked list
	DO i = 1, np
		!ALLOCATE(CurPar%next)
		!CurPar%next%prev => CurPar
		!!CurPar%next%next => ParListHead
		!CurPar%next%next => null()
		!CurPar => CurPar%next
		READ(60,*) parid,xp,yp,zp,par_radius ! read particle.dat file

		! Search the partition this particle belongs to
		DO ipartition = 1_lng,NumSubsTotal 

			IF ((xp.GE.REAL(iMinDomain(ipartition),dbl)-1.0_dbl).AND.&
			(xp.LT.(REAL(iMaxDomain(ipartition),dbl)+0.0_dbl)).AND. &
			(yp.GE.REAL(jMinDomain(ipartition),dbl)-1.0_dbl).AND. &
			(yp.LT.(REAL(jMaxDomain(ipartition),dbl)+0.0_dbl)).AND. &
			(zp.GE.REAL(kMinDomain(ipartition),dbl)-1.0_dbl).AND. &
			(zp.LT.(REAL(kMaxDomain(ipartition),dbl)+0.0_dbl))) THEN

				particle_partition = ipartition
			END IF
		END DO
		! Create a particle element in the linked list only if the particles belongs to this partition
		IF (particle_partition.EQ.mySub) THEN
			CALL list_init(CurPar%next)		
			CurPar%next%prev => CurPar
			CurPar%next%next => null()
			CurPar%next%pardata%parid = parid
			CurPar%next%pardata%xp = xp
			CurPar%next%pardata%yp = yp
			CurPar%next%pardata%zp = zp
			CurPar%next%pardata%up = 0.0_dbl
			CurPar%next%pardata%vp = 0.0_dbl
			CurPar%next%pardata%wp = 0.0_dbl
			CurPar%next%pardata%rp = par_radius!R0!0.005_dbl
			CurPar%next%pardata%xpold = CurPar%next%pardata%xp
			CurPar%next%pardata%ypold = CurPar%next%pardata%yp
			CurPar%next%pardata%zpold = CurPar%next%pardata%zp
			CurPar%next%pardata%upold = CurPar%next%pardata%up
			CurPar%next%pardata%vpold = CurPar%next%pardata%vp
			CurPar%next%pardata%wpold = CurPar%next%pardata%wp
			CurPar%next%pardata%rpold = CurPar%next%pardata%rp
			CurPar%next%pardata%par_conc = Cs_mol!3.14854e-6
			CurPar%next%pardata%gamma_cont = 0.0000_dbl
			CurPar%next%pardata%sh = 1.0000_dbl/(1.0_dbl-CurPar%next%pardata%gamma_cont)
			CurPar%next%pardata%S = 0.0_dbl
			CurPar%next%pardata%Sst = 0.0_dbl
			CurPar%next%pardata%Veff = 0.0_dbl
			CurPar%next%pardata%Nbj = 0.0_dbl
			CurPar%next%pardata%bulk_conc = 0.0000_dbl
			CurPar%next%pardata%delNBbyCV= 0.00000_dbl
			CurPar%next%pardata%cur_part= mySub
			CurPar%next%pardata%new_part= mySub
!			!WRITE(*,*) "Particle Initializing ",i,xp(i),yp(i),zp(i)
!	 		!ss(:,:)=uu(:,:,(nz+1)/2)
!		        !CALL interp(xp(i),yp(i),ss,nx,ny,up(i))
!		        !ss(:,:)=vv(:,:,(nz+1)/2)
!		        !CALL interp(xp(i),yp(i),ss,nx,ny,vp(i))
			! point to next node in the list
			CurPar => CurPar%next
		END IF
	END DO
     !END IF
	
	CLOSE(60)
ENDIF
!------------------------------------------------
END SUBROUTINE IniParticles
!------------------------------------------------

!------------------------------------------------
SUBROUTINE IniParticles_Old
!-----------------------------------------------
IMPLICIT NONE
INTEGER(lng)   :: i, parid
REAL(dbl) :: xp,yp,zp,par_radius
TYPE(ParRecord), POINTER	:: CurPar
IF (restart) THEN
	! Read particle number and position along with it's radius,concentration.
	! Interpolate to calculate particle velocities.
	! Interpolate particle concentration to nodes into delphi_particle.

ELSE
	! Linked list approach
	OPEN(60,FILE='particle.dat')
	READ(60,*) np
	num_particles = np

	! Initialize Header Pointer
	
	!ALLOCATE(ParListHead)
	!ParListHead%next => null()!ParListHead
	!ParListHead%prev => null()!ParListHead
	CALL list_init(ParListHead)
	CurPar => ParListHead

    IF (myid .EQ. master) THEN
	! Recursively allocate all the particle records and build the linked list
	DO i = 1, np
		!ALLOCATE(CurPar%next)
		!CurPar%next%prev => CurPar
		!!CurPar%next%next => ParListHead
		!CurPar%next%next => null()
		!CurPar => CurPar%next
		READ(60,*) parid,xp,yp,zp,par_radius

		CALL list_init(CurPar%next)		
		CurPar%next%prev => CurPar
		CurPar%next%next => null()
		CurPar%next%pardata%parid = parid
		CurPar%next%pardata%xp = xp
		CurPar%next%pardata%yp = yp
		CurPar%next%pardata%zp = zp
		CurPar%next%pardata%up = 0.0_dbl
		CurPar%next%pardata%vp = 0.0_dbl
		CurPar%next%pardata%wp = 0.0_dbl
		CurPar%next%pardata%rp = par_radius!R0!0.005_dbl
		CurPar%next%pardata%xpold = CurPar%next%pardata%xp
		CurPar%next%pardata%ypold = CurPar%next%pardata%yp
		CurPar%next%pardata%zpold = CurPar%next%pardata%zp
		CurPar%next%pardata%upold = CurPar%next%pardata%up
		CurPar%next%pardata%vpold = CurPar%next%pardata%vp
		CurPar%next%pardata%wpold = CurPar%next%pardata%wp
		CurPar%next%pardata%rpold = CurPar%next%pardata%rp
		CurPar%next%pardata%par_conc = Cs_mol!3.14854e-6
		CurPar%next%pardata%gamma_cont = 0.0000_dbl
		CurPar%next%pardata%sh = 1.0000_dbl/(1.0_dbl-CurPar%next%pardata%gamma_cont)
		CurPar%next%pardata%S = 0.0_dbl
		CurPar%next%pardata%Sst = 0.0_dbl
		CurPar%next%pardata%Veff = 0.0_dbl
		CurPar%next%pardata%Nbj = 0.0_dbl
		CurPar%next%pardata%bulk_conc = 0.0000_dbl
		CurPar%next%pardata%delNBbyCV= 0.00000_dbl
		CurPar%next%pardata%cur_part= mySub
		CurPar%next%pardata%new_part= mySub
!		!WRITE(*,*) "Particle Initializing ",i,xp(i),yp(i),zp(i)
! 		!ss(:,:)=uu(:,:,(nz+1)/2)
!	        !CALL interp(xp(i),yp(i),ss,nx,ny,up(i))
!	        !ss(:,:)=vv(:,:,(nz+1)/2)
!	        !CALL interp(xp(i),yp(i),ss,nx,ny,vp(i))
		! point to next node in the list
		CurPar => CurPar%next
		!write(*,*) i
	END DO
     END IF
	
	CLOSE(60)
ENDIF
!------------------------------------------------
END SUBROUTINE IniParticles_Old
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarDistribution		! Sets/Maintains initial distributions of scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,ii,jj		! lattice indices
	!write(*,*) iter
	!pause
! INTRODUCTION OF SCALAR
IF(iter .EQ. phiStart) THEN
	!write(*,*) sclrIC
	!pause
  SELECT CASE(sclrIC) 
      
    CASE(BLOB)							! blob of scalar at the center of the domain
      
      DO k=0,nzSub+1
        DO j=0,nySub+1
          DO i=0,nxSub+1

            !phi(i,j,k) = phiIC*ee**(-((x(i)**2 + y(j)**2 + (z(k)-0.5_dbl*L)**2)/(2.0_dbl*sigma**2)))		! 3D Gaussian Distribution
            IF (((i.GE.16).AND.(i.LE.25)).AND.((j.GE.16).AND.(j.LE.25)).AND.((k.GE.16).AND.(k.LE.25))) THEN
		    !IF (node(i,j,k).EQ.FLUID) THEN
	            	phi(i,j,k) = phiIC!i+j+k!phiIC		! 3D Cube test
		    !ENDIF
	    ENDIF

          END DO
        END DO
      END DO

    CASE(LINE) 						! line of scalar along axis
  
      DO k=0,nzSub+1
        DO j=0,nySub+1
          DO i=0,nxSub+1

            phi(i,j,k) = phiIC*ee**(-((x(i)**2 + y(j)**2)/(2.0_dbl*sigma**2)))									! 2D Gaussian Distribution in x and y

          END DO
        END DO
      END DO

    CASE(INLET) 						! line of scalar at the inlet
  
      DO k=0,nzSub+1
        DO j=0,nySub+1
          DO i=0,nxSub+1

            phi(i,j,k) = phiIC*ee**(-((z(k)**2)/(2.0_dbl*sigma**2)))													! 1D Gaussian Distribution in z

          END DO
        END DO
      END DO

    CASE(UNIFORM)						! uniform initial distribution
	    !write(*,*) "balaji"
	    !pause
      phi(:,:,:) = phiIC			! set the full scalar field to phiIC

    CASE DEFAULT
   
      OPEN(1000,FILE="error.txt")
      WRITE(1000,*) "Error in ScalarIC in Setup.f90: sclrIC is not BLOB(1), LINE(2) or INLET(3)..."
      WRITE(1000,*) "sclrIC=", sclrIC
      CLOSE(1000)
      STOP

  END SELECT


  ! Calculate the intial amount of scalar
  phiTotal = 0.0_dbl
  DO k=1,nzSub
    DO j=1,nySub
      DO i=1,nxSub

        IF(node(i,j,k) .EQ. FLUID) THEN
          phiTotal = phiTotal + phi(i,j,k)
        END IF

      END DO
    END DO
  END DO

ELSE

  ! MAINTAINENCE OF SCALAR
  SELECT CASE(sclrIC) 
      
    CASE(BLOB)							! blob of scalar at the center of the domain

      ! scalar is not maintained
  
    CASE(LINE) 						! line of scalar along axis

      IF((SubID(2) .EQ. 0) .AND. (SubID(4) .EQ. 0)) THEN			! if no neighboring subdomains exist in the 2nd and 4th directions, then they lie at the centerline

!        ! maintain scalar at centerline
!        DO k=0,nzSub+1
!          phi(1,1,k) = phiIC*ee**(-((x(1)**2 + y(1)**2)/(2.0_dbl*sigma**2)))
!        END DO

        DO k=0,nzSub+1
          DO j=0,nySub+1
            DO i=0,nxSub+1

              IF((ABS(x(i)) .LE. 2.51_dbl*xcf) .AND. (ABS(y(j)) .LE. 2.51_dbl*ycf)) THEN				! (5.01 in case it is slightly higher than 5 due to round-off)
                phi(i,j,k) = phiIC																						! 2D Gaussian Distribution in x and y (maintain phi=1 at r=2.5*zcf)
              END IF

            END DO
          END DO
        END DO  

      END IF

    CASE(INLET) 						! constant scalar at the inlet
 
      ! maintain scalar at inlet (1.0)
      IF(kMin .EQ. 1) THEN	
        DO j=1,ny
          DO i=1,nx
            DO k=0,1

	           phi(i,j,k) = phiIC*ee**(-((z(k)**2)/(2.0_dbl*sigma**2)))										! 1D Gaussian Distribution in z for the inlet nodes (maintain peak at z=0)

            END DO
          END DO
        END DO   
      END IF
 
    CASE(UNIFORM)						! uniform initial distribution 

      ! scalar is not maintained

    CASE DEFAULT
   
       OPEN(1000,FILE="error.txt")
       WRITE(1000,*) "Error in ScalarIC in Setup.f90: sclrIC is not BLOB(1), LINE(2) or INLET(3)..."
       WRITE(1000,*) "sclrIC=", sclrIC
       CLOSE(1000)
       STOP

  END SELECT


END IF	

!------------------------------------------------
END SUBROUTINE ScalarDistribution
!------------------------------------------------




!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			! implements the (moving) bounceback boundary conditions (1st order accurate - Ladd)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb											! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk													! radius of current node

!rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
rijk = x(im1)								! height at current location

cosTheta = x(im1)/rijk											! COS(theta)
sinTheta = y(jm1)/rijk											! SIN(theta)

IF (rijk .GE. rOut(k)) THEN
	ub = 0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
	vb = 0.0_dbl!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
	wb = velOut(km1)!vel(km1)!0.0_dbl						! only z-component in this case			
	!ub = -velOut(km1)*sinTheta!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
	!vb = velOut(km1)*cosTheta!0.0!vel(km1)*sinTheta							! y-component of the velocity at i,j,k
	!wb = 0.0_dbl!vel(km1)!0.0_dbl									! no z-component in this case			
	!ub = vel(km1)*cosTheta										! x-component of the velocity at i,j,k
	!vb = vel(km1)*sinTheta										! y-component of the velocity at i,j,k
	!wb = 0.0_dbl											! no z-component in this case
ELSE IF (rijk .LE. rIn(k)) THEN
	ub = 0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
	vb = 0.0_dbl!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
	wb = velIn(km1)!vel(km1)!0.0_dbl						! only z-component in this case	
	!ub = -velIn(km1)*sinTheta!0.0!vel(km1)*cosTheta					! x-component of the velocity at i,j,k
	!vb = velIn(km1)*cosTheta!0.0!vel(km1)*sinTheta					! y-component of the velocity at i,j,k
	!wb = 0.0_dbl!vel(km1)!0.0_dbl							! no z-component in this case			
	!ub = vel(km1)*cosTheta											! x-component of the velocity at i,j,k
	!vb = vel(km1)*sinTheta											! y-component of the velocity at i,j,k
	!wb = 0.0_dbl														! no z-component in this case
END IF				

!fbb = fplus(bb(m),i,j,k) + 6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m))	! bounced back distribution function with added momentum
fbb = fplus(bb(m),i,j,k) + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m))	! bounced back distribution function with added momentum
!fbb = fplus(bb(m),i,j,k) + 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))	! bounced back distribution function with added momentum
fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
!------------------------------------------------
END SUBROUTINE BounceBackL
!------------------------------------------------





!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBack2(m,i,j,k,im1,jm1,km1,fbb)	! implements the (moving) bounceback boundary conditions (2nd order accurate - Lallemand)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2					! index variables
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb											! wall velocity (x-, y-, z- components)
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]
REAL(dbl) :: rijk													! radius of current node

ip1 = i + ex(m)													! i location of 1st neighbor in the m direction
jp1 = j + ey(m)													! j location of 1st neighbor in the m direction
kp1 = k + ez(m)													! k location of 1st neighbor in the m direction

ip2 = i + 2_lng*ex(m)											! i location of 2nd neighbor in the m direction
jp2 = j + 2_lng*ey(m)											! j location of 2nd neighbor in the m direction
kp2 = k + 2_lng*ez(m)											! k location of 2nd neighbor in the m direction


IF((node(ip1,jp1,kp1) .EQ. FLUID) .AND. (node(ip2,jp2,kp2) .EQ. FLUID)) THEN		! continue with 2nd order BB if the two positive neighbors are in the fluid (most cases)

	!rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
	rijk = x(im1)								! height at current location

  cosTheta = x(im1)/rijk										! COS(theta)
  sinTheta = y(jm1)/rijk										! SIN(theta)
  !cosTheta = x(im1)/r(km1)									! COS(theta)
  !sinTheta = y(jm1)/r(km1)									! SIN(theta)

	IF (rijk .GE. rOut(k)) THEN
		ub = 0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
		vb = 0.0_dbl!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
		wb = velOut(km1)!vel(km1)!0.0_dbl						! only z-component in this case			
		!ub = -velOut(km1)*sinTheta!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
		!vb = velOut(km1)*cosTheta!0.0!vel(km1)*sinTheta							! y-component of the velocity at i,j,k
		!wb = 0.0_dbl!vel(km1)!0.0_dbl									! no z-component in this case			
		!ub = vel(km1)*cosTheta										! x-component of the velocity at i,j,k
		!vb = vel(km1)*sinTheta										! y-component of the velocity at i,j,k
		!wb = 0.0_dbl											! no z-component in this case
	ELSE IF (rijk .LE. rIn(k)) THEN
		ub = 0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
		vb = 0.0_dbl!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
		wb = velIn(km1)!vel(km1)!0.0_dbl						! only z-component in this case	
		!ub = -velIn(km1)*sinTheta!0.0!vel(km1)*cosTheta					! x-component of the velocity at i,j,k
		!vb = velIn(km1)*cosTheta!0.0!vel(km1)*sinTheta					! y-component of the velocity at i,j,k
		!wb = 0.0_dbl!vel(km1)!0.0_dbl							! no z-component in this case			
		!ub = vel(km1)*cosTheta											! x-component of the velocity at i,j,k
		!vb = vel(km1)*sinTheta											! y-component of the velocity at i,j,k
		!wb = 0.0_dbl														! no z-component in this case
	END IF	

  CALL qCalc(m,i,j,k,im1,jm1,km1,q)							! calculate q					

  ! bounced back distribution function with added momentum
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
    fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 															&
        + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) 													& 
        - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 													&
        !+ 6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m))
        + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)) ! Set rho = 1.0
        !+ 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m))) ! use actual rho that fluctates

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
    fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q																&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)							&
        !+ (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
        + (6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))	! Set rho = 1.0
        !+ (6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m))))/(q*(2.0_dbl*q + 1.0_dbl)) ! Use actual rho that fluctuateѕ

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
  ELSE
    OPEN(1000,FILE='error-'//sub//'.txt')
    WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
    CLOSE(1000)
    STOP
  END IF

ELSE

  CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)

END IF

!------------------------------------------------
END SUBROUTINE BounceBack2
!------------------------------------------------



!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBack2New(m,i,j,k,im1,jm1,km1,fbb)	! implements the (moving) bounceback boundary conditions (2nd order accurate - Lallemand)
! Implemented by Balaji 10/28/2014 using a method similar to Yanxing
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2					! index variables
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb											! wall velocity (x-, y-, z- components)
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]
REAL(dbl) :: rijk													! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt				! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 
INTEGER(lng) :: it			! loop index variables

ip1 = i + ex(m)													! i location of 1st neighbor in the m direction
jp1 = j + ey(m)													! j location of 1st neighbor in the m direction
kp1 = k + ez(m)													! k location of 1st neighbor in the m direction

ip2 = i + 2_lng*ex(m)											! i location of 2nd neighbor in the m direction
jp2 = j + 2_lng*ey(m)											! j location of 2nd neighbor in the m direction
kp2 = k + 2_lng*ez(m)											! k location of 2nd neighbor in the m direction


IF((node(ip1,jp1,kp1) .EQ. FLUID) .AND. (node(ip2,jp2,kp2) .EQ. FLUID)) THEN		! continue with 2nd order BB if the two positive neighbors are in the fluid (most cases)

!*****************************************************************************
	!rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
	rijk = x(im1)								! height at current location
		 ! Initial fluid node guess
                 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                
		 ! Initial solid node guess
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
                 
	 IF (k.NE.km1) THEN
                 DO it=1,qitermax
		   ! guess of boundary location 
                   xt=(x1+x2)/2.0_dbl
                   yt=(y1+y2)/2.0_dbl
                   zt=(z1+z2)/2.0_dbl

      		   !rt = SQRT(xt*xt + yt*yt)
      		   rt = xt
		   !Write(*,*) 'test'
                   IF (rijk .GE. rOut(k)) THEN
		   	!ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
		   	ht = ((zt-z(k))*rOut(km1)+(z(km1)-zt)*rOut(k))/(z(km1)-z(k))
		   	!ht = (r(km1)+r(k))/2.0_dbl
		   ELSE
		   	!ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
		   	ht = ((zt-z(k))*rIn(km1)+(z(km1)-zt)*rIn(k))/(z(km1)-z(k))
		   	!ht = (r(km1)+r(k))/2.0_dbl
		   END IF

                   IF(rt.GT.ht) then
                     x2=xt
                     y2=yt
                     z2=zt
                   ELSE
                     x1=xt
                     y1=yt
                     z1=zt
                   END IF
				   
                 END DO
		 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                 
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
 
                 q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
		 !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
	 ELSE
		 DO it=1,qitermax
		   ! guess of boundary location 
                   xt=(x1+x2)/2.0_dbl
                   yt=(y1+y2)/2.0_dbl
                   zt=(z1+z2)/2.0_dbl

      		   !rt = SQRT(xt*xt + yt*yt)
      		   rt = xt
		   !Write(*,*) 'test'
		   IF (rijk .GE. rOut(k)) THEN
			   !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
			   !ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
			   ht = (rOut(km1)+rOut(k))/2.0_dbl
		   ELSE
			   !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
			   !ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
			   ht = (rIn(km1)+rIn(k))/2.0_dbl
		   END IF

                   IF(rt.GT.ht) then
                     x2=xt
                     y2=yt
                     z2=zt
                   ELSE
                     x1=xt
                     y1=yt
                     z1=zt
                   END IF
                 END DO

		 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                 
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
 
                 q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
		 !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
	 ENDIF
		 cosTheta=xt/rt
		 sinTheta=yt/rt
	 IF (k.NE.km1) THEN
		   IF (rijk .GE. rOut(k)) THEN
			 !vt = (ABS(zt-z(k))*vel(km1)+ABS(z(km1)-zt)*vel(k))/ABS(z(km1)-z(k))
			 vt = ((zt-z(k))*velOut(km1)+(z(km1)-zt)*velOut(k))/(z(km1)-z(k))
		   ELSE
			 !vt = (ABS(zt-z(k))*vel(km1)+ABS(z(km1)-zt)*vel(k))/ABS(z(km1)-z(k))
			 vt = ((zt-z(k))*velIn(km1)+(z(km1)-zt)*velIn(k))/(z(km1)-z(k))
		   END IF
	 ELSE
		   IF (rijk .GE. rOut(k)) THEN
			 vt = (velOut(k)+velOut(km1))*0.5_dbl
		   ELSE
			 vt = (velIn(k)+velIn(km1))*0.5_dbl
		   END IF
	 ENDIF
		ub = 0.0_dbl!0.0!vel(km1)*cosTheta								! x-component of the velocity at i,j,k
		vb = 0.0_dbl!0.0!vel(km1)*sinTheta								! y-component of the velocity at i,j,k
		wb = vt!vel(km1)!0.0_dbl									! only z-component in this case)
		!ub = -vt*sinTheta!0.0!vel(km1)*cosTheta											! x-component of the velocity at i,j,k
		!vb = vt*cosTheta!0.0!vel(km1)*sinTheta											! y-component of the velocity at i,j,k
		!wb = 0.0_dbl!vel(km1)!0.0_dbl											! no z-component in this case)
		!ub = vt*cosTheta										! x-component of the velocity at i,j,k
		!vb = vt*sinTheta										! y-component of the velocity at i,j,k
		!wb = 0.0_dbl											! no z-component in this case)
		!write(*,*) 'ht',ht,rt
		!write(*,*) 'q-yanxing',q,i,j,k,im1,jm1,km1
	! make sure 0<q<1
        IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
          OPEN(1000,FILE="error.txt")
          WRITE(1000,*) "q=",q
          WRITE(1000,*) "m=",m
          WRITE(1000,*) "i=",i,"j=",j,"k=",k
          WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
          CLOSE(1000)
          STOP
        END IF	
	!q = max(q, 0.25_dbl)
!*****************************************************************************
!*****************************************************************************
!! Original method used by Gino
!  rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
!
!  cosTheta = x(im1)/rijk										! COS(theta)
!  sinTheta = y(jm1)/rijk										! SIN(theta)
!  !cosTheta = x(im1)/r(km1)									! COS(theta)
!  !sinTheta = y(jm1)/r(km1)									! SIN(theta)
!
!  ub = vel(km1)*cosTheta										! x-component of the velocity at i,j,k
!  vb = vel(km1)*sinTheta										! y-component of the velocity at i,j,k
!  wb = 0.0_dbl														! no z-component in this case)
!
!  CALL qCalc(m,i,j,k,im1,jm1,km1,q)							! calculate q					
!*****************************************************************************
		 !write(*,*) 'q-Gino',q,i,j,k,im1,jm1,km1


  ! bounced back distribution function with added momentum
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
    fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 															&
        + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) 													& 
        - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 													&
        !+ 6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)) ! use actual rho
        + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)) ! use rho = 1.0
        !+ 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))

!    fbb = (2.0_dbl*q)*fplus(bb(m),i,j,k)+ &
!	(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip1,jp1,kp1) &
!        + 6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m))
!        !+ 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
    fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q																&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)							&
        !+ (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))	! Use actual rho
        + (6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))		! Use rho = 1.0
        !+ (6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m))))/(q*(2.0_dbl*q + 1.0_dbl))

!    fbb = fplus(bb(m),i,j,k)/(2.0_dbl*q)+ &
!		(2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k)/(2.0_dbl*q) &
!        + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m))) &
!	/(2.0_dbl*q)
!        !+ (6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m))))/(2.0_dbl*q)

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))

  ELSE
    OPEN(1000,FILE='error-'//sub//'.txt')
    WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
    CLOSE(1000)
    STOP
  END IF
  !IF (k.EQ.42) THEN
  !	write(9,*) i,j,k,q
  !ENDIF

ELSE

  CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)

END IF

!------------------------------------------------
END SUBROUTINE BounceBack2New
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalc(m,i,j,k,im1,jm1,km1,q)			! calculates q (boundary distance ratio) using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
REAL(dbl), INTENT(OUT) :: q							! distance ratio
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Bx,By,Bz									! solid node
REAL(dbl) :: AB,AP										! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: r1,r2,z1,z2,slope,intercept			! radius and z location at k and km1, slope of line connecting those two points, z-intercept of the r-equation
REAL(dbl) :: slope2,term1,term2						! terms used in calculation

! RAY
! point A (current node)
Ax = x(i)
Ay = y(j)
Az = z(k)

! point B (solid node)
Bx = x(im1)
By = y(jm1)
Bz = z(km1)

! distance from A to B
AB = SQRT((Bx - Ax)*(Bx - Ax) + (By - Ay)*(By - Ay) + (Bz - Az)*(Bz - Az))

! unit vector (d) from point A to point B
dx = (x(im1)-x(i))/AB									! i direction
dy = (y(jm1)-y(j))/AB									! j direction
dz = (z(km1)-z(k))/AB									! k direction

! SURFACE
IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	r1 = rOut(k)													! radius at k (distance from CL)
	r2 = rOut(km1)													! radius at km1 (distance from CL)
ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
	r1 = rIn(k)													! radius at k (distance from CL)
	r2 = rIn(km1)													! radius at km1 (distance from CL)
END IF
z1 = z(k)													! z-coordinate at k
z2 = z(km1)													! z-coordinate at km1

IF(k .NE. km1) THEN
  slope = (r2-r1)/(z2-z1)								! approximate the surface as a conincal shell (linear between k values)
ELSE
  slope = 0.0_dbl
END IF

intercept = r1 - slope*z1								! z-intercept of the linearly approximated r-equation

! terms used in calculation
slope2 = slope*slope															! slope^2
term1 = Ax*dx + Ay*dy - Az*dz*slope2 - intercept*dz*slope		! reoccuring term
term2 = dx*dx + dy*dy - dz*dz*slope*slope								! reoccuring term

! calculate the distance from the current node (point A) to the wall (point P)
AP = (1.0_dbl/(2.0_dbl*term2)) * &
     (-2.0_dbl*term1					&
   + SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))

q = AP/AB													! distance ratio
q = max(q,0.001_dbl)

! balaji added
!q=0.5

! make sure 0<q<1
IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) "q=",q
  WRITE(1000,*) "m=",m
  WRITE(1000,*) "i=",i,"j=",j,"k=",k
  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
  WRITE(1000,*) "Ax=",Ax,"Ay=",Ay,"Az=",Az
  WRITE(1000,*) "Bx=",Bx,"By=",By,"Bz=",Bz
  WRITE(1000,*) "dx=",dx,"dy=",dy,"dz=",dz
  WRITE(1000,*) "r1=",r1,"r2=",r2
  WRITE(1000,*) "z1=",z1,"z2=",z2
  WRITE(1000,*) "slope=",slope
  WRITE(1000,*) "term1=",term1,"term2=",term2
  WRITE(1000,*) "intercept=",intercept
  WRITE(1000,*) "AB=",AB,"AP=",AP
  CLOSE(1000)
  STOP
END IF																																

!------------------------------------------------
END SUBROUTINE qCalc
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBackVL(m,i,j,k,im1,jm1,km1,vNum,fbb)	! implements the (moving) bounceback boundary conditions (1st order accurate - Ladd)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
INTEGER(lng), INTENT(IN) :: vNum								! number of the current villus
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
REAL(dbl) :: Cx,Cy,Cz,Vx,Vy,Vz								! vector between villous base and current node, vector between villous base and villous tip
REAL(dbl) :: uV,vV,wV											! villus wall velocity at actual coordinate system
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]

! find C and V
CALL qCalcV(m,i,j,k,im1,jm1,km1,vNum,q,Cx,Cy,Cz,Vx,Vy,Vz,1)

! find the influence of villous velocity on the current point
CALL VilliVelocity(vNum,Cx,Cy,Cz,uV,vV,wV)

! calculate bounced back distribution function with added momentum
fbb = fplus(bb(m),i,j,k) + 6.0_dbl*wt(m)*rho(i,j,k)*(uV*ex(m) + vV*ey(m) + wV*ez(m))	

!------------------------------------------------
END SUBROUTINE BounceBackVL
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBackV2(m,i,j,k,im1,jm1,km1,vNum,fbb)	! implements the (moving) bounceback boundary conditions (2nd order accurate - Lallemand)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
INTEGER(lng), INTENT(IN) :: vNum								! number of the current villus
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2					! index variables
REAL(dbl) :: Cx,Cy,Cz,Vx,Vy,Vz								! vector between villous base and current node, vector between villous base and villous tip
REAL(dbl) :: uV,vV,wV											! villus wall velocity at actual coordinate system
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]

! neighboring nodes
ip1 = i + ex(m)													! i location of 1st neighbor in the m direction
jp1 = j + ey(m)													! j location of 1st neighbor in the m direction
kp1 = k + ez(m)													! k location of 1st neighbor in the m direction

ip2 = i + 2_lng*ex(m)											! i location of 2nd neighbor in the m direction
jp2 = j + 2_lng*ey(m)											! j location of 2nd neighbor in the m direction
kp2 = k + 2_lng*ez(m)											! k location of 2nd neighbor in the m direction

IF((node(ip1,jp1,kp1) .EQ. FLUID) .AND. (node(ip2,jp2,kp2) .EQ. FLUID)) THEN		! continue with 2nd order BB if the two positive neighbors are in the fluid (most cases)

  ! find q (and C and V)
  CALL qCalcV(m,i,j,k,im1,jm1,km1,vNum,q,Cx,Cy,Cz,Vx,Vy,Vz,1)

  ! find the influence of villous velocity on the current point
  CALL VilliVelocity(vNum,Cx,Cy,Cz,uV,vV,wV)

  ! bounced back distribution function with added momentum
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
    fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 															&
        + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) 													& 
        - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 													&
        + 6.0_dbl*wt(m)*rho(i,j,k)*(uV*ex(m) + vV*ey(m) + wV*ez(m))
  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
    fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q																&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)							&
        + (6.0_dbl*wt(m)*rho(i,j,k)*(uV*ex(m) + vV*ey(m) + wV*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
  ELSE
    OPEN(1000,FILE='error-'//sub//'.txt')
    WRITE(1000,*) "Error in BounceBackV2() in ICBC.f90 (line 516): q is not (0<=q<=1)...? Aborting."
    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
    CLOSE(1000)
    STOP
  END IF

ELSE			! must be a case close to two boundaries (villi and surface) - use 1st order BB

  CALL BounceBackVL(m,i,j,k,im1,jm1,km1,vNum,fbb)

END IF

!------------------------------------------------
END SUBROUTINE BounceBackV2
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE VilliVelocity(vNum,Cx,Cy,Cz,uV,vV,wV)					! calculates the influence of the velocity of the villi on point (i,j,k)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: vNum											! number of the current villus
REAL(dbl), INTENT(IN) :: Cx,Cy,Cz										! vector between villous base and current node
REAL(dbl), INTENT(OUT) :: uV,vV,wV										! villus velocity
REAL(dbl) :: ubx,uby,rijk													! wall velocity, radius at current villous base location
REAL(dbl) :: omegaX, omegaY, omegaZ, omegaT							! angular velocity components - x,y,z, and azimuthal directions
REAL(dbl) :: alpha															! angle made between the x-axis and the base of the villus													

! find the angular velocity components
omegaT	= (villiLoc(vNum,5) - villiLoc(vNum,9))/tcf				! angular velocity in the aziumthal direction
alpha		= ATAN(villiLoc(vNum,2)/villiLoc(vNum,1))					! angle made between the x-axis and the base of the villus
omegaX	= -omegaT*SIN(alpha)												! angular velocity in the x-direction
omegaY	= omegaT*COS(alpha)												! angular velocity in the y-direction
omegaZ	= (villiLoc(vNum,4) - villiLoc(vNum,10))/tcf				! angular velocity in the z-direction

! find the villous velocity components from the cross product of omega and C
uV = (omegaY*Cz - omegaZ*Cy)/vcf											! x-component
vV = (omegaZ*Cx - omegaX*Cz)/vcf											! y-component
wV = (omegaX*Cy - omegaY*Cx)/vcf											! z-component

! add the components from the wall to the u and v velocity
!rijk = SQRT(villiLoc(vNum,1)**2 + villiLoc(vNum,2)**2)						! radius at the villous base
!ubx = (velDom(NINT(villiLoc(vNum,3)/zcf))/vcf)*(villiLoc(vNum,1)/rijk)	! x-component of the wall velocity
!uby = (velDom(NINT(villiLoc(vNum,3)/zcf))/vcf)*(villiLoc(vNum,2)/rijk)	! y-component of the wall velocity

!uV = uV + (velDom(NINT(villiLoc(vNum,3)/zcf))/vcf)*COS(alpha)	! active villous velocity component + wall velocity component (x-direction)
!vV = vV + (velDom(NINT(villiLoc(vNum,3)/zcf))/vcf)*SIN(alpha)	! active villous velocity component + wall velocity component (y-direction)
uV = uV + (velDomOut(NINT(villiLoc(vNum,3)/zcf))/vcf)*COS(alpha)	! active villous velocity component + wall velocity component (x-direction)
vV = vV + (velDomOut(NINT(villiLoc(vNum,3)/zcf))/vcf)*SIN(alpha)	! active villous velocity component + wall velocity component (y-direction)

!uV = uV + ubx																	! active villous velocity component + wall velocity component (x-direction)
!vV = vV + uby																	! active villous velocity component + wall velocity component (y-direction)

!------------------------------------------------
END SUBROUTINE VilliVelocity
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalcV(m,i,j,k,im1,jm1,km1,vNum,q,Cx,Cy,Cz,Vx,Vy,Vz,scalarORfluid)	! calculates q (boundary distance ratio) using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
INTEGER(lng), INTENT(IN) :: vNum						! number of the current villus
REAL(dbl), INTENT(OUT) :: q,Cx,Cy,Cz,Vx,Vy,Vz	! distance ratio
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Bx,By,Bz									! solid node
REAL(dbl) :: ABx,ABy,ABz								! vector between nodes
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: AB,AP										! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: term1,term2,term3						! terms used in the solution of AP
REAL(dbl) :: t1,t2										! two solutions for AP (AP is the smaller one)
REAL(dbl) :: dotAV,AMag,cosThetaAV,sinThetaAV	! dot product of A and V, magnitudes of A and V, cosine and sine of the angle between A and V
REAL(dbl) :: dotBV,BMag,cosThetaBV					! dot product of B and V, magnitudes of B, angle between B and V
REAL(dbl) :: VMag											! magnitudes of V
INTEGER(lng) :: cylORsphere							! flag to tell whether the solid node is in the cylinder or the sphere (1 for cyl, 2 for sphere)
INTEGER(lng) :: scalarORfluid							! flag to tell whether qCalcV was called for scalar or fluid BC
INTEGER(lng) :: mpierr									! MPI standard error variable

! find the vector between the villous base and the current node (shift coordinate system)
Ax = (x(i) - villiLoc(vNum,1))						! x-coordinate
Ay = (y(j) - villiLoc(vNum,2))						! y-coordinate
Az = (z(k) - villiLoc(vNum,3))						! z-coordinate

! find the vector between the villous base and the solid node (shift coordinate system)
Bx = (x(im1) - villiLoc(vNum,1))						! x-coordinate
By = (y(jm1) - villiLoc(vNum,2))						! y-coordinate
Bz = (z(km1) - villiLoc(vNum,3))						! z-coordinate

! find the vector between the current node and the solid node
ABx = x(im1) - x(i)										! x-coordinate 
ABy = y(jm1) - y(j)										! y-coordinate 				
ABz = z(km1) - z(k)										! z-coordinate 
AB = SQRT(ABx*ABx + ABy*ABy + ABz*ABz)				! magnitude of AB

! find the vector between the villous base and the villous tip
Vx = (villiLoc(vNum,6)-villiLoc(vNum,1))			! x-coordinate
Vy = (villiLoc(vNum,7)-villiLoc(vNum,2))			! y-coordinate
Vz = (villiLoc(vNum,8)-villiLoc(vNum,3))			! z-coordinate

! find the (unit) vector between the current node and solid node
dx = ABx/AB													! x-coordinate
dy = ABy/AB													! y-coordinate
dz = ABz/AB													! z-coordinate

! determine if the point is above or below the top of the villous cylinder (touching a solid node within the hemi-sphere or the within the cylinder)
!		and calculate the distance from the current node to the wall accordingly
dotAV	= Ax*Vx + Ay*Vy + Az*Vz							! dot product of A and V
dotBV	= Bx*Vx + By*Vy + Bz*Vz							! dot product of A and V
AMag	= SQRT(Ax*Ax + Ay*Ay + Az*Az)					! magnitude of A
BMag	= SQRT(Bx*Bx + By*By + Bz*Bz)					! magnitude of B
VMag	= SQRT(Vx*Vx + Vy*Vy + Vz*Vz)					! magnitude of V
cosThetaAV = dotAV/(AMag*VMag)						! cosine of angle between A and V
sinThetaAV = SQRT(1-cosThetaAV*cosThetaAV)		! sine of the angle between A and V
cosThetaBV = dotBV/(BMag*VMag)						! cosine of angle between B and V
IF((AMag*cosThetaAV .GE. VMag) .AND. ((BMag*cosThetaBV .GE. VMag) .OR. (AMag*sinThetaAV .LE. Rv))) THEN	! the solid node is in the sphere
  ! Calculate the distance from the current node (point A) to the wall (point P) (mathematica solution)  						
  term1 = (Ax - Vx)*dx + (Ay - Vy)*dy + (Az - Vz)*dz
  term2 = term1*term1
  term3 = (Ax - Vx)*(Ax - Vx) + (Ay - Vy)*(Ay - Vy) + (Az - Vz)*(Az - Vz) - Rv*Rv
  t1 = -term1 + SQRT(term2 - term3)
  t2 = -term1 - SQRT(term2 - term3)
  AP = MIN(t1,t2)											! distance from the current node (point A) to the villus (point P)							
  q = AP/AB													! distance ratio  
  cylORsphere = 2_lng
ELSE															! solid node is in the cylinder
  ! Calculate the distance from the current node (point A) to the wall (point P) (mathematica solution)
  term1 = Ay*dy*Vx*Vx + Az*dz*Vx*Vx - Ay*dx*Vx*Vy - Ax*dy*Vx*Vy + Ax*dx*Vy*Vy + Az*dz*Vy*Vy - Az*dx*Vx*Vz - Ax*dz*Vx*Vz	&
        - Az*dy*Vy*Vz - Ay*dz*Vy*Vz + Ax*dx*Vz*Vz + Ay*dy*Vz*Vz
  term2 = 0.5_dbl*SQRT((2.0_dbl*Ay*dy*Vx*Vx + 2.0_dbl*Az*dz*Vx*Vx - 2.0_dbl*Ay*dx*Vx*Vy - 2.0_dbl*Ax*dy*Vx*Vy + 2.0_dbl*Ax*dx*Vy*Vy &
        + 2.0_dbl*Az*dz*Vy*Vy - 2.0_dbl*Az*dx*Vx*Vz - 2.0_dbl*Ax*dz*Vx*Vz - 2.0_dbl*Az*dy*Vy*Vz - 2.0_dbl*Ay*dz*Vy*Vz &
        + 2.0_dbl*Ax*dx*Vz*Vz + 2.0_dbl*Ay*dy*Vz*Vz)**2 - 4.0_dbl*(dy*dy*Vx*Vx + dz*dz*Vx*Vx - 2.0_dbl*dx*dy*Vx*Vy + dx*dx*Vy*Vy &
        + dz*dz*Vy*Vy - 2.0_dbl*dx*dz*Vx*Vz - 2.0_dbl*dy*dz*Vy*Vz + dx*dx*Vz*Vz + dy*dy*Vz*Vz)*(Ay*Ay*Vx*Vx + Az*Az*Vx*Vx &
        - 2.0_dbl*Ax*Ay*Vx*Vy + Ax*Ax*Vy*Vy + Az*Az*Vy*Vy - 2.0_dbl*Ax*Az*Vx*Vz - 2.0_dbl*Ay*Az*Vy*Vz + Ax*Ax*Vz*Vz + Ay*Ay*Vz*Vz &
        - Rv*Rv*Vx*Vx - Rv*Rv*Vy*Vy - Rv*Rv*Vz*Vz))
  term3 = dy*dy*Vx*Vx + dz*dz*Vx*Vx - 2.0_dbl*dx*dy*Vx*Vy + dx*dx*Vy*Vy + dz*dz*Vy*Vy - 2.0_dbl*dx*dz*Vx*Vz - 2.0_dbl*dy*dz*Vy*Vz 	&
        + dx*dx*Vz*Vz + dy*dy*Vz*Vz
  t1 = -(term1 + term2)/term3
  t2 = -(term1 - term2)/term3
  AP = MIN(t1,t2)											! distance from the current node (point A) to the villus (point P)							
  q = AP/AB													! distance ratio
  cylORsphere = 1_lng
END IF

! fix slight discrepencies in calculation of q and node tracking
IF((q .LT. 0_lng) .AND. (q .GT. -0.3_dbl)) THEN
  q = 0.00000001_dbl
ELSE IF ((q .GT. 1_lng) .AND. (q .LT. 1.3_dbl)) THEN
  q = 0.99999999_dbl
!ELSE IF(((q .LT. -0.3_dbl) .AND. (q .GT. -2.0_dbl)) .OR. ((q .GT. 1.3_dbl) .AND. (q .LT. 2.0_dbl))) THEN
ELSE
  q = 0.5_dbl
END IF 

!! make sure 0<q<1
!IF((q .LT. 0_lng) .OR. (q .GT. 1_lng)) THEN 
!
!  OPEN(1000,FILE='error.'//sub//'.txt')
!  WRITE(1000,'(A50)') "error in ICBC.f90 at Line 753: q is out of range"
!  WRITE(1000,*) "iter",iter
!  WRITE(1000,*) "vNum",vNum
!  WRITE(1000,*) "m=",m
!  WRITE(1000,*) "i=",i,"j=",j,"k=",k
!  WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
!  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
!  WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
!  WRITE(1000,*) "node(i,j,k)=",node(i,j,k)
!  WRITE(1000,*) "node(im1,jm1,km1)=",node(im1,jm1,km1)
!  WRITE(1000,*) "Ax=",Ax,"Ay=",Ay,"Az=",Az
!  WRITE(1000,*) "Bx=",Bx,"By=",By,"Bz=",Bz
!  WRITE(1000,*) "ABx=",ABx,"ABy=",ABy,"ABz=",ABz
!  WRITE(1000,*) "Vx=",Vx,"Vy=",Vy,"Vz=",Vz
!  WRITE(1000,*) "villiLoc(vNum,1)=",villiLoc(vNum,1),"villiLoc(vNum,2)=",villiLoc(vNum,2),"villiLoc(vNum,3)=",villiLoc(vNum,3)  
!  WRITE(1000,*) "villiLoc(vNum,6)=",villiLoc(vNum,6),"villiLoc(vNum,7)=",villiLoc(vNum,7),"villiLoc(vNum,8)=",villiLoc(vNum,8)
!  WRITE(1000,*) "dx=",dx,"dy=",dy,"dz=",dz
!  WRITE(1000,*) "term1=",term1,"term2=",term2,"term3=",term3
!  WRITE(1000,*) "t1=",t1,"t2=",t2
!  WRITE(1000,*) "AB=",AB,"AP=",AP
!  WRITE(1000,*) "q=",q
!  WRITE(1000,*) "cylORsphere=",cylORsphere
!  WRITE(1000,*) "scalarORfluid=",scalarORfluid
!  CLOSE(1000)
!
!  OPEN(1001,FILE='ABV.'//sub//'.dat')
!  WRITE(1001,'(A60)') 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
!  WRITE(1001,'(3E15.5,6I6)') x(i), y(j), z(k), 0, 0, 0, 0, 0, -vNum
!  WRITE(1001,'(3E15.5,6I6)') x(im1), y(jm1), z(km1), 0, 0, 0, 0, 0, -vNum
!  WRITE(1001,'(3E15.5,6I6)') villiLoc(vNum,1), villiLoc(vNum,2), villiLoc(vNum,3), 0, 0, 0, 0, 0, -vNum
!  WRITE(1001,'(3E15.5,6I6)') villiLoc(vNum,6), villiLoc(vNum,7), villiLoc(vNum,8), 0, 0, 0, 0, 0, -vNum
!  CLOSE(1001)
!
!  CALL PrintFieldsTEST										! output the velocity, density, and scalar fields [MODULE: Output]
!
!  STOP
!
!END IF						

! find the vector from the base of the villus to the point P		
Cx = Ax + AP*dx											! x-coordinate
Cy = Ay + AP*dy											! y-coordinate
Cz = Az + AP*dz											! z-coordinate								

!------------------------------------------------
END SUBROUTINE qCalcV
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE PrintFieldsTEST	! print velocity, density, and scalar to output files
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng)	:: i,j,k,ii,jj,kk,n		! index variables (local and global)
CHARACTER(7)	:: iter_char				! iteration stored as a character

  ! scale the iteration by 1/10 such that the numbers used in the output file aren't too large
  WRITE(iter_char(1:7),'(I7.7)') iter

  ! store the current iteration in "filenum"
  filenum(fileCount) = iter
  fileCount = fileCount + 1_lng

  ! open the proper output file
  OPEN(60,FILE='out-TEST-'//iter_char//'-'//sub//'.dat')
  WRITE(60,*) 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
  WRITE(60,'(A10,E15.5,A5,I4,A5,I4,A5,I4,A8)') 'ZONE T="',iter/(nt/nPers),'" I=',nxSub,' J=',nySub,' K=',nzSub,'F=POINT'

  DO k=1,nzSub
    DO j=1,nySub
      DO i=1,nxSub

         ! convert local i,j,k, to global ii,jj,kk
         ii = ((iMin - 1_lng) + i)
         jj = ((jMin - 1_lng) + j)
         kk = ((kMin - 1_lng) + k)

         WRITE(60,'(8E15.5,I6)') xx(ii), yy(jj), zz(kk), u(i,j,k)*vcf, v(i,j,k)*vcf, w(i,j,k)*vcf, (rho(i,j,k)-denL)*dcf*pcf,	&
                                 phi(i,j,k), node(i,j,k)

      END DO
    END DO
  END DO

  CLOSE(60)

!  ! print villi locations
!  IF(myid .EQ. master) THEN
!
!    OPEN(607,FILE='villi-'//iter_char//'.dat')
!    OPEN(608,FILE='villi2-'//iter_char//'.dat')
!    WRITE(607,'(A60)') 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
!    WRITE(608,'(A60)') 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
!
!    DO n=1,numVilli
!
!      WRITE(607,'(3E15.5,6I4)') villiLoc(n,1), villiLoc(n,2), villiLoc(n,3), 0, 0, 0, 0, 0, 0
!      WRITE(608,'(3E15.5,6I4)') villiLoc(n,6), villiLoc(n,7), villiLoc(n,8), 0, 0, 0, 0, 0, 0
!
!    END DO
!
!    CLOSE(607)
!    CLOSE(608)
!
!  END IF

!------------------------------------------------
END SUBROUTINE PrintFieldsTEST
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SymmetryBC							! implements symmetry boundary conditions
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) 	:: i,j,k,m,ii,jj,iComm		! index variables
INTEGER(lng)	:: mpierr						! MPI standard error variable

! Loop through the subdomain faces
! -YZ Face
iComm = 2
IF(SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   
  i = YZ_RecvIndex(iComm)	   				! i location of phantom nodes
  ii = YZ_SendIndex(iComm) + 1_lng			! i location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng
    DO j=0,nySub+1_lng
      DO m=1,NumFs_Face

        f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),ii,j,k)			! symmetry BC for 'f'

      END DO

      rho(i,j,k) = rho(ii,j,k)    			! symmetry BC for density 
      phi(i,j,k) = phi(ii,j,k)    			! symmetry BC for scalar

    END DO
  END DO

END IF

! -ZX Face
iComm = 4
IF(SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   
  j = ZX_RecvIndex(iComm)						! j location of phantom nodes
  jj = ZX_SendIndex(iComm) + 1_lng			! j location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng
    DO i=0,nxSub+1_lng
      DO m=1,NumFs_Face

        f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),i,jj,k)        ! symmetry BC

      END DO
 
      rho(i,j,k) = rho(i,jj,k)    			! symmetry BC for density 
      phi(i,j,k) = phi(i,jj,k)    			! symmetry BC for scalar

    END DO
  END DO

END IF

! Z Axis
iComm = 8
IF((SubID(2) .EQ. 0) .AND. (SubID(4) .EQ. 0) .AND. (SubID(iComm) .EQ. 0)) THEN									! if no neighbor in the iComm communication direction exists, iplement symmetry BC

  i = Z_RecvIndex(iComm,1)						! i location of phantom nodes	
  ii = Z_SendIndex(iComm,1) + 1_lng			! i location of 1 row in from the boundary nodes

  j = Z_RecvIndex(iComm,2)						! j location of phantom nodes	
  jj = Z_SendIndex(iComm,2) + 1_lng			! j location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng

    DO m=1,NumFs_Side

        f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),ii,jj,k)        ! symmetry BC

    END DO

    rho(i,j,k) = rho(ii,jj,k)  		  		! symmetry BC for density 
    phi(i,j,k) = phi(ii,jj,k)    			! symmetry BC for scalar

  END DO

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)	! synchronize all processing units

!------------------------------------------------
END SUBROUTINE SymmetryBC
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SymmetryBC_NODE						! implements symmetry boundary conditions
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) 	:: i,j,k,m,ii,jj,iComm		! index variables
INTEGER(lng)	:: mpierr						! MPI standard error variable

! Loop through the subdomain faces
! -YZ Face
iComm = 2
IF(SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   
  i = YZ_RecvIndex(iComm)	   				! i location of phantom nodes
  ii = YZ_SendIndex(iComm) + 1_lng			! i location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng
    DO j=0,nySub+1_lng

      node(i,j,k) = node(ii,j,k)    		! symmetry BC for node flag 

    END DO
  END DO

END IF

! -ZX Face
iComm = 4
IF(SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   
  j = ZX_RecvIndex(iComm)						! j location of phantom nodes
  jj = ZX_SendIndex(iComm) + 1_lng			! j location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng
    DO i=0,nxSub+1_lng
 
      node(i,j,k) = node(i,jj,k)    		! symmetry BC for node flag 

    END DO
  END DO

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)	! synchronize all processing units

!------------------------------------------------
END SUBROUTINE SymmetryBC_NODE
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarBC(m,i,j,k,im1,jm1,km1,phiBC)								! implements the scalar BCs 
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1								! index variables
REAL(dbl), INTENT(OUT) :: phiBC     											! scalar contribution from the boundary condition
INTEGER(lng) :: ip1,jp1,kp1 														! neighboring nodes (2 away from the wall)
REAL(dbl) :: q																			! distance ratio from the current node to the solid node
REAL(dbl) :: rhoB,phiB																! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: feq_m																	! equilibrium distribution function in the mth direction
REAL(dbl) :: phiijk_m																! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)
REAL(dbl) :: cosTheta, sinTheta													! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb																! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk ! radius of the solid node

CALL qCalc(m,i,j,k,im1,jm1,km1,q)												! calculate q	
!rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
rijk = x(im1)								! height at current location

cosTheta = x(im1)/rijk!r(km1)	! COS(theta)
sinTheta = y(jm1)/rijk!r(km1)	! SIN(theta)


IF (rijk .GE. rOut(k)) THEN
	ub = 0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
	vb = 0.0_dbl!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
	wb = velOut(km1)!vel(km1)!0.0_dbl						! only z-component in this case			
	!ub = -velOut(km1)*sinTheta!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
	!vb = velOut(km1)*cosTheta!0.0!vel(km1)*sinTheta							! y-component of the velocity at i,j,k
	!wb = 0.0_dbl!vel(km1)!0.0_dbl									! no z-component in this case			
	!ub = vel(km1)*cosTheta										! x-component of the velocity at i,j,k
	!vb = vel(km1)*sinTheta										! y-component of the velocity at i,j,k
	!wb = 0.0_dbl											! no z-component in this case
ELSE IF (rijk .LE. rIn(k)) THEN
	ub = 0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
	vb = 0.0_dbl!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
	wb = velIn(km1)!vel(km1)!0.0_dbl						! only z-component in this case	
	!ub = -velIn(km1)*sinTheta!0.0!vel(km1)*cosTheta					! x-component of the velocity at i,j,k
	!vb = velIn(km1)*cosTheta!0.0!vel(km1)*sinTheta					! y-component of the velocity at i,j,k
	!wb = 0.0_dbl!vel(km1)!0.0_dbl							! no z-component in this case			
	!ub = vel(km1)*cosTheta											! x-component of the velocity at i,j,k
	!vb = vel(km1)*sinTheta											! y-component of the velocity at i,j,k
	!wb = 0.0_dbl														! no z-component in this case
END IF		

! neighboring node (fluid side)	
ip1 = i + ex(m) 																		! i + 1
jp1 = j + ey(m)																		! j + 1
kp1 = k + ez(m)																		! k + 1

! if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF	

! assign values to boundary (density, scalar, f)
rhoB = (rho(i,j,k) - rho(ip1,jp1,kp1))*(1+q) + rho(ip1,jp1,kp1)		! extrapolate the density
CALL Equilibrium_LOCAL(m,rhoB,ub,vb,wb,feq_m)			        ! calculate the equibrium distribution function in the mth direction

!! Balaji added for sero flux BC. Otherwise set to constant value for Dirichlet BC
!phiWall = (phi(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) - (phi(ip1,jp1,kp1)*q*q/(1.0+2.0*q)) 	! calculate phiWall for flux BC (eq. 28 in paper)
!phiWall = (phiTemp(i,j,k)*(1.0_dbl+q-q*q) + phiTemp(ip1,jp1,kp1)*(q*q-2.0_dbl*q+1.0_dbl))/(2.0_dbl-q) 	! calculate phiWall for flux BC - Balaji's idea
!phiWall = (phiTemp(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) - (phiTemp(ip1,jp1,kp1)*q*q/(1.0+2.0*q)) 	! calculate phiWall for flux BC (eq. 28 in paper)

! find the contribution of scalar streamed from the wall to the current node (i,j,k), and from the current node to the next neighboring node (ip1,jp1,kp1)
phiB		= (feq_m/rhoB - wt(m)*Delta)*phiWall								! contribution from the wall in the mth direction (zero if phiWall=0)
phiijk_m	= (fplus(m,i,j,k)/rho(i,j,k) - wt(m)*Delta)*phiTemp(i,j,k)	! contribution from the current node to the next node in the mth direction

! if q is too small, the extrapolation to phiBC can create a large error...
IF(q .LT. 0.25) THEN
  q = 0.25_dbl  																		! approximate the distance ratio as 0.25
END IF

! extrapolate using phiB and phijk_m to obtain contribution from the solid node to the current node
phiBC		= ((phiB - phiijk_m)/q) + phiB										! extrapolated scalar value at the solid node, using q

!------------------------------------------------
END SUBROUTINE ScalarBC
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarBC2(m,i,j,k,im1,jm1,km1,phiBC,phiOut,phiIn)								! implements the scalar BCs by Balaji Dec 2014.
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1								! index variables
REAL(dbl), INTENT(OUT) :: phiBC,phiOut,phiIn     											! scalar contribution from the boundary condition
INTEGER(lng) :: ip1,jp1,kp1 														! neighboring nodes (2 away from the wall)
INTEGER(lng) :: ip2,jp2,kp2 														! neighboring nodes (3 away from the wall)
REAL(dbl) :: q																! distance ratio from the current node to the solid node
REAL(dbl) :: rhoB,phiB															! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: rhoAst,PAstToBSt,ScAst													! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes to an interior node 1 lattice unit away. 
REAL(dbl) :: rhoBst,PBstToCSt,ScBst,fBst												! values of density and at the boundary, and contribution of scalar from the B* node to C* node (wee yanxing wang's 2010 paper)
REAL(dbl) :: PBstToASt,fBstopp														! values of density and at the boundary, and contribution of scalar from the B* node to A* node (see yanxing wang's 2010 paper)
REAL(dbl) :: PAtoB, PAtoO, PBtoA, PCtoB, PBtoC, PCtoD,POtoA, Pmovingwall,fbb,fmoving															! values of density and at the boundary, and contribution of scalar from the B* node to A* node (see yanxing wang's 2010 paper)
REAL(dbl) :: rhoX,ScX
REAL(dbl) :: feq_m,feq_bbm,feq1_m,feq2_m,fnoneq_m,fnoneq1_m,fnoneq2_m									! equilibrium distribution function in the mth direction
REAL(dbl) :: phiijk_m	! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)
REAL(dbl) :: phiAst	! Scalar value at the boundary - For both Dirichlet and Neumann BC. For Dirichlet BC it is phiWall, but for Neumann BC it is computed from the flux BC.  
REAL(dbl) :: dphidn	! Neumann flux BC.  
REAL(dbl) :: cosTheta, sinTheta													! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb																! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk													! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt				! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 
INTEGER(lng) :: it			! loop index variables

LOGICAL :: BC2FLAG

! neighboring node (fluid side)	
ip1 = i + ex(m) 	! i + 1
jp1 = j + ey(m)		! j + 1
kp1 = k + ez(m)		! k + 1
ip2 = ip1 + ex(m) 	! i + 2
jp2 = jp1 + ey(m)	! j + 2
kp2 = kp1 + ez(m)	! k + 2

BC2FLAG = .FALSE.
IF ((node(ip1,jp1,kp1) .EQ. FLUID).AND.(node(ip2,jp2,kp2) .EQ. FLUID)) THEN
	BC2FLAG = .TRUE.
END IF
BC2FLAG = .FALSE.

!! if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
!IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
!  ip1 = i
!  jp1 = j
!  kp1 = k
!END IF
!
!!if (rho(ip1,jp1,kp1).lt.0.00001_dbl) THEN
!!	rho(ip1,jp1,kp1) = rho(i,j,k)!1.0_dbl
!!	phiTemp(ip1,jp1,kp1) = phiTemp(i,j,k)!1.0_dbl
!!	fplus(bb(m),ip1,jp1,kp1) = fplus(bb(m),i,j,k)
!!	fplus(m,ip1,jp1,kp1) = fplus(m,i,j,k)
!!ENDIF
!
!IF(node(ip2,jp2,kp2) .NE. FLUID) THEN
!  ip2 = ip1
!  jp2 = jp1
!  kp2 = kp1
!END IF	
!!if (rho(ip2,jp2,kp2).lt.0.00001_dbl) THEN
!!	rho(ip2,jp2,kp2) = rho(ip1,jp1,kp1)!1.0_dbl
!!	phiTemp(ip2,jp2,kp2) = phiTemp(ip1,jp1,kp1)!1.0_dbl
!!	fplus(bb(m),ip2,jp2,kp2) = fplus(bb(m),ip1,jp1,kp1)
!!	fplus(m,ip2,jp2,kp2) = fplus(m,ip1,jp1,kp1)
!!ENDIF



!! Original Code to estimate q
!CALL qCalc(m,i,j,k,im1,jm1,km1,q) ! calculate q
!!! if q is too small, the extrapolation to phiBC can create a large error...
!!IF(q .LT. 0.5) THEN
!!  q = 0.5_dbl  	! approximate the distance ratio as 0.25
!!END IF
!
!!q=1.0_dbl
!rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))! radius at current location
!
!cosTheta = x(im1)/rijk!r(km1)	! COS(theta)
!sinTheta = y(jm1)/rijk!r(km1)	! SIN(theta)
!
!
!IF (rijk .GE. rOut(k)) THEN
!	ub = -velOut(km1)*sinTheta!0.0!vel(km1)*cosTheta											! x-component of the velocity at i,j,k
!	vb = velOut(km1)*cosTheta!0.0!vel(km1)*sinTheta											! y-component of the velocity at i,j,k
!	wb = 0.0_dbl!vel(km1)!0.0_dbl														! no z-component in this case			
!	!ub = vel(km1)*cosTheta											! x-component of the velocity at i,j,k
!	!vb = vel(km1)*sinTheta											! y-component of the velocity at i,j,k
!	!wb = 0.0_dbl														! no z-component in this case
!ELSE IF (rijk .LE. rIn(k)) THEN
!	ub = -velIn(km1)*sinTheta!0.0!vel(km1)*cosTheta											! x-component of the velocity at i,j,k
!	vb = velIn(km1)*cosTheta!0.0!vel(km1)*sinTheta											! y-component of the velocity at i,j,k
!	wb = 0.0_dbl!vel(km1)!0.0_dbl														! no z-component in this case			
!	!ub = vel(km1)*cosTheta											! x-component of the velocity at i,j,k
!	!vb = vel(km1)*sinTheta											! y-component of the velocity at i,j,k
!	!wb = 0.0_dbl														! no z-component in this case
!END IF		
	

!BC2FLAG = .TRUE.
IF (BC2FLAG) THEN
! Calculate q using Yanxing's method
!*****************************************************************************
		!rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
		rijk = x(im1)								! height at current location
		 ! Initial fluid node guess
                 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                
		 ! Initial solid node guess
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
                 
	 IF (k.NE.km1) THEN
                 DO it=1,qitermax
		   ! guess of boundary location 
                   xt=(x1+x2)/2.0_dbl
                   yt=(y1+y2)/2.0_dbl
                   zt=(z1+z2)/2.0_dbl

      		   !rt = SQRT(xt*xt + yt*yt)
      		   rt = xt
		   !Write(*,*) 'test'
                   IF (node(im1,jm1,km1).EQ.SOLID2) THEN
		   	!ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
		   	ht = ((zt-z(k))*rOut(km1)+(z(km1)-zt)*rOut(k))/(z(km1)-z(k))
		   	!ht = (r(km1)+r(k))/2.0_dbl
		   ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
		   	!ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
		   	ht = ((zt-z(k))*rIn(km1)+(z(km1)-zt)*rIn(k))/(z(km1)-z(k))
		   	!ht = (r(km1)+r(k))/2.0_dbl
		   END IF

                   IF(rt.GT.ht) then
                     x2=xt
                     y2=yt
                     z2=zt
                   ELSE
                     x1=xt
                     y1=yt
                     z1=zt
                   END IF
				   
                 END DO
		 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                 
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
 
                 q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
		 !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
	 ELSE
		  DO it=1,qitermax
		   ! guess of boundary location 
                   xt=(x1+x2)/2.0_dbl
                   yt=(y1+y2)/2.0_dbl
                   zt=(z1+z2)/2.0_dbl

      		   !rt = SQRT(xt*xt + yt*yt)
      		   rt = xt
		   !Write(*,*) 'test'
		   IF (node(im1,jm1,km1).EQ.SOLID2) THEN
			   !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
			   !ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
			   ht = (rOut(km1)+rOut(k))/2.0_dbl
		   ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
			   !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
			   !ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
			   ht = (rIn(km1)+rIn(k))/2.0_dbl
		   END IF

                   IF(rt.GT.ht) then
                     x2=xt
                     y2=yt
                     z2=zt
                   ELSE
                     x1=xt
                     y1=yt
                     z1=zt
                   END IF
				   
                 END DO
		 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                 
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
 
                 q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
		 !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
	 ENDIF
		 cosTheta=xt/rt
		 sinTheta=yt/rt
	 IF (k.NE.km1) THEN
		   IF (node(im1,jm1,km1).EQ.SOLID2) THEN
			 !vt = (ABS(zt-z(k))*vel(km1)+ABS(z(km1)-zt)*vel(k))/ABS(z(km1)-z(k))
			 vt = ((zt-z(k))*velOut(km1)+(z(km1)-zt)*velOut(k))/(z(km1)-z(k))
		   ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
			 !vt = (ABS(zt-z(k))*vel(km1)+ABS(z(km1)-zt)*vel(k))/ABS(z(km1)-z(k))
			 vt = ((zt-z(k))*velIn(km1)+(z(km1)-zt)*velIn(k))/(z(km1)-z(k))
		   END IF
	 ELSE
		   IF (node(im1,jm1,km1).EQ.SOLID2) THEN
			 vt = (velOut(k)+velOut(km1))*0.5_dbl
		   ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
			 vt = (velIn(k)+velIn(km1))*0.5_dbl
		   END IF
	 ENDIF
		ub = 0.0_dbl!0.0!vel(km1)*cosTheta											! x-component of the velocity at i,j,k
		vb = 0.0_dbl!0.0!vel(km1)*sinTheta											! y-component of the velocity at i,j,k
		wb = vt!vel(km1)!0.0_dbl												! only z-component of velocity	
		!ub = -vt*sinTheta!0.0!vel(km1)*cosTheta											! x-component of the velocity at i,j,k
		!vb = vt*cosTheta!0.0!vel(km1)*sinTheta											! y-component of the velocity at i,j,k
		!wb = 0.0_dbl!vel(km1)!0.0_dbl			

	        ! make sure 0<q<1
        IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
          OPEN(1000,FILE="error.txt")
          WRITE(1000,*) "q=",q
          WRITE(1000,*) "m=",m
          WRITE(1000,*) "i=",i,"j=",j,"k=",k
          WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
          CLOSE(1000)
          STOP
        END IF	

	!q = max(q, 0.25_dbl)
ELSE

        ! Original Code to estimate q
        !CALL qCalc(m,i,j,k,im1,jm1,km1,q) ! calculate q
        !! if q is too small, the extrapolation to phiBC can create a large error...
        !IF(q .LT. 0.5) THEN
        !  q = 0.5_dbl  	! approximate the distance ratio as 0.25
        !END IF
        
        q=0.5_dbl
	!rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
	rijk = x(im1)								! height at current location
        
        cosTheta = x(im1)/rijk!r(km1)	! COS(theta)
        sinTheta = y(jm1)/rijk!r(km1)	! SIN(theta)
        
        
	IF (rijk .GE. rOut(k)) THEN
		ub = 0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
		vb = 0.0_dbl!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
		wb = velOut(km1)!vel(km1)!0.0_dbl						! only z-component in this case			
		!ub = -velOut(km1)*sinTheta!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
		!vb = velOut(km1)*cosTheta!0.0!vel(km1)*sinTheta							! y-component of the velocity at i,j,k
		!wb = 0.0_dbl!vel(km1)!0.0_dbl									! no z-component in this case			
		!ub = vel(km1)*cosTheta										! x-component of the velocity at i,j,k
		!vb = vel(km1)*sinTheta										! y-component of the velocity at i,j,k
		!wb = 0.0_dbl											! no z-component in this case
	ELSE IF (rijk .LE. rIn(k)) THEN
		ub = 0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
		vb = 0.0_dbl!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
		wb = velIn(km1)!vel(km1)!0.0_dbl						! only z-component in this case	
		!ub = -velIn(km1)*sinTheta!0.0!vel(km1)*cosTheta					! x-component of the velocity at i,j,k
		!vb = velIn(km1)*cosTheta!0.0!vel(km1)*sinTheta					! y-component of the velocity at i,j,k
		!wb = 0.0_dbl!vel(km1)!0.0_dbl							! no z-component in this case			
		!ub = vel(km1)*cosTheta											! x-component of the velocity at i,j,k
		!vb = vel(km1)*sinTheta											! y-component of the velocity at i,j,k
		!wb = 0.0_dbl														! no z-component in this case
	END IF			


END IF
!BC2FLAG = .FALSE.


!! neighboring node (fluid side)	
!ip1 = i + ex(m) 	! i + 1
!jp1 = j + ey(m)		! j + 1
!kp1 = k + ez(m)		! k + 1
!ip2 = ip1 + ex(m) 	! i + 2
!jp2 = jp1 + ey(m)	! j + 2
!kp2 = kp1 + ez(m)	! k + 2
!
!! if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
!IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
!  ip1 = i
!  jp1 = j
!  kp1 = k
!END IF
!
!if (rho(ip1,jp1,kp1).lt.0.00001_dbl) THEN
!	rho(ip1,jp1,kp1) = rho(i,j,k)!1.0_dbl
!	phiTemp(ip1,jp1,kp1) = phiTemp(i,j,k)!1.0_dbl
!	fplus(bb(m),ip1,jp1,kp1) = fplus(bb(m),i,j,k)
!	fplus(m,ip1,jp1,kp1) = fplus(m,i,j,k)
!ENDIF
!
!IF(node(ip2,jp2,kp2) .NE. FLUID) THEN
!  ip2 = ip1
!  jp2 = jp1
!  kp2 = kp1
!END IF	
!if (rho(ip2,jp2,kp2).lt.0.00001_dbl) THEN
!	rho(ip2,jp2,kp2) = rho(ip1,jp1,kp1)!1.0_dbl
!	phiTemp(ip2,jp2,kp2) = phiTemp(ip1,jp1,kp1)!1.0_dbl
!	fplus(bb(m),ip2,jp2,kp2) = fplus(bb(m),ip1,jp1,kp1)
!	fplus(m,ip2,jp2,kp2) = fplus(m,ip1,jp1,kp1)
!ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Yanxing's BC Scheme
!!q = max(q,0.01_dbl)
!! assign values to boundary (density, scalar, f)
!rhoAst = (rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)		! extrapolate the density
!
!CALL Equilibrium_LOCAL(m,rhoAst,ub,vb,wb,feq_m)			        ! calculate the equibrium distribution function in the mth direction
!CALL Equilibrium_LOCAL(m,rho(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),feq1_m)			        ! calculate the equibrium distribution function in the mth direction for node A
!CALL Equilibrium_LOCAL(m,rho(ip1,jp1,kp1),u(ip1,jp1,kp1),v(ip1,jp1,kp1),w(ip1,jp1,kp1),feq2_m)			        ! calculate the equibrium distribution function in the mth direction for node B
!
!!! Comment non-eq. effects for the time being. 
!fnoneq1_m = fplus(m,i,j,k)-feq1_m!max(fplus(m,i,j,k)-feq1_m,0.0_dbl)
!fnoneq2_m = fplus(m,ip1,jp1,kp1)-feq2_m!max(fplus(m,ip1,jp1,kp1)-feq2_m,0.0_dbl)
!fnoneq_m = fnoneq1_m+q*(fnoneq1_m-fnoneq2_m)
!feq_m = feq_m + fnoneq_m
!!write(9,*) feq_m, fnoneq_m, m
!!feq_m = (fplus(m,i,j,k) - fplus(m,ip1,jp1,kp1))*(q) + fplus(m,i,j,k)		! extrapolate the dist. function
!
!!! Balaji added for sero flux BC. Otherwise set to constant value for Dirichlet BC
!
!!phiWall = (phiTemp(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) - (phiTemp(ip1,jp1,kp1)*q*q/(1.0+2.0*q)) 	! calculate phiWall for flux BC (eq. 28 in paper)
!
!!phiWall = phiTemp(i,j,k) 	! calculate phiWall for flux BC (eq. 28 in paper)
!
!!phiWall = (phiTemp(i,j,k)*(1.0_dbl+q-q*q) + phiTemp(ip1,jp1,kp1)*(q*q-2.0_dbl*q+1.0_dbl))/(2.0_dbl-q) 	! calculate phiWall for flux BC - Balaji's idea
!
!!phiWall = min (phiWall,phiTemp(i,j,k))
!
!phiAst = phiWall
!!dphidn = 0.0_dbl
!!CALL GetPhiWall(i,j,k,m,q,dphidn,phiAst)
!dphidn = 0.0_dbl
!CALL GetPhiWallNew(i,j,k,m,q,dphidn,phiAst,rhoAst)
!rhoAst = (rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)		! extrapolate the density
!!phiAst = (phiTemp(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) - (phiTemp(ip1,jp1,kp1)*q*q/(1.0+2.0*q)) 	! calculate phiWall for flux BC (eq. 28 in paper)
!!phiAst = (phiTemp(i,j,k) - phiTemp(ip1,jp1,kp1))*(q)+ phiTemp(i,j,k) 	! calculate phiWall for flux BC (eq. 28 in paper)
!!phiAst = 2.0_dbl*q*phiTemp(i,j,k) + (1.0_dbl-2.0_dbl*q)*phiTemp(ip1,jp1,kp1) 	! calculate phiWall for flux BC (eq. 28 in paper)
!
!ScAst = phiAst
!
!! find the contribution of scalar streamed from the wall to the current node (i,j,k), and from the current node to the next neighboring node (ip1,jp1,kp1)
!PAstToBst = (feq_m/rhoAst - wt(m)*Delta)*ScAst								! contribution from the wall in the mth direction (zero if phiWall=0)
!
!rhoBst = rho(i,j,k) +(1.0_dbl-q)*(rho(ip1,jp1,kp1)-rho(i,j,k))		! extrapolate the density
!fBst = fplus(m,i,j,k) +(1.0_dbl-q)*(fplus(m,ip1,jp1,kp1)-fplus(m,i,j,k)) ! extrapolate the distribution function
!ScBst = phiTemp(i,j,k) +(1.0_dbl-q)*(phiTemp(ip1,jp1,kp1)-phiTemp(i,j,k)) ! extrapolate the scalar
!PBstToCst = (fBst/rhoBst - wt(m)*Delta)*ScBst	! contribution from the wall in the mth direction (zero if phiWall=0)
!
!fBstopp = fplus(bb(m),ip1,jp1,kp1) - q*(fplus(bb(m),ip1,jp1,kp1)-fplus(bb(m),i,j,k)) ! extrapolate the dist. function
!PBstToAst = (fBstopp/rhoBst - wt(bb(m))*Delta)*ScBst	! contribution to the wall in the bb(m) th direction 
!PAtoB = (fplus(m,i,j,k)/rho(i,j,k) - wt(m)*Delta)*phiTemp(i,j,k) ! contribution from A to B				
!
!! contribution from the wall in the mth direction (zero if phiWall=0)
!
!! if q is too small, the extrapolation to phiBC can create a large error...
!!zIF(q .LT. 0.25) THEN
!!  q = 0.25_dbl  ! approximate the distance ratio as 0.25
!!END IF
!
!! extrapolate using phiB and phijk_m to obtain contribution from the solid node to the current node
!phiBC	= PAstToBst - (PBstToCst - PAstToBst)*(1.0_dbl-q)
!!phiBC	= PAstToBst - (PAtoB - PAstToBst)*(1.0_dbl-q)/q
!phiOut = phiBC!PAstToBst
!phiIn = phiBC!PBstToAst
!
!!phiBtoA = ((fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k))!PBstToAst
!!phiAtoB = PAstToBst+wt(m)*Delta*ScAst
!!phiBtoA = PBstToAst+wt(bb(m))*Delta*ScBst
!
!!phiBC	= PBstToAst - (PBstToCst - PBstToAst)*(1.0_dbl-q)
!!!phiBC	= PAstToBst - (PAtoB - PAstToBst)*(1.0_dbl-q)/q
!!phiAtoB = PBstToAst
!!phiBtoA = PBstToAst
!
!!IF (m.EQ.1) THEN
!!	!write(9,*) i,j,k,ScAst,phiTemp(i,j,k)
!!	write(9,*) i,j,k,phiAtoB,phiBtoA
!!ENDIF
!
!!phiBC = pBstToAst
!!phiAtoB = PBstToAst
!!phiBC = 2.0_dbl*q*((fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k))+(1.0_dbl-2.0*dbl*q)*((fplus(bb(m),ip1,jp1,kp1)/rho(ip1,jp1,kp1) - wt(bb(m))*Delta)*phiTemp(ip1,jp1,kp1))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Using P to bounce back
!!rho(i,j,k) = 1.0_dbl
!!rho(ip1,jp1,kp1) = 1.0_dbl
!if (rho(ip2,jp2,kp2).lt.0.00001_dbl) THEN
!	rho(ip2,jp2,kp2) = rho(ip1,jp1,kp1)!1.0_dbl
!	phiTemp(ip2,jp2,kp2) = phiTemp(ip1,jp1,kp1)!1.0_dbl
!	fplus(bb(m),ip2,jp2,kp2) = fplus(bb(m),ip1,jp1,kp1)
!	fplus(m,ip2,jp2,kp2) = fplus(m,ip1,jp1,kp1)
!ENDIF
!!if (phiTemp(i,j,k).ne.0.0_dbl) then
!!write(*,*) i,j,k,iter
!!PAUSE
!!endif
!ScAst = (phiTemp(i,j,k) - phiTemp(ip1,jp1,kp1))*(q) + phiTemp(i,j,k)		! extrapolate the scalar
!!delta = SQRT(ex(bb(m))**2 + ex(bb(m))**2 + ex(bb(m))**2)*xcf
!!ScAst = ((q+1.0_dbl)*(q+1.0_dbl)*phiTemp(i,j,k)-q*q*phiTemp(ip1,jp1,kp1))!-dphidk*delta*q*(q+1.0_dbl))/(2.0_dbl*q+1.0_dbl)
!
!q = max(q,0.01_dbl)
!IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
!
!	! Interpolate the flux that leaves - pre streaming
!	PAtoO = 1.0_dbl*(fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta) &
!					*phiTemp(i,j,k)
!	PBtoA = 1.0_dbl*(fplus(bb(m),ip1,jp1,kp1)/rho(ip1,jp1,kp1) &
!			- wt(bb(m))*Delta)*phiTemp(ip1,jp1,kp1)
!	PCtoB = 1.0_dbl*(fplus(bb(m),ip2,jp2,kp2)/rho(ip2,jp2,kp2) - &
!			wt(bb(m))*Delta)*phiTemp(ip2,jp2,kp2)
!	Pmovingwall = 6.0_dbl*wt(m)*ScAst*(ub*ex(m) + vb*ey(m) + &
!			wb*ez(m))
!!	Pmovingwall = 6.0_dbl*wt(bb(m))*ScAst*(ub*ex(bb(m)) + vb*ey(bb(m)) + &
!!			wb*ez(bb(m)))
!	POtoA = q*(1.0_dbl+2.0_dbl*q)*PAtoO + (1.0_dbl-4.0_dbl*q*q)*PBtoA -  &
!			q*(1.0_dbl-2.0_dbl*q)*PCtoB &
!		+ Pmovingwall
!!write(*,*) phiBC, q, PAtoO,PBtoA,PCtoB
!!write(*,*) phiBC,q,phiTemp(i,j,k),phiTemp(ip1,jp1,kp1),phiTemp(ip2,jp2,kp2)
!!write(*,*) phiBC,q,fplus(bb(m),i,j,k),fplus(bb(m),ip1,jp1,kp1),fplus(bb(m),ip2,jp2,kp2)
!!write(*,*) phiBC,q,rho(i,j,k),rho(ip1,jp1,kp1),rho(ip2,jp2,kp2),node(ip2,jp2,kp2), FLUID
!ELSE IF((q.GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
!	! Interpolate the flux that arrives - pre streaming
!	PAtoO = 1.0_dbl*(fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta) &
!					*phiTemp(i,j,k)
!	PAtoB = 1.0_dbl*(fplus(m,i,j,k)/rho(i,j,k) &
!			- wt(m)*Delta)*phiTemp(i,j,k)
!	PBtoC = 1.0_dbl*(fplus(m,ip1,jp1,kp1)/rho(ip1,jp1,kp1) - &
!			wt(m)*Delta)*phiTemp(ip1,jp1,kp1)
!	Pmovingwall = 6.0_dbl*wt(m)*ScAst*(ub*ex(m) + vb*ey(m) + &
!			wb*ez(m))/(q*(1.0_dbl+2.0_dbl*q))
!!	Pmovingwall = 6.0_dbl*wt(bb(m))*ScAst*(ub*ex(bb(m)) + vb*ey(bb(m)) + &
!!			wb*ez(bb(m)))/(q*(1.0_dbl+2.0_dbl*q))
!	POtoA = (1.0_dbl/(q*(1.0_dbl+2.0_dbl*q)))*PAtoO + (2.0_dbl*q-1.0_dbl) &
!		*PAtoB/q -	(2.0_dbl*q-1.0_dbl)*PBtoC/(2.0_dbl*q+1.0_dbl) &
!		+ Pmovingwall 
!!write(*,*) phiBC, q, PAtoO,PBtoC,PCtoD
!ELSE
!    OPEN(1000,FILE='error-'//sub//'.txt')
!    WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
!    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
!    CLOSE(1000)
!    STOP
!END IF
!phiBC = POtoA
!
!PBtoC = 1.0_dbl*(fplus(m,ip1,jp1,kp1)/rho(ip1,jp1,kp1) &
!		- wt(m)*Delta)*phiTemp(ip1,jp1,kp1)
!PAtoB = 1.0_dbl*(fplus(m,i,j,k)/rho(i,j,k) &
!		- wt(m)*Delta)*phiTemp(i,j,k)
!PBtoA = 1.0_dbl*(fplus(bb(m),ip1,jp1,kp1)/rho(ip1,jp1,kp1) &
!		- wt(bb(m))*Delta)*phiTemp(ip1,jp1,kp1)
!
!phiOut = phiBC!POtoA!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!phiIn =  phiBC!PAtoO!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface
!
!!phiOut = POtoA!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!!phiIn =  PAtoO!PAtoO!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface
!!phiOut = POtoA!POtoA*(1.0+q)-PAtoB!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!!phiIn =  PAtoO - Pmovingwall!(PAtoO*q + PBtoA*(1.0_dbl-q))!PAtoO!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface
!!write(9,*) Pmovingwall,ex(m),ey(m),ez(m),PAtoO,POtoA
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!q = max(q,0.01_dbl)
! Using fbb to bounce back
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN

    IF (BC2FLAG) THEN

            rhoX = 2.0_dbl*q*rho(i,j,k) + (1.0_dbl-2.0_dbl*q)*rho(ip1,jp1,kp1)
            ScX  =  2.0_dbl*q*phiTemp(i,j,k) + (1.0_dbl-2.0_dbl*q)*phiTemp(ip1,jp1,kp1)
            rhoAst =  (rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)
            ScAst =  (phiTemp(i,j,k) - phiTemp(ip1,jp1,kp1))*(q) + phiTemp(i,j,k)
        
        ! 2nd order Lallemand and Luo
            !fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            !fmoving = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            !fmoving = 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))
            fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) &
                + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) & 
                - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) &
                + fmoving

    ELSE

        !! first order Balaji
        !    fmoving = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
        !    !fmoving = 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))
        !    fbb = (2.0_dbl*q)*fplus(bb(m),i,j,k) &
        !	  +(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip1,jp1,kp1) &
        !  	  + fmoving
        
        ! zeroth order half-way bounce back
            rhoX = rho(i,j,k)
            ScX  = phiTemp(i,j,k)
            rhoAst = rho(i,j,k)
            ScAst =  phiTemp(i,j,k)
            !fmoving = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))!*(rhoAst/ScAst)
            !fmoving = 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))
            fbb = fplus(bb(m),i,j,k) &
          	  + fmoving

    END IF

  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN

    IF (BC2FLAG) THEN

            rhoX = 2.0_dbl*q*rho(i,j,k) + (1.0_dbl-2.0_dbl*q)*rho(ip1,jp1,kp1)
            ScX  =  2.0_dbl*q*phiTemp(i,j,k) + (1.0_dbl-2.0_dbl*q)*phiTemp(ip1,jp1,kp1)
            !rhoX = rho(i,j,k)
            !ScX  = phiTemp(i,j,k)
            rhoAst =  (rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)
            ScAst =  (phiTemp(i,j,k) - phiTemp(ip1,jp1,kp1))*(q) + phiTemp(i,j,k)
        
        ! 2nd order Lallemand and Luo
            !fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
            !fmoving = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
            fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
            !fmoving = 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))/(q*(2.0_dbl*q + 1.0_dbl))
            fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 	&
                + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q	&
                - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1) &
                + fmoving

    ELSE

        !! first order Balaji
        !    fmoving = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(2.0_dbl*q)
        !    !fmoving = 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))/(2.0_dbl*q)
        !    fbb = fplus(bb(m),i,j,k)/(2.0_dbl*q) &
        !	+ (2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k)/(2.0_dbl*q) &
        !	+ fmoving
        
        ! zeroth order half-way bounce back
            rhoX = rho(i,j,k)
            ScX  = phiTemp(i,j,k)
            rhoAst = rho(i,j,k)
            ScAst =  phiTemp(i,j,k)
            !fmoving = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))!*(rhoAst/ScAst)
            !fmoving = 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))
            fbb = fplus(bb(m),i,j,k) &
        	+ fmoving

    END IF

  ELSE
    OPEN(1000,FILE='error-'//sub//'.txt')
    WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
    CLOSE(1000)
    STOP
  END IF

!rhoAst = rho(ip1,jp1,kp1)
!ScAst = phiTemp(ip1,jp1,kp1)
!rhoAst = (rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)		! extrapolate the density
!ScAst = (phiTemp(i,j,k) - phiTemp(ip1,jp1,kp1))*(q) + phiTemp(i,j,k)		! extrapolate the density
!rhoAst = 2.0_dbl*q*rho(i,j,k) + (1.0_dbl-2.0_dbl*q)*rho(ip1,jp1,kp1)		! extrapolate the density
!ScAst =  2.0_dbl*q*phiTemp(i,j,k) + (1.0_dbl-2.0_dbl*q)*phiTemp(ip1,jp1,kp1)		! extrapolate the density
!rhoAst = rhoX
!ScAst = ScX

!dphidn = 0.0_dbl
!CALL GetPhiWallNew(i,j,k,m,q,dphidn,phiAst,rhoAst)
!rhoAst = rhoX!(rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)		! extrapolate the density
!ScAst = ScX!phiAst!phiTemp(ip1,jp1,kp1)
!ScAst = max((2.0_dbl*q*phiAst + (1.0_dbl-q)*phiTemp(ip1,jp1,kp1))/(1.0_dbl+q),0.0_dbl) 
!write(9,*) ScX-ScAst,phiAst,ScX,ScAst

!phiBC = 1.0_dbl*(1.0*fbb*ScX/rhoX - wt(bb(m))*Delta*ScX) !+ wt(bb(m))*Delta*phiTemp(i,j,k)
!phiBC = 1.0_dbl*(1.0_dbl*fbb/rhoAst - wt(bb(m))*Delta)*ScAst !+ wt(bb(m))*Delta*phiTemp(i,j,k)
!phiBC = 1.0_dbl*(1.0_dbl*fplus(bb(m),i,j,k)/rhoAst - wt(bb(m))*Delta)*ScAst + fmoving*1.0_dbl !+ wt(bb(m))*Delta*phiTemp(i,j,k)
phiBC = 1.0_dbl*(1.0_dbl*(fbb-fmoving)/rhoX - wt(bb(m))*Delta)*ScX+(fmoving/rhoAst)*ScAst !+ wt(bb(m))*Delta*phiTemp(i,j,k)
!phiBC = 1.0_dbl*(1.0*fbb/rhoAst + wt(bb(m))*Delta)*ScAst - wt(bb(m))*Delta*phiTemp(i,j,k)
!phiBC = 1.0_dbl*(fbb/rhoX - wt(bb(m))*Delta)*ScX
phiOut = 1.0_dbl*(1.0_dbl*fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!phiOut = phiBC!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
phiIn =  phiBC!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface

!phiIn = 0.0_dbl
!phiOut = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))*(phiTemp(i,j,k)/(rho(i,j,k)+0.0e-10))


!!!!!! First order bounceback for scalar
!fbb = fplus(bb(m),i,j,k)+(6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
!rhoAst = rho(i,j,k)!(rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)		! extrapolate the density
!ScAst = phiTemp(i,j,k)!phiAst!phiTemp(ip1,jp1,kp1)
!phiBC = 1.0_dbl*(1.0_dbl*fbb/rhoAst - wt(bb(m))*Delta)*ScAst !+ wt(bb(m))*Delta*phiTemp(i,j,k)
!phiOut = 1.0_dbl*(1.0_dbl*fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!!phiOut = 1.0_dbl*(1.0_dbl*fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*ScAst!phiTemp(i,j,k)!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!phiIn =  phiBC!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface
!!phiOut = 0.0_dbl!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!!phiIn =  (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))*phiTemp(i,j,k)/rho(i,j,k)!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Using P to bounce back
!!rho(i,j,k) = 1.0_dbl
!!rho(ip1,jp1,kp1) = 1.0_dbl
!if (rho(ip2,jp2,kp2).lt.0.00001_dbl) THEN
!	rho(ip2,jp2,kp2) = rho(ip1,jp1,kp1)!1.0_dbl
!	phiTemp(ip2,jp2,kp2) = phiTemp(ip1,jp1,kp1)!1.0_dbl
!	fplus(bb(m),ip2,jp2,kp2) = fplus(bb(m),ip1,jp1,kp1)
!	fplus(m,ip2,jp2,kp2) = fplus(m,ip1,jp1,kp1)
!ENDIF
!!if (phiTemp(i,j,k).ne.0.0_dbl) then
!!write(*,*) i,j,k,iter
!!PAUSE
!!endif
!!ScAst = (phiTemp(i,j,k) - phiTemp(ip1,jp1,kp1))*(q) + phiTemp(i,j,k)		! extrapolate the scalar
!!delta = SQRT(ex(bb(m))**2 + ex(bb(m))**2 + ex(bb(m))**2)*xcf
!!ScAst = ((q+1.0_dbl)*(q+1.0_dbl)*phiTemp(i,j,k)-q*q*phiTemp(ip1,jp1,kp1))!-dphidk*delta*q*(q+1.0_dbl))/(2.0_dbl*q+1.0_dbl)
!dphidn = 0.0_dbl
!CALL GetPhiWallNew(i,j,k,m,q,dphidn,phiAst,rhoAst)
!ScAst = phiAst
!rhoAst = (rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)		! extrapolate the density
!
!q = max(q,0.01_dbl)
!IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
!
!	! Interpolate the flux that leaves - pre streaming
!	PAtoO = 1.0_dbl*(fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta) &
!					*phiTemp(i,j,k)
!	PBtoA = 1.0_dbl*(fplus(bb(m),ip1,jp1,kp1)/rho(ip1,jp1,kp1) &
!			- wt(bb(m))*Delta)*phiTemp(ip1,jp1,kp1)
!	PCtoB = 1.0_dbl*(fplus(bb(m),ip2,jp2,kp2)/rho(ip2,jp2,kp2) - &
!			wt(bb(m))*Delta)*phiTemp(ip2,jp2,kp2)
!	Pmovingwall = 6.0_dbl*wt(m)*ScAst*(ub*ex(m) + vb*ey(m) + &
!			wb*ez(m)) !-wt(m)*Delta*ScAst
!!	Pmovingwall = 6.0_dbl*wt(bb(m))*ScAst*(ub*ex(bb(m)) + vb*ey(bb(m)) + &
!!			wb*ez(bb(m))) !-wt(m)*Delta*ScAst
!	POtoA = q*(1.0_dbl+2.0_dbl*q)*PAtoO + (1.0_dbl-4.0_dbl*q*q)*PBtoA -  &
!			q*(1.0_dbl-2.0_dbl*q)*PCtoB &
!		+ Pmovingwall
!!write(*,*) phiBC, q, PAtoO,PBtoA,PCtoB
!!write(*,*) phiBC,q,phiTemp(i,j,k),phiTemp(ip1,jp1,kp1),phiTemp(ip2,jp2,kp2)
!!write(*,*) phiBC,q,fplus(bb(m),i,j,k),fplus(bb(m),ip1,jp1,kp1),fplus(bb(m),ip2,jp2,kp2)
!!write(*,*) phiBC,q,rho(i,j,k),rho(ip1,jp1,kp1),rho(ip2,jp2,kp2),node(ip2,jp2,kp2), FLUID
!ELSE IF((q.GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
!	! Interpolate the flux that arrives or received - pre streaming
!	PAtoO = 1.0_dbl*(fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta) & ! flux recd. by O or X after bounceback
!					*phiTemp(i,j,k)
!	PAtoB = 1.0_dbl*(fplus(m,i,j,k)/rho(i,j,k) & ! flux recd by B
!			- wt(m)*Delta)*phiTemp(i,j,k)
!	PBtoC = 1.0_dbl*(fplus(m,ip1,jp1,kp1)/rho(ip1,jp1,kp1) - & ! flux recd. by C
!			wt(m)*Delta)*phiTemp(ip1,jp1,kp1)
!	Pmovingwall = 6.0_dbl*wt(m)*ScAst*(ub*ex(m) + vb*ey(m) + & ! contribution form moving wall
!			wb*ez(m))/(q*(1.0_dbl+2.0_dbl*q)) !-wt(m)*Delta*ScAst
!!	Pmovingwall = 6.0_dbl*wt(bb(m))*ScAst*(ub*ex(bb(m)) + vb*ey(bb(m)) + &
!!			wb*ez(bb(m)))/(q*(1.0_dbl+2.0_dbl*q)) !-wt(m)*Delta*ScAst
!	POtoA = (1.0_dbl/(q*(1.0_dbl+2.0_dbl*q)))*PAtoO + (2.0_dbl*q-1.0_dbl) & ! flux recd. by A
!		*PAtoB/q -	(2.0_dbl*q-1.0_dbl)*PBtoC/(2.0_dbl*q+1.0_dbl) &
!		+ Pmovingwall 
!!write(*,*) phiBC, q, PAtoO,PBtoC,PCtoD
!ELSE
!    OPEN(1000,FILE='error-'//sub//'.txt')
!    WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
!    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
!    CLOSE(1000)
!    STOP
!END IF
!phiBC = POtoA
!
!!PBtoC = 1.0_dbl*(fplus(m,ip1,jp1,kp1)/rho(ip1,jp1,kp1) &
!!		- wt(m)*Delta)*phiTemp(ip1,jp1,kp1)
!!PAtoB = 1.0_dbl*(fplus(m,i,j,k)/rho(i,j,k) &
!!		- wt(m)*Delta)*phiTemp(i,j,k)
!!PBtoA = 1.0_dbl*(fplus(bb(m),ip1,jp1,kp1)/rho(ip1,jp1,kp1) &
!!		- wt(bb(m))*Delta)*phiTemp(ip1,jp1,kp1)
!
!phiOut = phiBC!POtoA!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!phiIn =  phiBC!PAtoO!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface
!
!!phiOut = POtoA!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!!phiIn =  PAtoO!PAtoO!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface
!!phiOut = POtoA!POtoA*(1.0+q)-PAtoB!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!!phiIn =  PAtoO - Pmovingwall!(PAtoO*q + PBtoA*(1.0_dbl-q))!PAtoO!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface
!!write(9,*) Pmovingwall,ex(m),ey(m),ez(m),PAtoO,POtoA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!phiBC = 1.0_dbl*(1.0*fplus(bb(m),i,j,k)/rho(i,j,k) - 1.0*wt(bb(m))*Delta)*phiTemp(i,j,k) &
!	+ 6.0_dbl*wt(m)*phiTemp(i,j,k)*(ub*ex(m) + vb*ey(m) + &
!			wb*ez(m))!- 1.0*wt(bb(m))*Delta*phiTemp(i,j,k)
!phiOut = phiBC!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
!phiIn =  phiBC!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface

!CALL Equilibrium_LOCAL(m,rhoB,ub,vb,wb,feq_m)			        ! calculate the equibrium distribution function in the mth direction
!phiBC = 0.9999_dbl*(fplus(bb(m),i,j,k)/rho(i,j,k) - 1.0*wt(bb(m))*Delta)*phiTemp(i,j,k)

!phiBC = min(phiBC,1.0_dbl*(fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k))

!------------------------------------------------
END SUBROUTINE ScalarBC2
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE GetPhiWallNew(ip0,jp0,kp0,m,q,dphidn,phiAst,rhoAst)	! Get phiAst for specified flux  New (Balaji, Jan 2015)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: ip0,jp0,kp0,m	! index variables
REAL(dbl), INTENT(IN) :: q     	! Boundary distance fraction
REAL(dbl) :: qlocal     	! Boundary distance fraction (local variable)
REAL(dbl), INTENT(IN) :: dphidn ! Specified boundary flux (gradient)
REAL(dbl), INTENT(OUT) :: phiAst,rhoAst ! scalar value at the wall
!REAL(dbl), INTENT(IN) :: phiA,phiB,phiC ! scalar values at the wall neighbours
REAL(dbl):: phiA,phiB,phiC ! scalar values at the wall neighbours
REAL(dbl):: rhoA,rhoB,rhoC ! scalar values at the wall neighbours
REAL(dbl) :: delta ! lattice distance in direction m
REAL(dbl) :: nx1,ny1,nz1 ! normal direction 1
REAL(dbl) :: px1,py1,pz1 ! surface point 1
REAL(dbl) :: px2,py2,pz2 ! surface point 2
REAL(dbl) :: px3,py3,pz3 ! surface point 3
REAL(dbl) :: px4,py4,pz4 ! surface point 4
REAL(dbl) :: magtvec,magnvec ! magnitude of vectors
REAL(dbl) :: dx,dy,dz,dist  ! direction vector components
INTEGER(lng):: i,j,k,im1,jm1,km1,ip1,jp1,kp1,ip2,jp2,kp2	! index variables
REAL(dbl),DIMENSION(3,3) :: arr,arrinv
REAL(dbl),DIMENSION(3) :: brhs
REAL(dbl) :: epsx,epsy,epsz,epsmin,veclen ! temp variables


! Solid node
im1 = ip0 - ex(m)
jm1 = jp0 - ey(m)
km1 = kp0 - ez(m)
! neighboring node (fluid side)	
ip1 = ip0 + ex(m) 	! i + 1
jp1 = jp0 + ey(m)	! j + 1
kp1 = kp0 + ez(m)	! k + 1
ip2 = ip1 + ex(m) 	! i + 2
jp2 = jp1 + ey(m)	! j + 2
kp2 = kp1 + ez(m)	! k + 2


!open(9,file='testoutput.dat',position='append')

! get normal direction 
! to do this, we will need at the least 3 points on the surface plane. 

!! Approach 1: We get the 3-4 surface points and normal vector from analytical surface equations
!!!! Compute normal from analytical surface
epsmin = 1.0e-15_dbl


! Get first surface point from Ray tracing
epsx=0.0_dbl
epsy=0.0_dbl
epsz=0.0_dbl
CALL IntersectSurface(m,ip0,jp0,kp0,im1,jm1,km1,px1,py1,pz1,qlocal,epsz,epsy,epsx)
!write(9,*) "px1,py1,pz1",px1,py1,pz1,qlocal
!write(9,*) q,qlocal
!CALL FLUSH(9)

nx1 = 2.0_dbl*px1
ny1 = 2.0_dbl*py1
IF (k.GT.1) THEN
	IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	nz1 = -(2.0_dbl*((rOut(k)-rOut(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*rOut(k-1)/zcf+(pz1-(z(k-1)/zcf))*rOut(k)/zcf)/((z(k)-z(k-1))/zcf))
	ELSE
	nz1 = -(2.0_dbl*((rIn(k)-rIn(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*rIn(k-1)/zcf+(pz1-(z(k-1)/zcf))*rIn(k)/zcf)/((z(k)-z(k-1))/zcf))
	ENDIF
ELSE 
	IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	nz1 = -(2.0_dbl*((rOut(k+1)-rOut(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*rOut(k)/zcf+(pz1-(z(k)/zcf))*rOut(k+1)/zcf)/((z(k+1)-z(k))/zcf))
	ELSE
	nz1 = -(2.0_dbl*((rIn(k+1)-rIn(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*rIn(k)/zcf+(pz1-(z(k)/zcf))*rIn(k+1)/zcf)/((z(k+1)-z(k))/zcf))
	END IF

ENDIF
magnvec = SQRT(nx1**2+ny1**2+nz1**2)
nx1 = nx1/magnvec
ny1 = ny1/magnvec
nz1 = nz1/magnvec

!! Make sure the normal vector is an outer normal
!IF (nx1*ex(bb(m))+ny1*ey(bb(m))+nz1*ez(bb(m)).LT.0.0_dbl) THEN
!nx1 = -nx1
!ny1 = -ny1
!nz1 = -nz1
!ENDIF
!write(9,*) "normal vector ",nx1,ny1,nz1,ex(bb(m)),ey(bb(m)),ez(bb(m))

!write(9,*) 'ip0,jp0,kp0',ip0,jp0,kp0,node(ip0,jp0,kp0)
! Get points along the normal vector in the fluid side of the surface. 
veclen = 0.1_dbl
px2 = px1 - veclen*nx1
py2 = py1 - veclen*ny1
pz2 = pz1 - veclen*nz1
!write(9,*) px2-REAL(iMin-2_lng)+REAL(Ci - 1_lng),py2-REAL(jMin-2_lng)+REAL(Cj - 1_lng),pz2+ 0.5_dbl- REAL(kMin-1_lng),ip1,jp1,kp1,ip0,jp0,kp0,px1-REAL(iMin-2_lng)+REAL(Ci - 1_lng),py1-REAL(jMin-2_lng)+REAL(Cj - 1_lng),pz1+ 0.5_dbl- REAL(kMin-1_lng)
CALL InterpolateProp(px2,py2,pz2,phiA,rhoA)

px3 = px2 - veclen*nx1
py3 = py2 - veclen*ny1
pz3 = pz2 - veclen*nz1
	!write(9,*) px3,py3,pz3
!write(9,*) px3-REAL(iMin-2_lng)+REAL(Ci - 1_lng),py3-REAL(jMin-2_lng)+REAL(Cj - 1_lng),pz3+ 0.5_dbl- REAL(kMin-1_lng)
CALL InterpolateProp(px3,py3,pz3,phiB,rhoB)
!check if these are fluid points. 

!phiAst = phiTemp(ip1,jp1,kp1)
!delta = SQRT(ex(bb(m))**2 + ex(bb(m))**2 + ex(bb(m))**2)*xcf
!phiAst = ((q+1.0_dbl)*(q+1.0_dbl)*phiA-(q)*(q)*phiB-dphidn*delta*q*(q+1.0_dbl))/(2.0_dbl*q+1.0_dbl)
!rhoAst = ((q+1.0_dbl)*(q+1.0_dbl)*rhoA-q*q*rhoB-0.0_dbl*delta*q*(q+1.0_dbl))/(2.0_dbl*q+1.0_dbl)
delta = veclen
!phiAst = (9.0_dbl*phiA - 4.0_dbl*phiB - 6.0*delta*dphidn)/5.0_dbl
!rhoAst = (9.0_dbl*rhoA - 4.0_dbl*rhoB - 6.0*delta*0.0_dbl)/5.0_dbl
phiAst = (4.0_dbl*phiA - 1.0_dbl*phiB - 2.0*delta*dphidn)/3.0_dbl
rhoAst = (4.0_dbl*rhoA - 1.0_dbl*rhoB - 2.0*delta*0.0_dbl)/3.0_dbl

!------------------------------------------------
END SUBROUTINE GetPhiWallNew
!------------------------------------------------

!------------------------------------------------
SUBROUTINE InterpolateProp(xp,yp,zp,phip,rhop) ! Interpolate solution field w to location xp,yp,zp 
!------------------------------------------------
IMPLICIT NONE
REAL(dbl), INTENT(IN) :: xp,yp,zp	! coordinates
REAL(dbl) :: xps,yps,zps	! coordinates
REAL(dbl), INTENT(OUT) :: phip,rhop   	! Interpolated field value being returned 
INTEGER(lng) :: ix0,iy0,iz0,ix1,iy1,iz1,numFLUIDs
REAL(dbl) :: xd,yd,zd
REAL(dbl) :: c00,c01,c10,c11,c,c0,c1
REAL(dbl) :: rhoSum,phiSum

	xps = xp-REAL(iMin-2_lng)+REAL(Ci - 1_lng)
	ix0=FLOOR(xps) 
	ix1=ix0 + 1_lng
        yps = yp-REAL(jMin-2_lng)+REAL(Cj - 1_lng)
	iy0=FLOOR(yps)
	iy1=iy0 + 1_lng
        zps = zp + 0.5_dbl- REAL(kMin-1_lng)
	iz0=FLOOR(zps)
	iz1=iz0+1_lng
	xd=(xps-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
	yd=(yps-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
	zd=(zps-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))

	!write(9,*) 'ix0,iy0,iz0',ix0,iy0,iz0
	!write(9,*) 'xd,yd,zd',xd,yd,zd
	!write(9,*) "xps,yps,zps",xps,yps,zps

	! w-interpolation
	! Do first level linear interpolation in x-direction
!	IF ((rho(ix0,iy0,iz0).GT.0.000001_dbl).OR.(rho(ix1,iy0,iz0).GT.0.000001_dbl) &
!	.OR.(rho(ix0,iy0,iz1).GT.0.000001_dbl).OR.(rho(ix1,iy0,iz1).GT.0.000001_dbl) &
!	.OR.(rho(ix0,iy1,iz0).GT.0.000001_dbl).OR.(rho(ix1,iy1,iz0).GT.0.000001_dbl) &
!	.OR.(rho(ix0,iy1,iz1).GT.0.000001_dbl).OR.(rho(ix1,iy1,iz1).GT.0.000001_dbl)) THEN
	IF ((node(ix0,iy0,iz0).EQ. FLUID).OR.(node(ix1,iy0,iz0).EQ. FLUID) &
	.OR.(node(ix0,iy0,iz1).EQ. FLUID).OR.(node(ix1,iy0,iz1).EQ. FLUID) &
	.OR.(node(ix0,iy1,iz0).EQ. FLUID).OR.(node(ix1,iy1,iz0).EQ. FLUID) &
	.OR.(node(ix0,iy1,iz1).EQ. FLUID).OR.(node(ix1,iy1,iz1).EQ. FLUID)) THEN
		c00 = phiTemp(ix0,iy0,iz0)*(1.0_dbl-xd)+phiTemp(ix1,iy0,iz0)*xd	
		c01 = phiTemp(ix0,iy0,iz1)*(1.0_dbl-xd)+phiTemp(ix1,iy0,iz1)*xd	
		c10 = phiTemp(ix0,iy1,iz0)*(1.0_dbl-xd)+phiTemp(ix1,iy1,iz0)*xd	
		c11 = phiTemp(ix0,iy1,iz1)*(1.0_dbl-xd)+phiTemp(ix1,iy1,iz1)*xd	
		! Do second level linear interpolation in y-direction
		c0  = c00*(1.0_dbl-yd)+c10*yd
		c1  = c01*(1.0_dbl-yd)+c11*yd
		! Do third level linear interpolation in z-direction
		c   = c0*(1.0_dbl-zd)+c1*zd
	        phip=c

		c00 = rho(ix0,iy0,iz0)*(1.0_dbl-xd)+rho(ix1,iy0,iz0)*xd	
		c01 = rho(ix0,iy0,iz1)*(1.0_dbl-xd)+rho(ix1,iy0,iz1)*xd	
		c10 = rho(ix0,iy1,iz0)*(1.0_dbl-xd)+rho(ix1,iy1,iz0)*xd	
		c11 = rho(ix0,iy1,iz1)*(1.0_dbl-xd)+rho(ix1,iy1,iz1)*xd	
		! Do second level linear interpolation in y-direction
		c0  = c00*(1.0_dbl-yd)+c10*yd
		c1  = c01*(1.0_dbl-yd)+c11*yd
		! Do third level linear interpolation in z-direction
		c   = c0*(1.0_dbl-zd)+c1*zd
	        rhop=c
	ELSE
		OPEN(1000,FILE='error-'//sub//'.txt')
		WRITE(1000,*) "Error in InterpolateProp in ICBC.f90 (line 1736): rhos are zero"
		WRITE(1000,*) "rho(i,j,k) is zero",rho(ix0,iy0,iz0),rho(ix1,iy0,iz0),rho(ix0,iy0,iz1) &
		,rho(ix1,iy0,iz1),rho(ix0,iy1,iz0),rho(ix1,iy1,iz0),rho(ix0,iy1,iz1),rho(ix1,iy1,iz1)
		WRITE(1000,*) "node(i,j,k) is ",node(ix0,iy0,iz0),node(ix1,iy0,iz0),node(ix0,iy0,iz1) &
		,node(ix1,iy0,iz1),node(ix0,iy1,iz0),node(ix1,iy1,iz0),node(ix0,iy1,iz1),node(ix1,iy1,iz1)
		WRITE(1000,*) "node index",ix0,iy0,iz0,ix1,iy1,iz1
		CLOSE(1000)
		STOP
	ENDIF
         
! 	 rhoSum = 0.0_dbl
!         phiSum = 0.0_dbl
!         numFLUIDs = 0_lng	
!	 !IF ((rho(ix0,iy0,iz0).GT.0.0000001_dbl)) THEN
!	 IF ((node(ix0,iy0,iz0).EQ. FLUID)) THEN
!	      rhoSum = rhoSum + rho(ix0,iy0,iz0)
!	      !rhoSum = max(rhoSum,rho(ii,jj,kk))
!	      phiSum = phiSum + phi(ix0,iy0,iz0)
!	      !phiSum = max(phiSum,phi(ii,jj,kk))!phiSum + phi(ii,jj,kk)
!	      numFLUIDs = numFLUIDs + 1_lng  
!	 END IF
!	 !IF ((rho(ix1,iy0,iz0).GT.0.0000001_dbl)) THEN
!	 IF ((node(ix1,iy0,iz0).EQ. FLUID)) THEN
!	      rhoSum = rhoSum + rho(ix1,iy0,iz0)
!	      !rhoSum = max(rhoSum,rho(ii,jj,kk))
!	      phiSum = phiSum + phi(ix1,iy0,iz0)
!	      !phiSum = max(phiSum,phi(ii,jj,kk))!phiSum + phi(ii,jj,kk)
!	      numFLUIDs = numFLUIDs + 1_lng  
!	 END IF
!	 !IF ((rho(ix0,iy0,iz1).GT.0.0000001_dbl)) THEN
!	 IF ((node(ix0,iy0,iz1).EQ. FLUID)) THEN
!	      rhoSum = rhoSum + rho(ix0,iy0,iz1)
!	      !rhoSum = max(rhoSum,rho(ii,jj,kk))
!	      phiSum = phiSum + phi(ix0,iy0,iz1)
!	      !phiSum = max(phiSum,phi(ii,jj,kk))!phiSum + phi(ii,jj,kk)
!	      numFLUIDs = numFLUIDs + 1_lng  
!	 END IF
!	 !IF ((rho(ix1,iy0,iz1).GT.0.0000001_dbl)) THEN
!	 IF ((node(ix1,iy0,iz1).EQ. FLUID)) THEN
!	      rhoSum = rhoSum + rho(ix1,iy0,iz1)
!	      !rhoSum = max(rhoSum,rho(ii,jj,kk))
!	      phiSum = phiSum + phi(ix1,iy0,iz1)
!	      !phiSum = max(phiSum,phi(ii,jj,kk))!phiSum + phi(ii,jj,kk)
!	      numFLUIDs = numFLUIDs + 1_lng  
!	 END IF
!	 !IF ((rho(ix0,iy1,iz0).GT.0.0000001_dbl)) THEN
!	 IF ((node(ix0,iy1,iz0).EQ. FLUID)) THEN
!	      rhoSum = rhoSum + rho(ix0,iy1,iz0)
!	      !rhoSum = max(rhoSum,rho(ii,jj,kk))
!	      phiSum = phiSum + phi(ix0,iy1,iz0)
!	      !phiSum = max(phiSum,phi(ii,jj,kk))!phiSum + phi(ii,jj,kk)
!	      numFLUIDs = numFLUIDs + 1_lng  
!	 END IF
!	 !IF ((rho(ix1,iy1,iz0).GT.0.0000001_dbl)) THEN
!	 IF ((node(ix1,iy1,iz0).EQ. FLUID)) THEN
!	      rhoSum = rhoSum + rho(ix1,iy1,iz0)
!	      !rhoSum = max(rhoSum,rho(ii,jj,kk))
!	      phiSum = phiSum + phi(ix1,iy1,iz0)
!	      !phiSum = max(phiSum,phi(ii,jj,kk))!phiSum + phi(ii,jj,kk)
!	      numFLUIDs = numFLUIDs + 1_lng  
!	 END IF
!	 !IF ((rho(ix0,iy1,iz1).GT.0.0000001_dbl)) THEN
!	 IF ((node(ix0,iy1,iz1).EQ. FLUID)) THEN
!	      rhoSum = rhoSum + rho(ix0,iy1,iz1)
!	      !rhoSum = max(rhoSum,rho(ii,jj,kk))
!	      phiSum = phiSum + phi(ix0,iy1,iz1)
!	      !phiSum = max(phiSum,phi(ii,jj,kk))!phiSum + phi(ii,jj,kk)
!	      numFLUIDs = numFLUIDs + 1_lng  
!	 END IF
!	 !IF ((rho(ix1,iy1,iz1).GT.0.0000001_dbl)) THEN
!	 IF ((node(ix1,iy1,iz1).EQ. FLUID)) THEN
!	      rhoSum = rhoSum + rho(ix1,iy1,iz1)
!	      !rhoSum = max(rhoSum,rho(ii,jj,kk))
!	      phiSum = phiSum + phi(ix1,iy1,iz1)
!	      !phiSum = max(phiSum,phi(ii,jj,kk))!phiSum + phi(ii,jj,kk)
!	      numFLUIDs = numFLUIDs + 1_lng  
!	 END IF
!
!	IF(numFLUIDs .NE. 0_lng) THEN
!	 rhop = rhoSum/numFLUIDs
!	 phip = phiSum/numFLUIDs
!	ELSE
!		OPEN(1000,FILE='error-'//sub//'.txt')
!		WRITE(1000,*) "Error in InterpolateProp in ICBC.f90 (line 1736): rhos are zero"
!		WRITE(1000,*) "rho(i,j,k) is zero",rho(ix0,iy0,iz0),rho(ix1,iy0,iz0),rho(ix0,iy0,iz1) &
!		,rho(ix1,iy0,iz1),rho(ix0,iy1,iz0),rho(ix1,iy1,iz0),rho(ix0,iy1,iz1),rho(ix1,iy1,iz1)
!		WRITE(1000,*) "node(i,j,k) is ",node(ix0,iy0,iz0),node(ix1,iy0,iz0),node(ix0,iy0,iz1) &
!		,node(ix1,iy0,iz1),node(ix0,iy1,iz0),node(ix1,iy1,iz0),node(ix0,iy1,iz1),node(ix1,iy1,iz1)
!		WRITE(1000,*) "node index",ix0,iy0,iz0,ix1,iy1,iz1
!		CLOSE(1000)
!		STOP
!	END IF



!------------------------------------------------
END SUBROUTINE InterpolateProp
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE GetPhiWall(ip0,jp0,kp0,m,q,dphidn,phiAst)	! Get phiAst for specified flux  New (Balaji, Dec 2014)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: ip0,jp0,kp0,m	! index variables
REAL(dbl), INTENT(IN) :: q     	! Boundary distance fraction
REAL(dbl) :: qlocal     	! Boundary distance fraction (local variable)
REAL(dbl), INTENT(IN) :: dphidn ! Specified boundary flux (gradient)
REAL(dbl), INTENT(OUT) :: phiAst ! scalar value at the wall
!REAL(dbl), INTENT(IN) :: phiA,phiB,phiC ! scalar values at the wall neighbours
REAL(dbl):: phiA,phiB,phiC ! scalar values at the wall neighbours
REAL(dbl) :: dphidtAst1,dphidtA1,dphidtB1,dphidtC1 ! tangential gradient of phi (1st tangent vector)
REAL(dbl) :: dphidtAst2,dphidtA2,dphidtB2,dphidtC2 ! tangential gradient of phi (2nd tangent vector)
REAL(dbl) :: delta ! lattice distance in direction m
REAL(dbl) :: tx1,ty1,tz1 ! tangent direction 1
REAL(dbl) :: tx2,ty2,tz2 ! tangent direction 2
REAL(dbl) :: nx1,ny1,nz1 ! normal direction 1
REAL(dbl) :: px1,py1,pz1 ! surface point 1
REAL(dbl) :: px2,py2,pz2 ! surface point 2
REAL(dbl) :: px3,py3,pz3 ! surface point 3
REAL(dbl) :: px4,py4,pz4 ! surface point 4
REAL(dbl) :: sx1,sy1,sz1 ! surface vector 1
REAL(dbl) :: sx2,sy2,sz2 ! surface vector 2
REAL(dbl) :: magtvec,magnvec ! magnitude of vectors
REAL(dbl) :: dx,dy,dz,dist  ! direction vector components
REAL(dbl) :: epsx,epsy,epsz,eps,epsmin  ! deviation vector components
REAL(dbl) :: dphidx,dphidy,dphidz,dphidk ! grad(phi) components
INTEGER(lng):: i,j,k,im1,jm1,km1,ip1,jp1,kp1,ip2,jp2,kp2	! index variables
REAL(dbl),DIMENSION(3,3) :: arr,arrinv
REAL(dbl),DIMENSION(3) :: brhs
REAL(dbl) :: temp,aa1,aa2,aa3,bb1,bb2,bb3,cc1,cc2,cc3,rb1,rb2,rb3,t1,t2,sz,sx,sy

!REAL(dbl) :: term1,term2,term3  ! temp variables
!REAL(dbl) :: randx1,randy1,randz1 ! temp variables
!REAL(dbl) :: randx2,randy2,randz2 ! temp variables
!REAL(dbl) :: s,r,h,theta,phiang,sinT,xc,yc,zc,newdx,newdy,newdz,randnum ! temp variables

! Solid node
im1 = ip0 - ex(m)
jm1 = jp0 - ey(m)
km1 = kp0 - ez(m)
! neighboring node (fluid side)	
ip1 = ip0 + ex(m) 	! i + 1
jp1 = jp0 + ey(m)	! j + 1
kp1 = kp0 + ez(m)	! k + 1
ip2 = ip1 + ex(m) 	! i + 2
jp2 = jp1 + ey(m)	! j + 2
kp2 = kp1 + ez(m)	! k + 2


!open(9,file='testoutput.dat',position='append')

! get normal direction 
! to do this, we will need at the least 3 points on the surface plane. 

!! Approach 1: We get the 3-4 surface points and normal vector from analytical surface equations
!!!! Compute normal from analytical surface
epsmin = 1.0e-15_dbl

! Get first surface point from Ray tracing
epsx=0.0_dbl
epsy=0.0_dbl
epsz=0.0_dbl
CALL IntersectSurface(m,ip0,jp0,kp0,im1,jm1,km1,px1,py1,pz1,qlocal,epsz,epsy,epsx)
!write(9,*) "px1,py1,pz1",px1,py1,pz1,qlocal
!CALL FLUSH(9)

nx1 = 2.0_dbl*px1
ny1 = 2.0_dbl*py1
IF (k.GT.1) THEN
	IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	nz1 = -(2.0_dbl*((rOut(k)-rOut(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*rOut(k-1)/zcf+(pz1-(z(k-1)/zcf))*rOut(k)/zcf)/((z(k)-z(k-1))/zcf))
	ELSE
	nz1 = -(2.0_dbl*((rIn(k)-rIn(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*rIn(k-1)/zcf+(pz1-(z(k-1)/zcf))*rIn(k)/zcf)/((z(k)-z(k-1))/zcf))
	ENDIF
ELSE 
	IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	nz1 = -(2.0_dbl*((rOut(k+1)-rOut(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*rOut(k)/zcf+(pz1-(z(k)/zcf))*rOut(k+1)/zcf)/((z(k+1)-z(k))/zcf))
	ELSE
	nz1 = -(2.0_dbl*((rIn(k+1)-rIn(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*rIn(k)/zcf+(pz1-(z(k)/zcf))*rIn(k+1)/zcf)/((z(k+1)-z(k))/zcf))
	END IF

ENDIF
magnvec = SQRT(nx1**2+ny1**2+nz1**2)
nx1 = nx1/magnvec
ny1 = ny1/magnvec
nz1 = nz1/magnvec

!! Make sure the normal vector is an outer normal
!IF (nx1*ex(bb(m))+ny1*ey(bb(m))+nz1*ez(bb(m)).LT.0.0_dbl) THEN
!nx1 = -nx1
!ny1 = -ny1
!nz1 = -nz1
!ENDIF
!write(9,*) "normal vector1 old",magnvec,nx1,ny1,nz1,nx1*sx1+ny1*sy1+nz1*sz1,nx1*sx2+ny1*sy2+nz1*sz2

! get 2nd point analytically
epsx = 0.0001_dbl
IF (abs(ny1+nz1).GT.epsmin) THEN
	epsy = -(nx1*epsx)/(ny1+nz1)
ELSE
	epsy = -(nx1*epsx)/sign(epsmin,ny1+nz1)
ENDIF
px2 = px1 + epsx
py2 = py1 + epsy
pz2 = pz1 + epsz

!! get 3rd point analytically
!epsx = -0.01_dbl
!IF (abs(ny1+nz1).GT.epsmin) THEN
!	epsy = -(nx1*epsx)/(ny1+nz1)
!ELSE
!	epsy = -(nx1*epsx)/sign(epsmin,ny1+nz1)
!ENDIF
!px3 = px1 + epsx
!py3 = py1 + epsy
!pz3 = pz1 + epsz

! Calculate first surface vector after obtaining the points
sx1 = px1-px2
sy1 = py1-py2
sz1 = pz1-pz2
magtvec = SQRT(sx1**2+sy1**2+sz1**2)
!write(9,*) sx1,sy1,sz1,magtvec
!CALL FLUSH(9)
sx1 = sx1/magtvec
sy1 = sy1/magtvec
sz1 = sz1/magtvec

!! Calculate second surface vector after obtaining the points
!sx2 = px3-px2
!sy2 = py3-py2
!sz2 = pz3-pz2
!magtvec = SQRT(sx2**2+sy2**2+sz2**2)
!!write(9,*) sx2,sy2,sz2,magtvec
!!CALL FLUSH(9)
!sx2 = sx2/magtvec
!sy2 = sy2/magtvec
!sz2 = sz2/magtvec




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Approach 2: We get the 3-4 surface points and normal vector from Ray Tracing Approach
!epsmin = 1.0e-8_dbl
!
!dist = SQRT(ex(bb(m))**2+ey(bb(m))**2+ez(bb(m))**2)
!dx=ex(bb(m))/dist	! i direction
!dy=ey(bb(m))/dist	! j direction
!dz=ez(bb(m))/dist	! k direction
!
!epsx=0.0_dbl
!epsy=0.0_dbl
!epsz=0.0_dbl
!CALL IntersectSurface(m,ip0,jp0,kp0,im1,jm1,km1,px1,py1,pz1,qlocal,epsz,epsy,epsx)
!!write(9,*) "px1,py1,pz1",px1,py1,pz1,qlocal
!!CALL FLUSH(9)
!
!CALL GetPerturbedVector(dx,dy,dz,epsx,epsy,epsz)
!CALL IntersectSurface(m,ip0,jp0,kp0,im1,jm1,km1,px2,py2,pz2,qlocal,epsz,epsy,epsx)
!!write(9,*) "px2,py2,pz2",px2,py2,pz2,qlocal
!!CALL FLUSH(9)
!
!
!CALL GetPerturbedVector(dx,dy,dz,epsx,epsy,epsz)
!CALL IntersectSurface(m,ip0,jp0,kp0,im1,jm1,km1,px3,py3,pz3,qlocal,epsz,epsy,epsx)
!!write(9,*) "px3,py3,pz3",px3,py3,pz3,qlocal
!!CALL FLUSH(9)
!
!!CALL GetPerturbedVector(dx,dy,dz,epsx,epsy,epsz)
!!CALL IntersectSurface(m,ip0,jp0,kp0,im1,jm1,km1,px4,py4,pz4,qlocal,epsz,epsy,epsx)
!!!!write(9,*) "px4,py4,pz4",px4,py4,pz4,qlocal
!!!!CALL FLUSH(9)
!
!! Calculate first surface vector after obtaining the points
!sx1 = px3-px1
!sy1 = py3-py1
!sz1 = pz3-pz1
!magtvec = SQRT(sx1**2+sy1**2+sz1**2)
!!write(9,*) sx1,sy1,sz1,magtvec
!!CALL FLUSH(9)
!sx1 = sx1/magtvec
!sy1 = sy1/magtvec
!sz1 = sz1/magtvec
!
!! Calculate second surface vector after obtaining the points
!sx2 = px3-px2
!sy2 = py3-py2
!sz2 = pz3-pz2
!magtvec = SQRT(sx2**2+sy2**2+sz2**2)
!!write(9,*) sx2,sy2,sz2,magtvec
!!CALL FLUSH(9)
!sx2 = sx2/magtvec
!sy2 = sy2/magtvec
!sz2 = sz2/magtvec
!
!!write(9,*) "tangent plane vectors", sx1,sy1,sz1,sx2,sy2,sz2
!
!!! Calculate surface normal direction
!! Approach 1: from cross product of surface vectors 
!!nx1 = (sy1*sz2-sz1*sy2)
!!ny1 = -(sx1*sz2-sz1*sx2)
!!nz1 = (sx1*sy2-sy1*sx2)
!!magnvec = SQRT(nx1**2+ny1**2+nz1**2)
!!!write(9,*) nx1,ny1,nz1,magnvec
!!!CALL FLUSH(9)
!!nx1 = nx1/magnvec
!!ny1 = ny1/magnvec
!!nz1 = nz1/magnvec
!!
!!! Make sure the normal vector is an outer normal
!!IF (nx1*ex(bb(m))+ny1*ey(bb(m))+nz1*ez(bb(m)).LT.0.0_dbl) THEN
!!nx1 = -nx1
!!ny1 = -ny1
!!nz1 = -nz1
!!ENDIF
!!!write(9,*) "normal vector1 old",magnvec,nx1,ny1,nz1,nx1*sx1+ny1*sy1+nz1*sz1,nx1*sx2+ny1*sy2+nz1*sz2
!
!
!!!!! Compute normal from analytical surface
!! 
!nx1 = 2.0_dbl*px1
!ny1 = 2.0_dbl*py1
!IF (k.GT.1) THEN
!	nz1 = (2.0_dbl*((r(k)-r(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*r(k-1)+(pz1-(z(k-1)/zcf))*r(k))/((z(k)-z(k-1))/zcf))/zcf
!ELSE 
!	nz1 = -(2.0_dbl*((r(k+1)-r(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*r(k)+(pz1-(z(k)/zcf))*r(k+1))/((z(k+1)-z(k))/zcf))/zcf
!ENDIF
!magnvec = SQRT(nx1**2+ny1**2+nz1**2)
!nx1 = nx1/magnvec
!ny1 = ny1/magnvec
!nz1 = nz1/magnvec
!!write(9,*) "normal vector1 analytical",magnvec,nx1,ny1,nz1,nx1*sx1+ny1*sy1+nz1*sz1,nx1*sx2+ny1*sy2+nz1*sz2
!!CALL FLUSH(9)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!! Set 1st tangent vector as the first surface vector
tx1 = sx1
ty1 = sy1
tz1 = sz1
!! Calculate 2nd tangent vector on tangent plane using surface normal direction and one of the surface vector
tx2 = (ny1*tz1-nz1*ty1)
ty2 = -(nx1*tz1-nz1*tx1)
tz2 = (nx1*ty1-ny1*tx1)
magtvec = SQRT(tx2**2+ty2**2+tz2**2)
tx2 = tx2/magtvec
ty2 = ty2/magtvec
tz2 = tz2/magtvec

!! Re-calculate 1st tangent vector on tangent plane using surface normal direction and the other surface vector
tx1 = -(ny1*tz2-nz1*ty2)
ty1 = -(-(nx1*tz2-nz1*tx2))
tz1 = -(nx1*ty2-ny1*tx2)
magtvec = SQRT(tx1**2+ty1**2+tz1**2)
tx1 = tx1/magtvec
ty1 = ty1/magtvec
tz1 = tz1/magtvec

!write(9,*) "tangent plane vectors", SQRT(tx1**2+ty1**2+tz1**2),sx1,sy1,sz1,tx1,ty1,tz1,sx2,sy2,sz2,tx2,ty2,tz2,SQRT(tx2**2+ty2**2+tz2**2)
!CALL FLUSH(9)
!write(9,*) "normal vector2",magnvec,nx1,ny1,nz1,nx1*tx1+ny1*ty1+nz1*tz1,nx1*tx2+ny1*ty2+nz1*tz2
!CALL FLUSH(9)

!write(*,*) "normal vector2",magnvec,nx1,ny1,nz1,nx1*tx1+ny1*ty1+nz1*tz1,nx1*tx2+ny1*ty2+nz1*tz2

!write(9,*) "normal vector error",nx1*tx1+ny1*ty1+nz1*tz1,nx1*tx2+ny1*ty2+nz1*tz2

!PAUSE

!
!! We have dphi/dn. Need to get dphi/dt from the fluid
!
!! get dphidtA
i = ip0
j = jp0
k = kp0
IF ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i-1,j,k))/(x(i+1)-x(i-1))
ELSEIF ((node(i+1,j,k).NE.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i,j,k)-phiTemp(i-1,j,k))/(x(i)-x(i-1))
ELSEIF  ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).NE.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i,j,k))/(x(i+1)-x(i))
ENDIF
IF ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j-1,k))/(y(j+1)-y(j-1))
ELSEIF ((node(i,j+1,k).NE.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j,k)-phiTemp(i,j-1,k))/(y(j)-y(j-1))
ELSEIF  ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).NE.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j,k))/(y(j+1)-y(j))
ENDIF
IF ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k-1))/(z(k+1)-z(k-1))
ELSEIF ((node(i,j,k+1).NE.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k)-phiTemp(i,j,k-1))/(z(k)-z(k-1))
ELSEIF  ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).NE.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k))/(z(k+1)-z(k))
ENDIF

dphidtA1 = dphidx*tx1 + dphidy*ty1 + dphidz*tz1
dphidtA2 = dphidx*tx2 + dphidy*ty2 + dphidz*tz2
!! get dphidtB 
!
i = ip1
j = jp1
k = kp1
IF ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i-1,j,k))/(x(i+1)-x(i-1))
ELSEIF ((node(i+1,j,k).NE.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i,j,k)-phiTemp(i-1,j,k))/(x(i)-x(i-1))
ELSEIF  ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).NE.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i,j,k))/(x(i+1)-x(i))
ENDIF
IF ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j-1,k))/(y(j+1)-y(j-1))
ELSEIF ((node(i,j+1,k).NE.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j,k)-phiTemp(i,j-1,k))/(y(j)-y(j-1))
ELSEIF  ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).NE.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j,k))/(y(j+1)-y(j))
ENDIF
IF ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k-1))/(z(k+1)-z(k-1))
ELSEIF ((node(i,j,k+1).NE.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k)-phiTemp(i,j,k-1))/(z(k)-z(k-1))
ELSEIF  ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).NE.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k))/(z(k+1)-z(k))
ENDIF

dphidtB1 = dphidx*tx1 + dphidy*ty1 + dphidz*tz1
dphidtB2 = dphidx*tx2 + dphidy*ty2 + dphidz*tz2
!! get dphidtC 
!
i = ip2
j = jp2
k = kp2
IF ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i-1,j,k))/(x(i+1)-x(i-1))
ELSEIF ((node(i+1,j,k).NE.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i,j,k)-phiTemp(i-1,j,k))/(x(i)-x(i-1))
ELSEIF  ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).NE.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i,j,k))/(x(i+1)-x(i))
ENDIF
IF ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j-1,k))/(y(j+1)-y(j-1))
ELSEIF ((node(i,j+1,k).NE.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j,k)-phiTemp(i,j-1,k))/(y(j)-y(j-1))
ELSEIF  ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).NE.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j,k))/(y(j+1)-y(j))
ENDIF
IF ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k-1))/(z(k+1)-z(k-1))
ELSEIF ((node(i,j,k+1).NE.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k)-phiTemp(i,j,k-1))/(z(k)-z(k-1))
ELSEIF  ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).NE.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k))/(z(k+1)-z(k))
ENDIF

dphidtC1 = dphidx*tx1 + dphidy*ty1 + dphidz*tz1
dphidtC2 = dphidx*tx2 + dphidy*ty2 + dphidz*tz2

!! get dphidt
dphidtAst1 = 0.5_dbl*(q+1.0_dbl)*(q+2.0_dbl)*dphidtA1+0.5_dbl*(q+1.0_dbl)*(q)*dphidtC1-(q+2.0_dbl)*(q)*dphidtB1
dphidtAst2 = 0.5_dbl*(q+1.0_dbl)*(q+2.0_dbl)*dphidtA2+0.5_dbl*(q+1.0_dbl)*(q)*dphidtC2-(q+2.0_dbl)*(q)*dphidtB2

! We already know dphidn. Now get dphidk.
! get grad(phi) vector
! Set up linear systems to get components of grad(phi) at Astar

!!write(*,*) nx1,ny1,nz1,tx1,ty1,tz1,tx2,ty2,tz2 
!arr(1,1) = nx1; 
!arr(1,2) = ny1; 
!arr(1,3) = nz1;
!arr(2,1) = tx1; 
!arr(2,2) = ty1; 
!arr(2,3) = tz1;
!arr(3,1) = tx2; 
!arr(3,2) = ty2;
!arr(3,3) = tz2;
!brhs(1) = dphidn;
!brhs(2) = dphidtAst1;
!brhs(3) = dphidtAst2;
!CALL Gauss(arr,arrinv,3)
!dphidx = arrinv(1,1)*brhs(1)+ arrinv(1,2)*brhs(2)+ arrinv(1,3)*brhs(3)
!dphidy = arrinv(2,1)*brhs(1)+ arrinv(2,2)*brhs(2)+ arrinv(2,3)*brhs(3)
!dphidz = arrinv(3,1)*brhs(1)+ arrinv(3,2)*brhs(2)+ arrinv(2,3)*brhs(3)

!write(*,*) nx1,ny1,nz1,tx1,ty1,tz1,tx2,ty2,tz2 
aa1 = nx1; 
aa2 = ny1; 
aa3 = nz1;
bb1 = tx1; 
bb2 = ty1; 
bb3 = tz1;
cc1 = tx2; 
cc2 = ty2;
cc3 = tz2;
rb1 = dphidn;
rb2 = dphidtAst1;
rb3 = dphidtAst2;

t1 = (cc1*aa2-cc2*aa1)*(cc1*bb3-cc3*bb1)-(cc1*aa3-cc3*aa1)*(cc1*bb2-cc2*bb1)
t2 = (cc1*rb2-cc2*rb1)*(cc1*bb3-cc3*bb1)-(cc1*rb3-cc3*rb1)*(cc1*bb2-cc2*bb1)
sx = t2/(MAX(abs(t1),epsmin)*sign(1.0_dbl,t1))
t1 = cc1*bb2-cc2*bb1
sy = ((cc1*rb2-cc2*rb1)-sx*(cc1*aa2-cc2*aa1))/(MAX(abs(t1),epsmin)*sign(1.0_dbl,t1))
t1 = cc1
sz = (rb1 - aa1*sx - bb1*sy)/(MAX(abs(t1),epsmin)*sign(1.0_dbl,t1))

dphidx = sx
dphidy = sy
dphidz = sz






! test - at the wall grad(phi).normalvec = dphidn
!write(*,*) dphidtA1,dphidtB1,dphidtC1,dphidtA2,dphidtB2,dphidtC2 
!write(*,*) dphidx,dphidy,dphidz,dphidn,dphidtAst1,dphidtAst2 

!temp = aa1*sx+aa2*sy+aa3*sz
!if (temp.gt.0.001) then
!write (9,*) temp,aa1,aa2,aa3,sx,sy,sz
!endif

!temp = dphidx*nx1+dphidy*ny1+dphidz*nz1
!If (abs(temp-dphidn).gt.0.001) then
!write(9,*) "Error in dphidn",temp,sx*nx1+sy*ny1+sz*nz1,ip0,jp0,kp0,m,nx1,ny1,nz1 
!write(9,*) "dphidx and sx",dphidx,dphidy,dphidz,sx,sy,sz
!endif


!PAUSE


!! get dphidk = gradphi.dot.k_vec
dist = SQRT(ex(bb(m))**2+ey(bb(m))**2+ez(bb(m))**2)
dx=ex(bb(m))/dist	! i direction
dy=ey(bb(m))/dist	! j direction
dz=ez(bb(m))/dist	! k direction
dphidk = (dphidx*dx+dphidy*dy+dphidz*dz)
!! now get phiAst
phiA = phiTemp(ip0,jp0,kp0)
phiB = phiTemp(ip1,jp1,kp1)
phiC = phiTemp(ip2,jp2,kp2)
delta = SQRT(ex(bb(m))**2 + ex(bb(m))**2 + ex(bb(m))**2)*xcf
phiAst = ((q+1.0_dbl)*(q+1.0_dbl)*phiA-q*q*phiB-dphidk*delta*q*(q+1.0_dbl))/(2.0_dbl*q+1.0_dbl)
!
!! return phiWall
!return phiAst

!CLOSE(9)
!phiAst = phiWall

!------------------------------------------------
END SUBROUTINE GetPhiWall
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE GetPerturbedVector(dx,dy,dz,epsx,epsy,epsz)	! Perturb a vector with components dx,dy and dz within a cone (using method from stackoverflow, Ray Tracing News)
!--------------------------------------------------------------------------------------------------
REAL(dbl), INTENT(IN) :: dx,dy,dz  ! direction vector components
REAL(dbl), INTENT(OUT) :: epsx,epsy,epsz ! deviation vector components
REAL(dbl) :: eps,epsmin  ! local variables
REAL(dbl) :: halfangle  ! cone halfangle
REAL(dbl) :: term1,term2,term3  ! temp variables
REAL(dbl) :: randx1,randy1,randz1 ! 1st random unit vector on tangent plane to (dx,dy,dz) vector
REAL(dbl) :: randx2,randy2,randz2 ! 2nd random unit vector on tangent plane to (dx,dy,dz) vector
REAL(dbl) :: s,rv,h,dist,theta,phiang,sinT,xc,yc,zc,newdx,newdy,newdz,randnum ! temp variables

halfangle = 3.0_dbl ! cone half angle in degrees
epsmin = 1.0e-8_dbl

! Algorithm to calculate perturbed vector within a cone

!dist = SQRT(dx**2+dy**2+dz**2)
!dx=dx/dist	! i direction
!dy=dy/dist	! j direction
!dz=dz/dist	! k direction


CALL RANDOM_NUMBER(randnum)
randx1 = 0.1_dbl*(0.0_dbl+randnum)+0.0_dbl
CALL RANDOM_NUMBER(randnum)
randy1 = 0.1_dbl*(0.0_dbl+randnum)+0.0_dbl
CALL RANDOM_NUMBER(randnum)
randz1 = 0.1_dbl*(0.0_dbl+randnum)+0.0_dbl
!write(*,*) randx1,randy1,randz1
if (abs(dx).GT.1.0e-8_dbl) then
randx1 = -(dz*randz1+dy*randy1)/(dx+epsmin)
elseif (abs(dy).GT.1.0e-8_dbl) then
randy1 = -(dz*randz1+dx*randx1)/(dy+epsmin)
else
randz1 = -(dx*randx1+dy*randy1)/(dz+epsmin)
endif
!randx1 = -(dz*randz1+dy*randy1)/(dx+epsmin)
!randy1 = -(dz*randz1+dx*randx1)/(dy+epsmin)
!randz1 = -(dx*randx1+dy*randy1)/(dz+epsmin)
dist = SQRT(randx1**2+randy1**2+randz1**2)
randx1=randx1/dist	! i direction
randy1=randy1/dist	! j direction
randz1=randz1/dist	! k direction

randx2 = (dy*randz1-dz*randy1)
randy2 = -(dx*randz1-dz*randx1)
randz2 = (dx*randy1-dy*randx1)

dist = SQRT(randx2**2+randy2**2+randz2**2)
randx2=randx2/dist	! i direction
randy2=randy2/dist	! j direction
randz2=randz2/dist	! k direction
!write(9,*) "testting",dx*randx1+dy*randy1+dz*randz1,dx*randx2+dy*randy2+dz*randz2 &
!		,randz1,dz,randx1,dx,randy1,dy
!CALL FLUSH(9)

!s = rand(0)
!r = rand(0)
CALL RANDOM_NUMBER(s)
CALL RANDOM_NUMBER(rv)
s = s*1.0_dbl + 0.0_dbl
rv = rv + 0.0_dbl
theta = halfangle*PI/(180.0_dbl) ! create a cone with coning angle 10 degrees
h = cos(theta)
phiang = 2.0_dbl*PI*s
zc = h + ( 1.0_dbl - h ) * rv
sinT = sqrt( 1.0_dbl - zc * zc )
xc = cos( phiang ) * sinT
yc = sin( phiang ) * sinT
! Estimate the perturbed vector
!perturbed = rand_vector * xc + cross_vector * yc + orig_vector * zc
newdx = randx1*xc + randx2*yc + dx*zc
newdy = randy1*xc + randy2*yc + dy*zc
newdz = randz1*xc + randz2*yc + dz*zc

dist = SQRT(newdx**2+newdy**2+newdz**2)
newdx = newdx/dist	! i direction
newdy = newdy/dist	! j direction
newdz = newdz/dist	! k direction


! We return the perturbations instead of the perturbed vector
epsx = newdx - dx
epsy = newdy - dy
epsz = newdz - dz
!write(9,*) "original vector",dx,dy,dz,dist
!write(9,*) "perturbed vector",newdx,newdy,newdz,acos(newdx*dx+newdy*dy+newdz*dz)*180.0/PI


!eps = 0.1_dbl
!epsx = 0.1_dbl
!term1 = 1.0_dbl + ((dy**2)/(dz**2+epsmin))
!term2 = 2.0_dbl*(dx*dy)/(dz**2+epsmin)
!term3 = (epsx**2)*(1.0_dbl + ((dx**2)/(dz**2+epsmin)))-(eps**2) 
!epsy = (-term2 + SQRT(term2**2 - 4.0_dbl*term1*term3))/(2.0_dbl*term1)
!epsz = -(dx*epsx+dy*epsy)/(dz*(1.0_dbl+epsmin))
!epsx=0.01_dbl
!epsy=0.01_dbl
!epsz=0.01_dbl


!------------------------------------------------
END SUBROUTINE GetPerturbedVector
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE RotateVector(dx,dy,dz,ang1,ang2,ang3)	! Rotate a vector with components dx,dy and dz by angles ang1,ang2,ang3 (in radians)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

REAL(dbl), INTENT(INOUT) :: dx,dy,dz	! original vector
REAL(dbl) :: tempdx,tempdy,tempdz	! teporary vector components
REAL(dbl), INTENT(IN) :: ang1,ang2,ang3	! rotation angles to rotate (perturb) the vector joining fluid and solid nodes for identifying nearby nodes on the solid surface

! ROtate a vector using a 3D ortation matrix. I directly compute the matrix multiplied vector

tempdx = dx*(cos(ang2)*cos(ang1)) + dy*(cos(ang3)*sin(ang1)+sin(ang3)*sin(ang2)*cos(ang1))  &
		+ dz*(sin(ang3)*sin(ang1)-cos(ang3)*sin(ang2)*cos(ang1))
tempdy = dx*(-cos(ang2)*sin(ang1)) + dy*(cos(ang3)*cos(ang1)-sin(ang3)*sin(ang2)*sin(ang1)) &
		+ dz*(sin(ang3)*cos(ang1)+cos(ang3)*sin(ang2)*sin(ang1))
tempdz = dx*(sin(ang2)) + dy*(-sin(ang3)*cos(ang2)) + dz*(cos(ang3)*cos(ang2))

dx=tempdx
dy=tempdy
dz=tempdz

!------------------------------------------------
END SUBROUTINE RotateVector
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE RotateVectorAboutAxis(dx,dy,dz,ang1,ang2,ang3)	! Rotate a vector with components dx,dy and dz by angles ang1,ang2,ang3 (in radians)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

REAL(dbl), INTENT(INOUT) :: dx,dy,dz	! original vector
REAL(dbl) :: tempdx,tempdy,tempdz	! teporary vector components
REAL(dbl), INTENT(IN) :: ang1,ang2,ang3	! rotation angles to rotate (perturb) the vector joining fluid and solid nodes for identifying nearby nodes on the solid surface

! ROtate a vector using a 3D ortation matrix. I directly compute the matrix multiplied vector

tempdx = dx*(cos(ang1)+dx*dx*(1.0_dbl-cos(ang1))) + dy*(dx*dy*(1.0_dbl-cos(ang1))-dz*sin(ang1))  &
		+ dz*(dx*dz*(1.0_dbl-cos(ang1))+dy*sin(ang1))
tempdy = dx*(dy*dx*(1.0_dbl-cos(ang1))+dz*sin(ang1)) + dy*(cos(ang1)+dy*dy*(1.0_dbl-cos(ang1))) &
		+ dz*(dy*dz*(1.0_dbl-cos(ang1))-dx*sin(ang1))
tempdz = dx*(dz*dx*(1.0_dbl-cos(ang1))-dy*sin(ang1)) + dy*(dz*dy*(1.0_dbl-cos(ang1))+dx*sin(ang1)) &
		+ dz*(cos(ang1)+dz*dz*(1.0_dbl-cos(ang1)))

dx=tempdx
dy=tempdy
dz=tempdz

!------------------------------------------------
END SUBROUTINE RotateVectorAboutAxis
!------------------------------------------------


!--------------------------------------------------------------------------------------------------
SUBROUTINE IntersectSurface(m,i,j,k,im1,jm1,km1,px,py,pz,q,ang1,ang2,ang3)			! Intersects boudnary surface using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
REAL(dbl), INTENT(INOUT) :: ang1,ang2,ang3	! rotation angles to rotate (perturb) the vector joining fluid and solid nodes for identifying nearby nodes on the solid surface
REAL(dbl), INTENT(OUT) :: px,py,pz			! boundary points
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Bx,By,Bz									! solid node
REAL(dbl) :: AB,AP1,AP2,AP										! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: dxorig,dyorig,dzorig				! unit vector pointing from A to B
REAL(dbl) :: r1,r2,z1,z2,slope,intercept			! radius and z location at k and km1, slope of line connecting those two points, z-intercept of the r-equation
REAL(dbl) :: slope2,term1,term2						! terms used in calculation
REAL(dbl), INTENT(OUT) :: q						! distance from fluid point to surface using Ray Tracing


!open(9,file='testoutput.dat',position='append')

! RAY
! point A (current node)
Ax = x(i)/xcf
Ay = y(j)/ycf
Az = z(k)/zcf

! point B (solid node)
Bx = x(im1)/xcf
By = y(jm1)/ycf
Bz = z(km1)/zcf

!! point A (current node)
!Ax = x(i)
!Ay = y(j)
!Az = z(k)
!
!! point B (solid node)
!Bx = x(im1)
!By = y(jm1)
!Bz = z(km1)

!! distance from A to B
!AB = SQRT((Bx - Ax)*(Bx - Ax) + (By - Ay)*(By - Ay) + (Bz - Az)*(Bz - Az))
!
!! unit vector (d) from point A to point B
!dx = (x(im1)-x(i))/AB									! i direction
!dy = (y(jm1)-y(j))/AB									! j direction
!dz = (z(km1)-z(k))/AB									! k direction

!!! unit vector (d) from point A to point B
!AB = SQRT(ex(m)**2+ey(m)**2+ez(m)**2)
!dxorig=ex(bb(m))/AB	! i direction
!dyorig=ey(bb(m))/AB	! j direction
!dzorig=ez(bb(m))/AB	! k direction

! distance from A to B
AB = SQRT((Bx - Ax)*(Bx - Ax) + (By - Ay)*(By - Ay) + (Bz - Az)*(Bz - Az))

! unit vector (d) from point A to point B
dxorig = (x(im1)-x(i))/AB									! i direction
dyorig = (y(jm1)-y(j))/AB									! j direction
dzorig = (z(km1)-z(k))/AB									! k direction

!Write(9,*) "test1",dx,dy,dz
!CALL FLUSH(9)
! Rotate vector (d) by (ang1,ang2,ang3) using a transformation matrix
!CALL RotateVector(dx,dy,dz,ang1,ang2,ang3)
!CALL RotateVectorAboutAxis(dx,dy,dz,ang1,ang2,ang3)
!Write(9,*) "test2",dx,dy,dz

10 dx=dxorig+ang3
dy=dyorig+ang2
dz=dzorig+ang1

AB = SQRT(dx**2+dy**2+dz**2)
dx=dx/AB	! i direction
dy=dy/AB	! j direction
dz=dz/AB	! k direction

!Write(9,*) "test2",dx,dy,dz
!CALL FLUSH(9)
! distance from A to B
AB = SQRT((Bx - Ax)*(Bx - Ax) + (By - Ay)*(By - Ay) + (Bz - Az)*(Bz - Az))

! SURFACE
IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	r1 = rOut(k)													! radius at k (distance from CL)
	r2 = rOut(km1)													! radius at km1 (distance from CL)
ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
	r1 = rIn(k)													! radius at k (distance from CL)
	r2 = rIn(km1)													! radius at km1 (distance from CL)
END IF
z1 = z(k)/zcf													! z-coordinate at k
z2 = z(km1)/zcf													! z-coordinate at km1

IF(k .NE. km1) THEN
  slope = (r2-r1)/(z2-z1)								! approximate the surface as a conincal shell (linear between k values)
ELSE
  slope = 0.0_dbl
END IF

intercept = r1 - slope*z1								! z-intercept of the linearly approximated r-equation

! terms used in calculation
slope2 = slope*slope															! slope^2
term1 = Ax*dx + Ay*dy - Az*dz*slope2 - intercept*dz*slope	! reoccuring term
term2 = dx*dx + dy*dy - dz*dz*slope*slope			! reoccuring term

! calculate the distance from the current node (point A) to the wall (point P)
AP1 = (1.0_dbl/(2.0_dbl*term2)) * &
     (-2.0_dbl*term1					&
   + SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))
AP2 = (1.0_dbl/(2.0_dbl*term2)) * &
     (-2.0_dbl*term1					&
   - SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))
IF ((AP1.GE.0.0_dbl).AND.(AP2.GE.0.0_dbl)) THEN
 AP = MIN(AP1,AP2)
ELSE
 AP = MAX(AP1,AP2)
ENDIF
q = AP/AB	! distance ratio

!IF ((q.GT.2.0_dbl).AND.(abs(dxorig).LE.0.0000001_dbl).AND.(abs(dyorig).LE.0.0000001_dbl)) THEN
!	ang3 = -ang3
!	ang2 = -ang2
!	GOTO 10
!ENDIF

!AP = (1.0_dbl/(2.0_dbl*term2)) * &
!     (-2.0_dbl*term1					&
!   + SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))
!q = AP/AB	! distance ratio

! calculate coordinates of the boudnayr interseciton point
px = Ax + AP*dx
py = Ay + AP*dy
pz = Az + AP*dz


! balaji added
!q=0.5

! make sure 0<q<1
IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 2.00000001_dbl)) THEN 
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) " ERROR in IntersectSurface"
  WRITE(1000,*) "q=",q
  WRITE(1000,*) "m=",m
  WRITE(1000,*) "i=",i,"j=",j,"k=",k
  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
  WRITE(1000,*) "Ax=",Ax,"Ay=",Ay,"Az=",Az
  WRITE(1000,*) "Bx=",Bx,"By=",By,"Bz=",Bz
  WRITE(1000,*) "Px=",px,"Py=",py,"Pz=",pz
  WRITE(1000,*) "dxorig=",dxorig,"dyorig=",dyorig,"dzorig=",dzorig
  WRITE(1000,*) "dx=",dx,"dy=",dy,"dz=",dz
  WRITE(1000,*) "epsx=",ang3,"epsy=",ang2,"epsz=",ang1
  WRITE(1000,*) "r1=",r1,"r2=",r2
  WRITE(1000,*) "z1=",z1,"z2=",z2
  WRITE(1000,*) "slope=",slope
  WRITE(1000,*) "term1=",term1,"term2=",term2
  WRITE(1000,*) "intercept=",intercept
  WRITE(1000,*) "AB=",AB,"AP=",AP
  CLOSE(1000)
  STOP
END IF																																
!CLOSE(9)

!------------------------------------------------
END SUBROUTINE IntersectSurface
!------------------------------------------------

! -------------------------------------------------------------------- 
SUBROUTINE Gauss (a,ainv,n)       ! Invert matrix by Gauss method 
! -------------------------------------------------------------------- 
IMPLICIT NONE 
INTEGER :: n 
REAL(dbl), INTENT(IN) :: a(n,n) 
REAL(dbl), INTENT(OUT) :: ainv(n,n) 
! - - - Local Variables - - - 
REAL(dbl) :: b(n,n), c, d, temp(n) 
INTEGER(lng) :: i, j, k, m, imx(1), ipvt(n) 
! - - - - - - - - - - - - - - 
b = a 
ipvt = (/ (i, i = 1, n) /) 
DO k = 1,n 
   imx = MAXLOC(ABS(b(k:n,k))) 
   m = k-1+imx(1) 
   IF (m /= k) THEN 
      ipvt( (/m,k/) ) = ipvt( (/k,m/) ) 
      b((/m,k/),:) = b((/k,m/),:) 
   END IF 
   d = 1/b(k,k) 
   temp = b(:,k) 
   DO j = 1, n 
      c = b(k,j)*d 
      b(:,j) = b(:,j)-temp*c 
      b(k,j) = c 
   END DO 
   b(:,k) = temp*(-d) 
   b(k,k) = d 
END DO 
ainv(:,ipvt) = b 
!--------------------------------------------------------------------------------------------------
END SUBROUTINE Gauss
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarBCV(m,i,j,k,im1,jm1,km1,vNum,phiBC)						! implements the scalar BCs for the villi
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1								! index variables
INTEGER(lng), INTENT(IN) :: vNum													! number of the current villus
REAL(dbl), INTENT(OUT) :: phiBC     											! scalar contribution from the boundary condition
INTEGER(lng) :: ip1,jp1,kp1 														! neighboring nodes (2 away from the wall)
REAL(dbl) :: Cx,Cy,Cz,Vx,Vy,Vz													! vector between villous base and current node, vector between villous base and villous tip
REAL(dbl) :: uV,vV,wV																! villus wall velocity at actual coordinate system
REAL(dbl) :: q																			! distance ratio from the current node to the solid node
REAL(dbl) :: rhoB,phiB																! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: feq_m																	! equilibrium distribution function in the mth direction
REAL(dbl) :: phiijk_m																! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)

! find C and V
CALL qCalcV(m,i,j,k,im1,jm1,km1,vNum,q,Cx,Cy,Cz,Vx,Vy,Vz,2)

! find the influence of villous velocity on the current point
CALL VilliVelocity(vNum,Cx,Cy,Cz,uV,vV,wV)

! neighboring node (fluid side)	
ip1 = i + ex(m) 																		! i + 1
jp1 = j + ey(m)																		! j + 1
kp1 = k + ez(m)																		! k + 1		

! if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF	

! assign values to boundary (density, scalar, f)
rhoB = (rho(i,j,k) - rho(ip1,jp1,kp1))*(1+q) + rho(ip1,jp1,kp1)		! extrapolate the density
CALL Equilibrium_LOCAL(m,rhoB,uV,vV,wV,feq_m)								! calculate the equibrium distribution function in the mth direction

! find the contribution of scalar streamed from the wall to the current node (i,j,k), and from the current node to the next neighboring node (ip1,jp1,kp1)
phiB		= (feq_m/rhoB - wt(m)*Delta)*phiWall								! contribution from the wall in the mth direction (zero if phiWall=0)
phiijk_m	= (fplus(m,i,j,k)/rho(i,j,k) - wt(m)*Delta)*phiTemp(i,j,k)	! contribution from the current node to the next node in the mth direction

! if q is too small, the extrapolation to phiBC can create a large error...
IF(q .LT. 0.25) THEN
  q = 0.25_dbl  																		! approximate the distance ratio as 0.25
END IF

! extrapolate using phiB and phijk_m to obtain contribution from the solid node to the current node
phiBC		= ((phiB - phiijk_m)/q) + phiB										! extrapolated scalar value at the solid node, using q

!------------------------------------------------
END SUBROUTINE ScalarBCV
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Equilibrium_LOCAL(m,rhoijk,uijk,vijk,wijk,feq_m)		! calculate and store the equilibrium distribution function
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN)	:: m											! distribution function direction		
REAL(dbl), INTENT(IN)		:: rhoijk,uijk,vijk,wijk				! density, and velocity components at the current node
REAL(dbl), INTENT(OUT)		:: feq_m										! density, and velocity components at the current node

REAL(dbl)						:: UU,ue,ve,we,Usum						! precalculated quantities for use in the feq equation

UU 	= uijk*uijk + vijk*vijk + wijk*wijk								! U . U
      
ue		= uijk*ex(m)															! u . e
ve		= vijk*ey(m)															! v . e
we		= wijk*ez(m)															! w . e

Usum	= ue + ve + we															! U . e
        
feq_m	= (wt(m)*rhoijk)*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*UU)	! equilibrium distribution function in the mth direction
        
!------------------------------------------------
END SUBROUTINE Equilibrium_LOCAL
!------------------------------------------------

!================================================
END MODULE ICBC
!================================================
