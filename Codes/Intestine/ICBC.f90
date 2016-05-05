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
        u(i,j,k)   = 0.0_dbl							! x-velocity
        v(i,j,k)   = 0.0_dbl							! y-velocity
        w(i,j,k)   = 0.0_dbl							! z-velocity
        rho(i,j,k) = denL							! density
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

  ! Initialize scalar values
  phiAbsorbed	= 0.0_dbl								! total amount of scalar absorbed
  phiAbsorbedS	= 0.0_dbl								! total amount of scalar absorbed through the macroscopic surface
  phiAbsorbedV	= 0.0_dbl								! total amount of scalar absorbed through the villi
  phiInOut	= 0.0_dbl								! total amount of scalar leaving the inlet/outlet
  delphi_particle = 0.0_dbl								! Initialize the scalar contirbution from particles to 0.0. Once the particle
END IF

!------------------------------------------------
END SUBROUTINE ICs
!------------------------------------------------



!------------------------------------------------
SUBROUTINE IniParticles
!-----------------------------------------------
IMPLICIT NONE
INTEGER(lng)   :: i, parid,particle_partition,ipartition
INTEGER(lng)   :: mpierr
REAL(dbl) :: xp,yp,zp,par_radius
TYPE(ParRecord), POINTER	:: CurPar
IF (restart) THEN
	! Read particle number and position along with it's radius,concentration.
	! Interpolate to calculate particle velocities.
	! Interpolate particle concentration to nodes into delphi_particle.

ELSE
	!----- Linked list approach
	OPEN(60,FILE='particle.dat')
	READ(60,*) np
	num_particles = np

	!----- Initialize Header Pointer
	CALL list_init(ParListHead)
	CurPar => ParListHead

	DO i = 1, np
	   READ(60,*) parid,xp,yp,zp,par_radius		! read particle.dat file

	   !----- Search the partition this particle belongs to
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
 
         !----- Create a particle element in the linked list only if the particles belongs to this partition
!        IF (particle_partition.EQ.mySub) THEN
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
			CurPar%next%pardata%rp = par_radius
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
			CurPar%next%pardata%cur_part= particle_partition
			CurPar%next%pardata%new_part= particle_partition
	   CurPar => CurPar%next
!	END IF
     END DO
	
 CLOSE(60)

ENDIF
!------------------------------------------------
END SUBROUTINE IniParticles
!------------------------------------------------





!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarDistribution		! Sets/Maintains initial distributions of scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,ii,jj		! lattice indices

IF(iter .EQ. phiStart) THEN
  SELECT CASE(sclrIC) 
      
    CASE(BLOB)							! blob of scalar at the center of the domain
      
      DO k=0,nzSub+1
        DO j=0,nySub+1
          DO i=0,nxSub+1
            phi(i,j,k) = phiIC*ee**(-((x(i)**2 + y(j)**2 + (z(k)-0.5_dbl*L)**2)/(2.0_dbl*sigma**2)))		! 3D Gaussian Distribution
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

    CASE(LINE) 						! line of scalar along axis
      IF((SubID(2) .EQ. 0) .AND. (SubID(4) .EQ. 0)) THEN			! if no neighboring subdomains exist in the 2nd and 4th directions, then they lie at the centerline
        DO k=0,nzSub+1
          DO j=0,nySub+1
            DO i=0,nxSub+1
              IF((ABS(x(i)) .LE. 2.51_dbl*xcf) .AND. (ABS(y(j)) .LE. 2.51_dbl*ycf)) THEN				! (5.01 in case it is slightly higher than 5 due to round-off)
                phi(i,j,k) = phiIC											! 2D Gaussian Distribution in x and y (maintain phi=1 at r=2.5*zcf)
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

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1						! index variables
REAL(dbl), INTENT(OUT) :: fbb								! bounced back distribution function
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb									! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk									! radius of current node

rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))						! radius at current location

cosTheta = x(im1)/rijk									! COS(theta)
sinTheta = y(jm1)/rijk									! SIN(theta)

ub = vel(km1)*cosTheta									! x-component of the velocity at i,j,k
vb = vel(km1)*sinTheta									! y-component of the velocity at i,j,k
wb = 0.0_dbl										! no z-component in this case			

fbb = fplus(bb(m),i,j,k) + 6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m))	! bounced back distribution function with added momentum

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

  rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))						! radius at current location

  cosTheta = x(im1)/rijk										! COS(theta)
  sinTheta = y(jm1)/rijk										! SIN(theta)

  ub = vel(km1)*cosTheta										! x-component of the velocity at i,j,k
  vb = vel(km1)*sinTheta										! y-component of the velocity at i,j,k
  wb = 0.0_dbl														! no z-component in this case)

  CALL qCalc(m,i,j,k,im1,jm1,km1,q)							! calculate q					

  ! bounced back distribution function with added momentum
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
    fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 															&
        + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) 													& 
        - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 													&
        + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)) ! Set rho = 1.0

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
    fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q																&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)							&
        + (6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))	! Set rho = 1.0

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
		 ! Initial fluid node guess
                 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                
		 ! Initial solid node guess
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
                 
	 IF (k.NE.km1) THEN
                 DO it=1,10
		   ! guess of boundary location 
                   xt=(x1+x2)/2.0_dbl
                   yt=(y1+y2)/2.0_dbl
                   zt=(z1+z2)/2.0_dbl

      		   rt = SQRT(xt*xt + yt*yt)
		   ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))

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
	 ELSE
		  DO it=1,10
		   ! guess of boundary location 
                   xt=(x1+x2)/2.0_dbl
                   yt=(y1+y2)/2.0_dbl
                   zt=(z1+z2)/2.0_dbl

      		   rt = SQRT(xt*xt + yt*yt)
		   ht = (r(km1)+r(k))/2.0_dbl

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
	 ENDIF
		 cosTheta=xt/rt
		 sinTheta=yt/rt
	 IF (k.NE.km1) THEN
		 vt = ((zt-z(k))*vel(km1)+(z(km1)-zt)*vel(k))/(z(km1)-z(k))
	 ELSE
		 vt = (vel(k)+vel(km1))*0.5_dbl
	 ENDIF
		 ub = vt*cosTheta										! x-component of the velocity at i,j,k
		 vb = vt*sinTheta										! y-component of the velocity at i,j,k
		 wb = 0.0_dbl											! no z-component in this case)
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
        + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)) ! use rho = 1.0

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
    fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q																&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)							&
        + (6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))		! Use rho = 1.0

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
r1 = r(k)													! radius at k (distance from CL)
r2 = r(km1)													! radius at km1 (distance from CL)
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

cosTheta = x(im1)/r(km1)	! COS(theta)
sinTheta = y(jm1)/r(km1)	! SIN(theta)

ub = vel(km1)*cosTheta																! x-component of the velocity at i,j,k
vb = vel(km1)*sinTheta																! y-component of the velocity at i,j,k
wb = 0.0_dbl	

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
