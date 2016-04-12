!==================================================================================================
MODULE IC		! Sets Initial and Boundary Conditions
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
	ELSE
		w(i,j,k) = 0.0_dbl
	ENDIF
        rho(i,j,k) = denL
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

            IF (((i.GE.16).AND.(i.LE.25)).AND.((j.GE.16).AND.(j.LE.25)).AND.((k.GE.16).AND.(k.LE.25))) THEN
	       	phi(i,j,k) = phiIC!i+j+k!phiIC		! 3D Cube test
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
END MODULE IC
!================================================
