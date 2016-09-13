!==================================================================================================
MODULE IC		! Sets Initial  Conditions
!==================================================================================================
USE SetPrecision
USE Setup  
USE MPI

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE ICs	! sets the initial conditions
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE  

INTEGER(lng) :: i, j, k, m, mpierr					! index variables
REAL(dbl)    :: feq
INTEGER(lng) :: imintemp, imaxtemp				! index variables
REAL(dbl)    :: x_center
REAL(dbl)    :: alpha,xmin, xmax,xmid
REAL(dbl)    :: Drug_Initial_l, Drug_Released_Total_l
REAL(dbl)    :: Drug_Absorbed_restart_l, Drug_Remained_in_Domain_l 
CHARACTER(7) :: iter_char                        ! iteration stored as a character

IF (Flag_Restart) THEN                           ! restart from  file 
 
   OPEN(55,FILE='Restart-iter.dat')              ! open initial iteration file
   READ(55,*) iter0                              ! read and set initial iteration
   CLOSE(55)

   WRITE(iter_char(1:7),'(I7.7)') iter0
   OPEN(50,FILE='Restart-Out-'//iter_char//'-'//sub//'.dat')
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
   READ(50,*) Drug_Initial_l
   READ(50,*) Drug_Released_Total
   READ(50,*) Drug_Absorbed_restart
   READ(50,*) Drug_Remained_in_Domain_l
   delphi_particle = 0.0_dbl                               	! Initialize the scalar contirbution from particles to 0.0. Once the particle
   CLOSE(50)
   CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   CALL MPI_ALLREDUCE(Drug_Initial_l, 		Drug_Initial,		1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
  !CALL MPI_ALLREDUCE(Drug_Released_Total_l, 	Drug_Released_Total,	1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
   CALL MPI_ALLREDUCE(Drug_Remained_in_Domain_l,Drug_Remained_in_Domain,1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

   iter0 = iter0 + 1						! set the initial iteration to one after last iteration from the previous run

ELSE								! clean start
  alpha = 0.4_dbl
  imintemp = -ANINT(alpha*(nx-1_lng))+ANINT(0.5_dbl*(nx+1))
  imaxtemp = +ANINT(alpha*(nx-1_lng))+ANINT(0.5_dbl*(nx+1))
  xmin = 0.5_dbl*(xx(imintemp)+xx(imintemp-1))
  xmax = 0.5_dbl*(xx(imaxtemp)+xx(imaxtemp+1))
  xmid = 0.5_dbl*(xmax+xmin)

  !----- Initial conditions on velocity, density, and scalar
  DO k=0,nzSub+1_lng
     DO j=0,nySub+1_lng
        DO i=0,nxSub+1_lng
           u(i,j,k)= 0.0_dbl 
           v(i,j,k)= 0.0_dbl
           IF (Flag_Couette) THEN
              w(i,j,k)= (xx(i+iMin-1) / (0.40_dbl*D)) * (s1/vcf)
           ELSE
              w(i,j,k)  = 0.0_dbl                                       ! z-velocity
           END IF   
           rho(i,j,k)= denL                                          ! density
           !------ distribution functions (set to equilibrium)
           DO m=0,NumDistDirs
              CALL Equilibrium_LOCAL(m,rho(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),feq)	
              f(m,i,j,k) = feq
           END DO
        END DO
     END DO
  END DO

  !----- Starting iteration
  iter0 = 1_lng

  !----- Initialize scalar values
  Drug_Absorbed_restart= 0.0_dbl
  phiAbsorbed    = 0.0_dbl								! total amount of scalar absorbed
  phiAbsorbedS   = 0.0_dbl								! total amount of scalar absorbed through the macroscopic surface
  phiAbsorbedV   = 0.0_dbl								! total amount of scalar absorbed through the villi
  phiInOut       = 0.0_dbl								! total amount of scalar leaving the inlet/outlet
  delphi_particle= 0.0_dbl								! Initialize the scalar contirbution from particles to 0.0. Once the particle
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
INTEGER(lng)   :: parid_r
REAL(dbl)      :: xp_r, yp_r, zp_r, up_r, vp_r, wp_r, rp_r, xpold_r, ypold_r, zpold_r			! parameters read form particle_restart.dat
REAL(dbl)      :: upold_r, vpold_r, wpold_r, rpold_r, par_conc_r, gamma_cont_r, Sh_r, S_r		! parameters read form particle_restart.dat
REAL(dbl)      :: Sst_r, Veff_r, Nbj_r, bulk_conc_r, delNBbyCV_r, cur_part_r, new_part_r		! parameters read form particle_restart.dat
REAL(dbl) :: xp,yp,zp,par_radius
CHARACTER(7) :: iter_char                                       ! iteration stored as a character
TYPE(ParRecord), POINTER	:: CurPar

IF (Flag_Restart) THEN								! restarting: read particle data from  particle_restart.dat
   OPEN(55,FILE='Restart-iter.dat')                                    ! open initial iteration file
   READ(55,*) iter0                                             ! read and set initial iteration
   CLOSE(55)

   WRITE(iter_char(1:7),'(I7.7)') iter0
   OPEN(59,FILE='Restart-Particles-'//iter_char//'.dat')
   iter0= iter0 + 1

   READ(59,*) np
   num_particles = np
   CALL list_init(ParListHead)
   CurPar => ParListHead

   DO i = 1, np
      READ(59,*) parid_r,xp_r,yp_r,zp_r,up_r,vp_r,wp_r,rp_r,xpold_r,ypold_r,zpold_r, &			
	         upold_r,vpold_r,wpold_r,rpold_r,par_conc_r,gamma_cont_r,Sh_r,S_r,   &
                 Sst_r,Veff_r,Nbj_r,bulk_conc_r,delNBbyCV_r,cur_part_r,new_part_r                                
      CALL list_init(CurPar%next)      
      CurPar%next%prev => CurPar
      CurPar%next%next => null()
      CurPar%next%pardata%parid 	= 	parid_r
      CurPar%next%pardata%xp 		=	xp_r
      CurPar%next%pardata%yp 		=	yp_r
      CurPar%next%pardata%zp 		= 	zp_r
      CurPar%next%pardata%up 		=	up_r
      CurPar%next%pardata%vp 		=	vp_r
      CurPar%next%pardata%wp    	=	wp_r
      CurPar%next%pardata%rp    	=	rp_r
      CurPar%next%pardata%xpold 	= 	xpold_r
      CurPar%next%pardata%ypold 	= 	ypold_r
      CurPar%next%pardata%zpold 	= 	zpold_r
      CurPar%next%pardata%upold 	= 	upold_r
      CurPar%next%pardata%vpold 	= 	vpold_r
      CurPar%next%pardata%wpold 	= 	wpold_r
      CurPar%next%pardata%rpold 	= 	rpold_r
      CurPar%next%pardata%par_conc 	= 	par_conc_r
      CurPar%next%pardata%gamma_cont 	= 	gamma_cont_r
      CurPar%next%pardata%sh 		=	Sh_r
      CurPar%next%pardata%S 		=	S_r
      CurPar%next%pardata%Sst 		= 	Sst_r
      CurPar%next%pardata%Veff 		= 	Veff_r
      CurPar%next%pardata%Nbj 		= 	Nbj_r
      CurPar%next%pardata%bulk_conc 	= 	bulk_conc_r 
      CurPar%next%pardata%delNBbyCV	= 	delNBbyCV_r
      CurPar%next%pardata%cur_part	= 	cur_part_r
      CurPar%next%pardata%new_part	= 	new_part_r
      CurPar => CurPar%next
   END DO
   CLOSE(59)
ELSE										! starting from scratch. No data to be read from restart files
   OPEN(60,FILE='particle.dat')
   READ(60,*) np
   num_particles = np

   CALL list_init(ParListHead)
   CurPar => ParListHead
   DO i = 1, np
      READ(60,*) parid,xp,yp,zp,par_radius					! read particle.dat file
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
   END DO
   CLOSE(60)

ENDIF
!------------------------------------------------
END SUBROUTINE IniParticles
!------------------------------------------------








!--------------------------------------------------------------------------------------------------
SUBROUTINE IC_Drug_Distribution		! Sets/Maintains initial distributions of scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,ii,jj		! lattice indices

IF(iter .EQ. phiStart) THEN
  SELECT CASE(sclrIC) 
      
    CASE(BLOB)							! blob of scalar at the center of the domain
      DO k=0,nzSub+1
        DO j=0,nySub+1
          DO i=0,nxSub+1
!            phi(i,j,k) = phiIC*ee**(-((x(i)**2 + y(j)**2 + (z(k)-0.5_dbl*L)**2)/(2.0_dbl*sigma**2)))		! 3D Gaussian Distribution
             IF ((i.LT.32) .AND. (i.GT.24) .AND. (j.LT.32) .AND. (j.GT.24) .AND. (k.LE.30) .AND. (k.GE.0) ) THEN
                phi(i,j,k) = Cs_mol
             END IF
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


ELSE  								  ! MAINTAINENCE OF SCALAR

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
END SUBROUTINE IC_Drug_Distribution
!------------------------------------------------








!--------------------------------------------------------------------------------------------------
SUBROUTINE Equilibrium_LOCAL(m,rhoijk,uijk,vijk,wijk,feq_m)		! calculate and store the equilibrium distribution function
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN)	:: m					! distribution function direction		
REAL(dbl), INTENT(IN)		:: rhoijk,uijk,vijk,wijk		! density, and velocity components at the current node
REAL(dbl), INTENT(OUT)		:: feq_m				! density, and velocity components at the current node
REAL(dbl)			:: UU,ue,ve,we,Usum			! precalculated quantities for use in the feq equation

UU 	= uijk*uijk + vijk*vijk + wijk*wijk				! U . U
ue	= uijk*ex(m)							! u . e
ve	= vijk*ey(m)							! v . e
we	= wijk*ez(m)							! w . e
Usum	= ue + ve + we							! U . e
        
feq_m	= (wt(m)*rhoijk)*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*UU)	! equilibrium distribution function in the mth direction
        
!------------------------------------------------
END SUBROUTINE Equilibrium_LOCAL
!------------------------------------------------



!================================================
END MODULE IC
!================================================
