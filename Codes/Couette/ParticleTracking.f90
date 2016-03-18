!==================================================================================================
MODULE ParticleTracking				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================
USE SetPrecision
USE Setup
USE ICBC
USE MPI
USE ParticleDrug

IMPLICIT NONE

CONTAINS

!===================================================================================================
SUBROUTINE Particle_Setup
!===================================================================================================
IMPLICIT NONE

IF (restart) THEN
ELSE
	CALL Interp_Parvel
ENDIF
!===================================================================================================
END SUBROUTINE Particle_Setup
!===================================================================================================



!===================================================================================================
SUBROUTINE Interp_Parvel 					     ! Using Trilinear interpolation
!===================================================================================================
IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
REAL(dbl)     :: xp,yp,zp,c00,c01,c10,c11,c0,c1,c,xd,yd,zd
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node
        IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	xp = current%pardata%xp - REAL(iMin-1_lng,dbl)
	yp = current%pardata%yp - REAL(jMin-1_lng,dbl)
	zp = current%pardata%zp - REAL(kMin-1_lng,dbl)

	ix0=FLOOR(xp)
	ix1=CEILING(xp)
	iy0=FLOOR(yp)
	iy1=CEILING(yp)
	iz0=FLOOR(zp)
	iz1=CEILING(zp)
!!!!!! MAKE SURE THE ABOVE NODES ARE FLUID NODES

	IF (ix1 /= ix0) THEN 
           xd=(xp-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
	ELSE
           xd = 0.0_dbl
	END IF

	IF (iy1 /= iy0) THEN 
           yd=(yp-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
	ELSE
           yd = 0.0_dbl
	END IF

	IF (iz1 /= iz0) THEN 
           zd=(zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
	ELSE
           zd = 0.0_dbl
	END IF

! u-interpolation
!------ 1st level linear interpolation in x-direction
	c00 = u(ix0,iy0,iz0)*(1.0_dbl-xd)+u(ix1,iy0,iz0)*xd	
	c01 = u(ix0,iy0,iz1)*(1.0_dbl-xd)+u(ix1,iy0,iz1)*xd	
	c10 = u(ix0,iy1,iz0)*(1.0_dbl-xd)+u(ix1,iy1,iz0)*xd	
	c11 = u(ix0,iy1,iz1)*(1.0_dbl-xd)+u(ix1,iy1,iz1)*xd	
!------ 2nd level linear interpolation in y-direction
	c0  = c00*(1.0_dbl-yd)+c10*yd
	c1  = c01*(1.0_dbl-yd)+c11*yd
!------ 3rd level linear interpolation in z-direction
	c   = c0*(1.0_dbl-zd)+c1*zd
        current%pardata%up=c

! v-interpolation
!------ 1st level linear interpolation in x-direction
	c00 = v(ix0,iy0,iz0)*(1.0_dbl-xd)+v(ix1,iy0,iz0)*xd
	c01 = v(ix0,iy0,iz1)*(1.0_dbl-xd)+v(ix1,iy0,iz1)*xd
	c10 = v(ix0,iy1,iz0)*(1.0_dbl-xd)+v(ix1,iy1,iz0)*xd
	c11 = v(ix0,iy1,iz1)*(1.0_dbl-xd)+v(ix1,iy1,iz1)*xd	
!------ 2nd level linear interpolation in y-direction
	c0  = c00*(1.0_dbl-yd)+c10*yd
	c1  = c01*(1.0_dbl-yd)+c11*yd
!------ 3rd level linear interpolation in z-direction
	c   = c0*(1.0_dbl-zd)+c1*zd
        current%pardata%vp=c

! w-interpolation
!------ 1st level linear interpolation in x-direction
	c00 = w(ix0,iy0,iz0)*(1.0_dbl-xd)+w(ix1,iy0,iz0)*xd	
	c01 = w(ix0,iy0,iz1)*(1.0_dbl-xd)+w(ix1,iy0,iz1)*xd	
	c10 = w(ix0,iy1,iz0)*(1.0_dbl-xd)+w(ix1,iy1,iz0)*xd	
	c11 = w(ix0,iy1,iz1)*(1.0_dbl-xd)+w(ix1,iy1,iz1)*xd	
!------ 2nd level linear interpolation in y-direction
	c0  = c00*(1.0_dbl-yd)+c10*yd
	c1  = c01*(1.0_dbl-yd)+c11*yd
!------ 3rd level linear interpolation in z-direction
	c   = c0*(1.0_dbl-zd)+c1*zd
        current%pardata%wp=c
      END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	current => next
ENDDO
!===================================================================================================
END SUBROUTINE Interp_Parvel 
!===================================================================================================




!===================================================================================================
SUBROUTINE Particle_Track
!===================================================================================================
IMPLICIT NONE

REAL(dbl)      		 :: xpold(1:np),ypold(1:np),zpold(1:np),z_tmp 		! old particle coordinates (working coordinates are stored in xp,yp,zp)
REAL(dbl)      		 :: upold(1:np),vpold(1:np),wpold(1:np) 		! old particle velocity components (new vales are stored in up, vp, wp)
REAL(dbl)                :: Cb_Domain, Cb_Local
REAL(dbl)                :: ZZZP
INTEGER(lng)   		 :: i,ipartition,ii,jj,kk
INTEGER(dbl)             :: RANK
INTEGER(lng) 		 :: mpierr
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

ParticleTransfer  = .FALSE. 							! AT this time we do not know if any particles need to be transferred.
delphi_particle   = 0.0_dbl 							! set delphi_particle to 0.0 before every time step, when the particle drug release happens. 
tausgs_particle_x = 0.0_dbl
tausgs_particle_y = 0.0_dbl
tausgs_particle_z = 0.0_dbl
	
IF (iter.GT.iter0+0_lng) THEN 	 						!At first step, the only part is finding the partitions the particles belong to without releasing any drug.  
!----Second order interpolation in time
!----Backup particle data from previous time step using a linked list of particle records

   current => ParListHead%next     
   DO WHILE (ASSOCIATED(current)) 
      next => current%next 	
      IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      current%pardata%xpold = current%pardata%xp
      current%pardata%ypold = current%pardata%yp
      current%pardata%zpold = current%pardata%zp
	
      current%pardata%upold = current%pardata%up
      current%pardata%vpold = current%pardata%vp
      current%pardata%wpold = current%pardata%wp
	
      current%pardata%xp=current%pardata%xpold+current%pardata%up
      current%pardata%yp=current%pardata%ypold+current%pardata%vp
      current%pardata%zp=current%pardata%zpold+current%pardata%wp
      END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      current => next
   ENDDO

   CALL Interp_Parvel

   current => ParListHead%next
   DO WHILE (ASSOCIATED(current))
      next => current%next 						! copy pointer of next node
      IF (mySub .EQ.current%pardata%cur_part) THEN !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      current%pardata%xp=current%pardata%xpold+0.5*(current%pardata%up+current%pardata%upold)
      current%pardata%yp=current%pardata%ypold+0.5*(current%pardata%vp+current%pardata%vpold)
      current%pardata%zp=current%pardata%zpold+0.5*(current%pardata%wp+current%pardata%wpold)
      END IF   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      current => next
   ENDDO

   CALL Interp_Parvel 							! interpolate final particle velocities after the final position is ascertained. 


!---- Parallel communication between all processors
   current => ParListHead%next
   DO WHILE (ASSOCIATED(current))
      next => current%next
           RANK= current%pardata%cur_part - 1
           CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
           CALL MPI_BCast(current%pardata%xp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
           CALL MPI_BCast(current%pardata%yp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
           CALL MPI_BCast(current%pardata%zp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
           CALL MPI_BCast(current%pardata%rp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
           CALL MPI_BCast(current%pardata%sh,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
      current => next
   ENDDO


   
!-- Particle tracking is done, now time for drug relaes calculations------------------------------------------------------------------------------------------------------------
!  CALL Interp_bulkconc(Cb_Local)  					! interpolate final bulk_concentration after the final position is ascertained.
!  CALL Calc_Global_Bulk_Scalar_Conc(Cb_Domain)
   CALL Compute_Cb  
   CALL Compute_Shear
   CALL Update_Sherwood							! Update the Sherwood number for each particle depending on the shear rate at the particle location. 
   CALL Particle_Drug_Release	  						! Updates particle radius, calculates new drug conc release rate delNBbyCV. 
   CALL Particle_Drug_To_Nodes   					! distributes released drug concentration to neightbouring nodes 
ENDIF

!---- Now update tausgs only for those cells that have non-zero values of tausgs
!DO kk=0,nzSub+1
!   DO jj=0,nySub+1
!      DO ii=0,nxSub+1
!         if (tausgs_particle_x(ii,jj,kk).ne.0.0_dbl) then
!            tausgs_particle_x(ii,jj,kk) = u(ii,jj,kk)*phi(ii,jj,kk)
!	 endif
!	 if (tausgs_particle_y(ii,jj,kk).ne.0.0_dbl) then
!            tausgs_particle_y(ii,jj,kk) = v(ii,jj,kk)*phi(ii,jj,kk)
!	 endif
!	 if (tausgs_particle_z(ii,jj,kk).ne.0.0_dbl) then
!            tausgs_particle_z(ii,jj,kk) = w(ii,jj,kk)*phi(ii,jj,kk)
!	 endif
!      ENDDO
!   ENDDO
!ENDDO


current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next ! copy pointer of next node
   IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !------- Wrappign around in z-direction for periodic BC in z
        IF (current%pardata%zp.GE.REAL(nz,dbl)) THEN
           current%pardata%zp = MOD(current%pardata%zp,REAL(nz,dbl))
        ENDIF
        IF (current%pardata%zp.LE.0.0_dbl) THEN
           current%pardata%zp = current%pardata%zp+REAL(nz,dbl)
        ENDIF

        !------- Wrappign around in y-direction for periodic BC in y
        IF (current%pardata%yp.GE.REAL(ny,dbl)) THEN
           current%pardata%yp = MOD(current%pardata%yp,REAL(ny,dbl))
        ENDIF
        IF (current%pardata%yp.LT.1.0_dbl) THEN
           current%pardata%yp = current%pardata%yp+REAL(ny,dbl)
        ENDIF

        !------- Estimate to which partition the updated position belongs to.
        DO ipartition = 1_lng,NumSubsTotal
           IF ((current%pardata%xp.GE.REAL(iMinDomain(ipartition),dbl)-1.0_dbl).AND.&
              (current%pardata%xp.LT.(REAL(iMaxDomain(ipartition),dbl)+0.0_dbl)).AND. &
              (current%pardata%yp.GE.REAL(jMinDomain(ipartition),dbl)-1.0_dbl).AND. &
              (current%pardata%yp.LT.(REAL(jMaxDomain(ipartition),dbl)+0.0_dbl)).AND. &
              (current%pardata%zp.GE.REAL(kMinDomain(ipartition),dbl)-1.0_dbl).AND. &
              (current%pardata%zp.LT.(REAL(kMaxDomain(ipartition),dbl)+0.0_dbl))) THEN
                 current%pardata%new_part = ipartition
            END IF
        END DO

        IF ((.NOT.ParticleTransfer).AND.(current%pardata%new_part .NE. current%pardata%cur_part)) THEN
           ParticleTransfer = .TRUE.
        END IF
   END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   current => next
ENDDO

!---- Parallel communication between all processors
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
	RANK= current%pardata%cur_part - 1
        current%pardata%cur_part = current%pardata%new_part 

	CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	CALL MPI_BCast(current%pardata%xp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%yp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%zp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%up,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%vp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%wp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%rp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%sh,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%xpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%ypold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%zpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%upold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%vpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%wpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%rpold,     1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%delNBbyCV, 1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%par_conc,  1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%bulk_conc, 1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%gamma_cont,1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%Nbj,       1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%S,         1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%Sst,       1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
	CALL MPI_BCast(current%pardata%cur_part,  1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%new_part,  1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
        
   current => next  
ENDDO
!===================================================================================================
END SUBROUTINE Particle_Track
!===================================================================================================





!================================================
END MODULE ParticleTracking 
!================================================
