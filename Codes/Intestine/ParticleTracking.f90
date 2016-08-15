!==================================================================================================
MODULE ParticleTracking				
!==================================================================================================
USE SetPrecision
USE Setup
USE IC
USE MPI
USE ParticleDrug

IMPLICIT NONE

CONTAINS

!===================================================================================================
SUBROUTINE Particle_Setup
!===================================================================================================
IMPLICIT NONE
IF (Flag_Restart) THEN
ELSE
   CALL Particle_Velocity
ENDIF
!===================================================================================================
END SUBROUTINE Particle_Setup
!===================================================================================================





!===================================================================================================
SUBROUTINE Particle_Velocity 					     ! Using Trilinear interpolation
!===================================================================================================
IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
INTEGER(lng)  :: ii,jj,kk
REAL(dbl)     :: xp,yp,zp,c00,c01,c10,c11,c0,c1,c,xd,yd,zd
REAL(dbl)     :: xaxis,yaxis,X_s,Y_s,R_s,CosTheta_s,SinTheta_s,Vel_s
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next

   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                                           !only calculate the drug release when particle radius is larger than 0.1 micron
      IF (mySub .EQ.current%pardata%cur_part) THEN
         xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
         yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
         zp= current%pardata%zp - REAL(kMin-1_lng,dbl)

         ix0= FLOOR(xp)
         ix1= CEILING(xp)
         iy0= FLOOR(yp)
         iy1= CEILING(yp)
         iz0= FLOOR(zp)
         iz1= CEILING(zp)

!------- Treating the solid nodes so that their velocity is not zero for interpolation purposes
!------- Velocity magnitude is calculated based on boundary velocity at that z location
!------- Velocity vector points to the center od the circle in that Z-location
         DO ii= ix0, ix1
            DO jj= iy0, iy1
               DO kk= iz0, iz1
                  IF (node(ii,jj,kk) .EQ. SOLID) THEN					! Solid nodes in the lattice cell encompassing the particle
                     xaxis 	   = ANINT(0.5_dbl*(nx+1))
                     yaxis 	   = ANINT(0.5_dbl*(ny+1))
                     X_s   	   = xcf* (ii+ (iMin-1_lng)- xaxis)
                     Y_s   	   = ycf* (jj+ (jMin-1_lng)- yaxis) 
                     R_s   	   = SQRT(X_s**2 + Y_s**2)
                     CosTheta_s    = X_s/R_s
                     SinTheta_s    = Y_s/R_s
                     Vel_s	   = vel(kk)
                     u_s(ii,jj,kk) = Vel_s * CosTheta_s
                     v_s(ii,jj,kk) = Vel_s * SinTheta_s
                     w_s(ii,jj,kk) = 0.0_dbl				
                  ELSE 									! Fluid nodes in the lattice cell encompassing the particle
                     u_s(ii,jj,kk) = u(ii,jj,kk) 
                     v_s(ii,jj,kk) = v(ii,jj,kk) 
                     w_s(ii,jj,kk) = w(ii,jj,kk)
                  END IF
               END DO
            END DO
         END DO

!------- Preparing for tri-linear interpolation
         IF (ix1 .NE. ix0) THEN 
            xd= (xp-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
         ELSE
            xd= 0.0_dbl
         END IF

         IF (iy1 .NE. iy0) THEN 
            yd= (yp-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
         ELSE
            yd= 0.0_dbl
         END IF

         IF (iz1 .NE. iz0) THEN 
            zd= (zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
         ELSE
            zd= 0.0_dbl
         END IF

!--------u-interpolation
!--------1st level linear interpolation in x-direction
         c00= u_s(ix0,iy0,iz0)*(1.0_dbl-xd)+u_s(ix1,iy0,iz0)*xd	
         c01= u_s(ix0,iy0,iz1)*(1.0_dbl-xd)+u_s(ix1,iy0,iz1)*xd	
         c10= u_s(ix0,iy1,iz0)*(1.0_dbl-xd)+u_s(ix1,iy1,iz0)*xd	
         c11= u_s(ix0,iy1,iz1)*(1.0_dbl-xd)+u_s(ix1,iy1,iz1)*xd	
!--------2nd level linear interpolation in y-direction
         c0 = c00*(1.0_dbl-yd)+c10*yd
         c1 = c01*(1.0_dbl-yd)+c11*yd
!--------3rd level linear interpolation in z-direction
         c  = c0*(1.0_dbl-zd)+c1*zd
         current%pardata%up=c

!--------v-interpolation
!--------1st level linear interpolation in x-direction
         c00= v_s(ix0,iy0,iz0)*(1.0_dbl-xd)+v_s(ix1,iy0,iz0)*xd
         c01= v_s(ix0,iy0,iz1)*(1.0_dbl-xd)+v_s(ix1,iy0,iz1)*xd
         c10= v_s(ix0,iy1,iz0)*(1.0_dbl-xd)+v_s(ix1,iy1,iz0)*xd
         c11= v_s(ix0,iy1,iz1)*(1.0_dbl-xd)+v_s(ix1,iy1,iz1)*xd	
!--------2nd level linear interpolation in y-direction
         c0 = c00*(1.0_dbl-yd)+c10*yd
         c1 = c01*(1.0_dbl-yd)+c11*yd
!--------3rd level linear interpolation in z-direction
         c  = c0*(1.0_dbl-zd)+c1*zd
         current%pardata%vp=c

!--------w-interpolation
!--------1st level linear interpolation in x-direction
         c00 = w_s(ix0,iy0,iz0)*(1.0_dbl-xd)+w_s(ix1,iy0,iz0)*xd	
         c01 = w_s(ix0,iy0,iz1)*(1.0_dbl-xd)+w_s(ix1,iy0,iz1)*xd	
         c10 = w_s(ix0,iy1,iz0)*(1.0_dbl-xd)+w_s(ix1,iy1,iz0)*xd	
         c11 = w_s(ix0,iy1,iz1)*(1.0_dbl-xd)+w_s(ix1,iy1,iz1)*xd	
!--------2nd level linear interpolation in y-direction
         c0  = c00*(1.0_dbl-yd)+c10*yd
         c1  = c01*(1.0_dbl-yd)+c11*yd
!--------3rd level linear interpolation in z-direction
         c   = c0*(1.0_dbl-zd)+c1*zd
         current%pardata%wp=c
      END IF 
   END IF
   current => next
ENDDO
!===================================================================================================
END SUBROUTINE Particle_Velocity 
!===================================================================================================





!===================================================================================================
SUBROUTINE Particle_Track   					!Second order interpolation in time
!===================================================================================================
IMPLICIT NONE

INTEGER(lng) :: mpierr
INTEGER(lng) :: ix0,ix1,iy0,iy1,iz0,iz1, Number_of_Solid_nodes
REAL(dbl)    :: xaxis,yaxis,CosTheta_p,SinTheta_p
REAL(dbl)    :: xpp,ypp
REAL(dbl)    :: Min_R_Acceptable
REAL(dbl)    :: xp,yp,zp,zd
REAL(dbl)    :: R_Particle, R_Boundary
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

!----- set delphi_particle to 0.0 before every time step, when the particle drug release happens. 
delphi_particle   = 0.0_dbl 
tausgs_particle_x = 0.0_dbl
tausgs_particle_y = 0.0_dbl
tausgs_particle_z = 0.0_dbl           
Min_R_Acceptable  = 1.0e-7						! 0.1 micron is the minimum acceptable particle size


!IF (iter.GT.iter0+0_lng) THEN 
   current => ParListHead%next     
   DO WHILE (ASSOCIATED(current)) 
      next => current%next 	
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN
         IF (mySub .EQ.current%pardata%cur_part) THEN 
            current%pardata%xpold = current%pardata%xp
            current%pardata%ypold = current%pardata%yp
            current%pardata%zpold = current%pardata%zp
            current%pardata%upold = current%pardata%up
            current%pardata%vpold = current%pardata%vp
            current%pardata%wpold = current%pardata%wp
            current%pardata%xp=current%pardata%xpold+current%pardata%up
            current%pardata%yp=current%pardata%ypold+current%pardata%vp
            current%pardata%zp=current%pardata%zpold+current%pardata%wp
            IF ((current%pardata%zp .LT. nz) .AND. (current%pardata%zp .GT. 1))THEN
               !----- If Particle leaves the fluid domain---------------------------------------
               xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
               yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
               zp= current%pardata%zp - REAL(kMin-1_lng,dbl)

               ix0= FLOOR(xp)
               ix1= CEILING(xp)
               iy0= FLOOR(yp)
               iy1= CEILING(yp)
               iz0= FLOOR(zp)
               iz1= CEILING(zp)
            
               IF (iz1 .NE. iz0) THEN 
                  zd= (zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
               ELSE 
                  zd= 0.0_dbl
               END IF
            
               R_boundary = r(iz0)*(1.0_dbl-zd) + r(iz1)*zd
               xaxis= 0.5_dbl*(nx+1)
               yaxis= 0.5_dbl*(ny+1)
               xpp = xcf*(current%pardata%xp- xaxis) 
               ypp = ycf*(current%pardata%yp- yaxis) 
               R_Particle    = SQRT(xpp**2 + ypp**2)
               CosTheta_p    = xpp/R_Particle
               SinTheta_p    = ypp/R_Particle
               IF (R_Particle .GT. R_Boundary) THEN  						! Check if particle is outside analytical boundary
                  write(*,*) '=========================================================================================='
                  write(*,*) 'Iter, parID,xp,yp,zp:', iter,current%pardata%parid, current%pardata%xp, current%pardata%yp, current%pardata%zp
                  write(*,*) 'R_Particle, R_Boundary:', R_Particle,R_Boundary
                  Number_of_Solid_nodes =   node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1) 
                  write(*,*) 'No of solid nodes', Number_of_Solid_nodes 
                  write(*,*) '---------------------------------------------------------------------------'
                  write(*,*) 'Treating the particle:'
                  xpp = R_Boundary * CosTheta_p
                  ypp = R_Boundary * SinTheta_p
                  R_Particle = SQRT(xpp**2 + ypp**2)
                  current%pardata%xp = (xpp/xcf) + xaxis
                  current%pardata%yp = (ypp/ycf) + yaxis 
                  write(*,*) 'Iter, parID,xp,yp,zp:', iter,current%pardata%parid, current%pardata%xp, current%pardata%yp, current%pardata%zp
                  write(*,*) 'R_Particle, R_Boundary:', R_Particle,R_Boundary
                  Number_of_Solid_nodes =   node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1) 
                  write(*,*) 'No of solid nodes', Number_of_Solid_nodes 
               END IF
               IF ((node(ix0,iy0,iz0) .EQ. SOLID) .AND. (node(ix1,iy0,iz0) .EQ. SOLID) .AND. &		! Check if all nodes around particle are solid
                  (node(ix0,iy1,iz0) .EQ. SOLID) .AND. (node(ix0,iy0,iz1) .EQ. SOLID) .AND. &
                  (node(ix1,iy1,iz0) .EQ. SOLID) .AND. (node(ix1,iy0,iz1) .EQ. SOLID) .AND. &
                  (node(ix0,iy1,iz1) .EQ. SOLID) .AND. (node(ix1,iy1,iz1) .EQ. SOLID)) THEN
                  write(*,*) '=========================================================================================='
                  write(*,*) '1st part of tracking: All nodes around are solid'    
                  write(*,*) 'Iter,ParID,x,y,z:', iter,current%pardata%parid,current%pardata%xp, current%pardata%yp, current%pardata%zp
               END IF !--------------------------------------------------------------------------------------   
            END IF
         END IF
      END IF
      current => next
   END DO

   CALL Particle_Transfer
   CALL Particle_Velocity

   current => ParListHead%next
   DO WHILE (ASSOCIATED(current))
      next => current%next 						
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN
         IF (mySub .EQ.current%pardata%cur_part) THEN 
            current%pardata%xp=current%pardata%xpold+0.5*(current%pardata%up+current%pardata%upold)
            current%pardata%yp=current%pardata%ypold+0.5*(current%pardata%vp+current%pardata%vpold)
            current%pardata%zp=current%pardata%zpold+0.5*(current%pardata%wp+current%pardata%wpold)
            
            IF ((current%pardata%zp .LT. nz) .AND. (current%pardata%zp .GT. 1) ) THEN 
               !----- If Particle leaves the fluid domain---------------------------------------
               xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
               yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
               zp= current%pardata%zp - REAL(kMin-1_lng,dbl)

               ix0= FLOOR(xp)
               ix1= CEILING(xp)
               iy0= FLOOR(yp)
               iy1= CEILING(yp)
               iz0= FLOOR(zp)
               iz1= CEILING(zp)

               IF (iz1 .NE. iz0) THEN
                  zd= (zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
               ELSE 
                  zd= 0.0_dbl
               END IF
               R_boundary = r(iz0)*(1.0_dbl-zd) + r(iz1)*zd
            
               xaxis= 0.5_dbl*(nx+1)
               yaxis= 0.5_dbl*(ny+1)
               xpp = xcf*(current%pardata%xp- xaxis)
               ypp = ycf*(current%pardata%yp- yaxis)
               R_Particle    = SQRT(xpp**2 + ypp**2)
               CosTheta_p    = xpp/R_Particle
               SinTheta_p    = ypp/R_Particle
               IF (R_Particle .GT. R_Boundary) THEN                                                ! Check if particle is outside analytical boundary
                  write(*,*) '=========================================================================================='
                  write(*,*) 'Iter, parID,xp,yp,zp:', iter,current%pardata%parid, current%pardata%xp, current%pardata%yp, current%pardata%zp
                  write(*,*) 'R_Particle, R_Boundary:', R_Particle,R_Boundary
                  write(*,*) 'No of solid nodes', node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1)  
                  write(*,*) '---------------------------------------------------------------------------'
                  write(*,*) 'Treating the particle:'
                  xpp = R_Boundary * CosTheta_p
                  ypp = R_Boundary * SinTheta_p
                  R_Particle = SQRT(xpp**2 + ypp**2)
                  current%pardata%xp = (xpp/xcf) + xaxis
                  current%pardata%yp = (ypp/ycf) + yaxis
                  write(*,*) 'Iter, parID,xp,yp,zp:', iter,current%pardata%parid, current%pardata%xp, current%pardata%yp, current%pardata%zp
                  write(*,*) 'R_Particle, R_Boundary:', R_Particle,R_Boundary
                  write(*,*) 'No of solid nodes', node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1)  
               END IF
               IF ((node(ix0,iy0,iz0) .EQ. SOLID) .AND. (node(ix1,iy0,iz0) .EQ. SOLID) .AND. &             ! Check if all nodes around particle are solid
                  (node(ix0,iy1,iz0) .EQ. SOLID) .AND. (node(ix0,iy0,iz1) .EQ. SOLID) .AND. &
                  (node(ix1,iy1,iz0) .EQ. SOLID) .AND. (node(ix1,iy0,iz1) .EQ. SOLID) .AND. &
                  (node(ix0,iy1,iz1) .EQ. SOLID) .AND. (node(ix1,iy1,iz1) .EQ. SOLID)) THEN
                  write(*,*) '=========================================================================================='
                  write(*,*) '2nd part of tracking: All nodes around are solid'
                  write(*,*) 'Iter,ParID,x,y,z:', iter,current%pardata%parid,current%pardata%xp, current%pardata%yp, current%pardata%zp
               END IF !--------------------------------------------------------------------------------------
            ENDIF 
         END IF
      END IF
      current => next
   ENDDO

   CALL Particle_Transfer
   CALL Particle_Velocity 		

!--Particle tracking is done, now time for drug relaes calculations---------------------------------
   IF (Flag_Shear_Effects) THEN
      CALL Compute_Shear
   END IF   

   CALL Compute_Sherwood                  ! Update the Sherwood number for each particle depending on the shear rate. 

   IF (Flag_Golbal_Bulk_Concentration) THEN
      CALL Compute_Cb_Global
   ELSE
      CALL Compute_Cb
   END IF

   CALL Particle_Drug_Release             ! Updates particle radius, calculates drug release rate delNBbyCV. 
   CALL Particle_Drug_To_Nodes            ! distributes released drug concentration to nodes in effective volume. 
!  CALL Particle_History                  ! Keep trak of a few particles
   CALL Particle_Transfer 
!ENDIF
!===================================================================================================
END SUBROUTINE Particle_Track
!===================================================================================================





!===================================================================================================
SUBROUTINE Particle_Transfer
!===================================================================================================
IMPLICIT NONE

INTEGER(lng)   		 :: i,ipartition,ii,jj,kk
INTEGER(lng)             :: RANK
INTEGER(lng)             :: mpierr
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 	

   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN
      IF (mySub .EQ.current%pardata%cur_part) THEN 
         !----- Wrappign around in z-direction for periodic BC in z
         IF (current%pardata%zp.GE.REAL(nz,dbl)) THEN
            current%pardata%zp = MOD(current%pardata%zp,REAL(nz,dbl))
         ENDIF
         IF (current%pardata%zp.LT.0.0_dbl) THEN
            current%pardata%zp = current%pardata%zp+REAL(nz,dbl)
         ENDIF

         !----- Wrappign around in y-direction for periodic BC in y
         IF (current%pardata%yp.GE.REAL(ny,dbl)) THEN
            current%pardata%yp = MOD(current%pardata%yp,REAL(ny,dbl))
         ENDIF
         IF (current%pardata%yp.LT.1.0_dbl) THEN
            current%pardata%yp = current%pardata%yp+REAL(ny,dbl)
         ENDIF

         !----- Estimate to which partition the updated position belongs to.
         DO ipartition = 1_lng,NumSubsTotal
            IF((current%pardata%xp.GE.REAL(iMinDomain(ipartition),dbl)-1.0_dbl).AND.&
               (current%pardata%xp.LT.(REAL(iMaxDomain(ipartition),dbl)+0.0_dbl)).AND. &
               (current%pardata%yp.GE.REAL(jMinDomain(ipartition),dbl)-1.0_dbl).AND. &
               (current%pardata%yp.LT.(REAL(jMaxDomain(ipartition),dbl)+0.0_dbl)).AND. &
               (current%pardata%zp.GE.REAL(kMinDomain(ipartition),dbl)-1.0_dbl).AND. &
               (current%pardata%zp.LT.(REAL(kMaxDomain(ipartition),dbl)+0.0_dbl))) THEN
               current%pardata%new_part = ipartition
            END IF
         END DO
      END IF
   END IF

   current => next
ENDDO

!---- Parallel communication between all processors
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 

   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN
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
   END IF
   current => next  
ENDDO
!===================================================================================================
END SUBROUTINE Particle_Transfer
!===================================================================================================





!===================================================================================================
SUBROUTINE Particle_History
!===================================================================================================
IMPLICIT NONE

TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next
 
   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN
      SELECT CASE(current%pardata%parid)
         CASE(1_lng)
         IF (mySub .EQ.current%pardata%cur_part) THEN
            open(72,file='History-Particle-1.dat',position='append')
            write(72,101) iter,iter*tcf,current%pardata%xp, current%pardata%yp, current%pardata%zp, &
            current%pardata%up, current%pardata%vp, current%pardata%wp, &
            current%pardata%sh, current%pardata%rp, current%pardata%bulk_conc, current%pardata%delNBbyCV
            close(72)
         END IF

         CASE(2_lng)
         IF (mySub .EQ.current%pardata%cur_part) THEN
            open(73,file='History-Particle-2.dat',position='append')
            write(73,101) iter,iter*tcf,current%pardata%xp, current%pardata%yp, current%pardata%zp, &
            current%pardata%up, current%pardata%vp, current%pardata%wp, &
            current%pardata%sh, current%pardata%rp, current%pardata%bulk_conc, current%pardata%delNBbyCV
            close(73)
         END IF

         CASE(3_lng)
         IF (mySub .EQ.current%pardata%cur_part) THEN
            open(74,file='History-Particle-3.dat',position='append')
            write(74,101) iter,iter*tcf,current%pardata%xp, current%pardata%yp, current%pardata%zp, &
            current%pardata%up, current%pardata%vp, current%pardata%wp, &
            current%pardata%sh, current%pardata%rp, current%pardata%bulk_conc, current%pardata%delNBbyCV
            close(74)
          END IF
                 
         CASE(4_lng)
         IF (mySub .EQ.current%pardata%cur_part) THEN
            open(75,file='History-Particle-4.dat',position='append')
            write(75,101) iter,iter*tcf,current%pardata%xp, current%pardata%yp, current%pardata%zp, &
            current%pardata%up, current%pardata%vp, current%pardata%wp, &
            current%pardata%sh, current%pardata%rp, current%pardata%bulk_conc, current%pardata%delNBbyCV
            close(75)
         END IF
      END SELECT
101 format (I6, F10.3, 7F20.14,3E20.12)
   END IF
   current => next  
ENDDO

!===================================================================================================
END SUBROUTINE Particle_History
!===================================================================================================





!================================================
END MODULE ParticleTracking 
!================================================
