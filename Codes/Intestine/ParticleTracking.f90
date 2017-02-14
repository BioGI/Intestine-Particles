!==================================================================================================
MODULE ParticleTracking				
!==================================================================================================
USE SetPrecision
USE Setup
USE IC
USE MPI
USE ParticleDrug
USE Output
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

   IF (mySub .EQ.current%pardata%cur_part) THEN
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                                           !only calculate the drug release when particle radius is larger than 0.1 micron
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
         IF (Flag_Couette) THEN !---- Couette simulation -------------------------------------------
            u_s= u 
            v_s= v 
            w_s= w
         ELSE !---------------------- Intestine Simulation -----------------------------------------
            DO ii= ix0, ix1
               DO jj= iy0, iy1
                  DO kk= iz0, iz1
                     IF (node(ii,jj,kk) .EQ. SOLID) THEN			! Solid nodes in the lattice cell encompassing the particle
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

         END IF   
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
REAL(dbl)    :: xp,yp,zp,zd
REAL(dbl)    :: R_Particle, R_Boundary
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

!----- set delphi_particle to 0.0 before every time step, when the particle drug release happens. 
delphi_particle   = 0.0_dbl 
tausgs_particle_x = 0.0_dbl
tausgs_particle_y = 0.0_dbl
tausgs_particle_z = 0.0_dbl           


current => ParListHead%next     
DO WHILE (ASSOCIATED(current)) 
   next => current%next 

   IF (mySub .EQ.current%pardata%cur_part) THEN 
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN
         current%pardata%xpold = current%pardata%xp
         current%pardata%ypold = current%pardata%yp
         current%pardata%zpold = current%pardata%zp
         current%pardata%upold = current%pardata%up
         current%pardata%vpold = current%pardata%vp
         current%pardata%wpold = current%pardata%wp
         current%pardata%xp=current%pardata%xpold+current%pardata%up
         current%pardata%yp=current%pardata%ypold+current%pardata%vp
         current%pardata%zp=current%pardata%zpold+current%pardata%wp

         IF (Flag_Couette) THEN !----- Couette simulation, no need to make sure particles do not leave the fluid domain
         ELSE !----------------------- Intestine simulation: make sure particles do not leave the fluid domain    
          IF ((current%pardata%zp .LT. nz) .AND. (current%pardata%zp .GT. 1))THEN              
             xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
             yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
             zp= current%pardata%zp - REAL(kMin-1_lng,dbl)
             ix0= FLOOR(xp)
             ix1= CEILING(xp)
             iy0= FLOOR(yp)
             iy1= CEILING(yp)
             iz0= FLOOR(zp)
             iz1= CEILING(zp)
             Number_of_Solid_nodes= node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1)
             IF (Number_of_Solid_nodes .GT. 0) THEN         !there is a solid node around the particle
                IF (iz1 .NE. iz0) THEN 
                  zd= (zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
                ELSE 
                  zd= 0.0_dbl
                END IF
                R_Boundary = r(iz0)*(1.0_dbl-zd) + r(iz1)*zd      !radius of solid boundary at z location of the particle 
                xaxis= 0.5_dbl*(nx+1)
                yaxis= 0.5_dbl*(ny+1)
                xpp = xcf*(current%pardata%xp- xaxis) 
                ypp = ycf*(current%pardata%yp- yaxis) 
                R_Particle    = SQRT(xpp**2 + ypp**2)              !radius at location of the particle
                CosTheta_p    = xpp/R_Particle
                SinTheta_p    = ypp/R_Particle
            
                IF ((R_Boundary.GT.xcf).AND.(R_Particle .GT. (R_Boundary-xcf))) THEN ! particle is outside analytical boundary
                   write(*,*) '=========================================================================================='
                   write(*,*) 'A:Iter, parID,xp,yp,zp,zd,Rz1,Rz2:', iter,current%pardata%parid, current%pardata%xp, current%pardata%yp, current%pardata%zp,zd,r(iz0),r(iz1)
                   write(*,*) 'R_Particle, R_Boundary:', R_Particle,R_Boundary
                   Number_of_Solid_nodes =   node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1) 
                   write(*,*) 'No of solid nodes', Number_of_Solid_nodes 
                   write(*,*) '---------------------------------------------------------------------------'
                   write(*,*) 'Treating the particle:'
                   xpp = (R_Boundary-1.20_dbl*xcf) * CosTheta_p
                   ypp = (R_Boundary-1.20_dbl*xcf) * SinTheta_p
                   R_Particle = SQRT(xpp**2 + ypp**2)
                   current%pardata%xp = (xpp/xcf) + xaxis
                   current%pardata%yp = (ypp/ycf) + yaxis 
                   write(*,*) 'A:Iter, parID,xp,yp,zp,zd,Rz1,Rz2:', iter,current%pardata%parid, current%pardata%xp, current%pardata%yp, current%pardata%zp, zd,r(iz0),r(iz1)
                   write(*,*) 'R_Particle, R_Boundary:', R_Particle,R_Boundary
                   Number_of_Solid_nodes =   node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1) 
                   write(*,*) 'No of solid nodes', Number_of_Solid_nodes 
                END IF
             END IF               ! If Number_of_Solid_nodes > 0 
          END IF                  ! If 1<zp<nz        
        END IF                    ! If this simulation is Couette 
      END IF                      ! If particle radius is over Min_R_Acceptable
   END IF                         ! If particle resides in this processor
   current => next
END DO

CALL Particle_Transfer
CALL Particle_Velocity

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next
   IF (mySub .EQ.current%pardata%cur_part) THEN 
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN
         current%pardata%xp=current%pardata%xpold+0.5*(current%pardata%up+current%pardata%upold)
         current%pardata%yp=current%pardata%ypold+0.5*(current%pardata%vp+current%pardata%vpold)
         current%pardata%zp=current%pardata%zpold+0.5*(current%pardata%wp+current%pardata%wpold)
         IF (Flag_Couette) THEN !----- Couette simulation, no need to make sure particles do not leave the fluid domain
         ELSE !----------------------- Intestine simulation: make sure particles do not leave the fluid domain    
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
               Number_of_Solid_nodes= node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1)
               IF (Number_of_Solid_nodes .GT. 0) THEN !a solid node is around the particle
                  IF (iz1 .NE. iz0) THEN
                     zd= (zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
                  ELSE 
                     zd= 0.0_dbl
                  END IF
                  R_Boundary = r(iz0)*(1.0_dbl-zd) + r(iz1)*zd
                  xaxis= 0.5_dbl*(nx+1)
                  yaxis= 0.5_dbl*(ny+1)
                  xpp = xcf*(current%pardata%xp- xaxis)
                  ypp = ycf*(current%pardata%yp- yaxis)
                  R_Particle    = SQRT(xpp**2 + ypp**2)
                  CosTheta_p    = xpp/R_Particle
                  SinTheta_p    = ypp/R_Particle
                  IF ((R_Boundary .GT. xcf).AND.(R_Particle .GT. (R_Boundary-xcf))) THEN  !particle is outside analytical boundary
                    write(*,*) '=========================================================================================='
                    write(*,*) 'B:Iter, parID,xp,yp,zp,zd,iz0,iz1,Rz1,Rz2:', iter,current%pardata%parid, current%pardata%xp, current%pardata%yp, current%pardata%zp,zd,iz0,iz1,r(iz0),r(iz1)
                    write(*,*) 'B: CPU,iMin,iMax,jMin,jMax,kMin,kMax',myid,iMin,iMax,jMin,jMax,kMin,kMax
                    write(*,*) 'R_Particle, R_Boundary:', R_Particle,R_Boundary
                    write(*,*) 'No of solid nodes', node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1)  
                    write(*,*) '---------------------------------------------------------------------------'
                    write(*,*) 'Treating the particle:'
                    xpp = (R_Boundary-1.20_dbl*xcf) * CosTheta_p
                    ypp = (R_Boundary-1.20_dbl*xcf) * SinTheta_p
                    R_Particle = SQRT(xpp**2 + ypp**2)
                    current%pardata%xp = (xpp/xcf) + xaxis
                    current%pardata%yp = (ypp/ycf) + yaxis
                    write(*,*) 'B:Iter, parID,xp,yp,zp,zd,iz0,iz1,Rz1,Rz2:', iter,current%pardata%parid, current%pardata%xp, current%pardata%yp, current%pardata%zp,zd,iz0,iz1,r(iz0),r(iz1)
                    write(*,*) 'B: CPU,iMin,iMax,jMin,jMax,kMin,kMax',myid,iMin,iMax,jMin,jMax,kMin,kMax
                    write(*,*) 'R_Particle, R_Boundary:', R_Particle,R_Boundary
                    write(*,*) 'No of solid nodes', node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1)  
                  END IF
               END IF   ! If Number_of_Solid_nodes > 0
            ENDIF       ! If 1<zp<nz
         END IF         ! If the simulation is Couette
      END IF            ! If particle's radius is larger than Min_R_Acceptable
   END IF               ! If particle resides inthis processor
   current => next
ENDDO

CALL PrintComputationalTime(2)
CALL Particle_Transfer
CALL PrintComputationalTime(3)
CALL Particle_Velocity 		
CALL PrintComputationalTime(4)
!-----Particle tracking is done, now time for drug relaes calculations---------------------------------
CALL Compute_C_bulk  
CALL PrintComputationalTime(5)
IF (Flag_Shear_Effects) THEN
   CALL Compute_Shear
END IF   
CALL Compute_Sherwood             ! Update the Sherwood number for each particle depending on the shear rate. 
CALL Compute_C_surface
CALL PrintComputationalTime(6)
CALL Particle_Drug_Release        ! Updates particle radius, calculates drug release rate delNBbyCV. 
CALL PrintComputationalTime(7)
CALL Particle_Transfer 
CALL PrintComputationalTime(8)
CALL Particle_Drug_To_Nodes       ! distributes released drug concentration to nodes in effective volume. 
CALL PrintComputationalTime(9)
!CALL Particle_Transfer 
!CALL Particle_History             ! Keep trak of a few particles
!===================================================================================================
END SUBROUTINE Particle_Track
!===================================================================================================





!===================================================================================================
SUBROUTINE Particle_Transfer
!===================================================================================================
IMPLICIT NONE

REAL(dbl)                :: Particle_Data_l(np,23), PArticle_Data_g(np,23)
INTEGER(lng)   		       :: i,ipartition,ii,jj,kk
INTEGER(lng)             :: PID, ID, RANK
INTEGER(lng)             :: mpierr
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 	
   ID = current%pardata%parid

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
  
      Particle_Data_l(ID,1) = current%pardata%xp
      Particle_Data_l(ID,2) = current%pardata%yp
      Particle_Data_l(ID,3) = current%pardata%zp
      Particle_Data_l(ID,4) = current%pardata%up
      Particle_Data_l(ID,5) = current%pardata%vp
      Particle_Data_l(ID,6) = current%pardata%wp
      Particle_Data_l(ID,7) = current%pardata%rp
      Particle_Data_l(ID,8) = current%pardata%sh_conf
      Particle_Data_l(ID,9) = current%pardata%sh_shear
      Particle_Data_l(ID,10)= current%pardata%sh_slip
      Particle_Data_l(ID,11)= current%pardata%xpold
      Particle_Data_l(ID,12)= current%pardata%ypold
      Particle_Data_l(ID,13)= current%pardata%zpold
      Particle_Data_l(ID,14)= current%pardata%upold
      Particle_Data_l(ID,15)= current%pardata%vpold
      Particle_Data_l(ID,16)= current%pardata%wpold
      Particle_Data_l(ID,17)= current%pardata%rpold
      Particle_Data_l(ID,18)= current%pardata%delNBbyCV
      Particle_Data_l(ID,19)= current%pardata%par_conc
      Particle_Data_l(ID,20)= current%pardata%bulk_conc
      Particle_Data_l(ID,21)= current%pardata%S
      Particle_Data_l(ID,22)= current%pardata%Sst
      Particle_Data_l(ID,23)= current%pardata%new_part
   ELSE
      Particle_Data_l(ID,1) = 1e5                          
      Particle_Data_l(ID,2) = 1e5                    
      Particle_Data_l(ID,3) = 1e5                     
      Particle_Data_l(ID,4) = 1e5                     
      Particle_Data_l(ID,5) = 1e5                     
      Particle_Data_l(ID,6) = 1e5                     
      Particle_Data_l(ID,7) = 1e5                     
      Particle_Data_l(ID,8) = 1e5                     
      Particle_Data_l(ID,9) = 1e5                     
      Particle_Data_l(ID,10)= 1e5                    
      Particle_Data_l(ID,11)= 1e5                     
      Particle_Data_l(ID,12)= 1e5                     
      Particle_Data_l(ID,13)= 1e5                     
      Particle_Data_l(ID,14)= 1e5                     
      Particle_Data_l(ID,15)= 1e5                     
      Particle_Data_l(ID,16)= 1e5                     
      Particle_Data_l(ID,17)= 1e5                     
      Particle_Data_l(ID,18)= 1e5                     
      Particle_Data_l(ID,19)= 1e5                     
      Particle_Data_l(ID,20)= 1e5                     
      Particle_Data_l(ID,21)= 1e5                     
      Particle_Data_l(ID,22)= 1e5                     
      Particle_Data_l(ID,23)= 1e5                     
   ENDIF
   current => next
ENDDO

!---- Parallel communication between all processors
CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
CALL MPI_ALLREDUCE(Particle_Data_l, Particle_Data_g, np*22, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)


current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
      ID= current%pardata%parid
      current%pardata%xp=         Particle_Data_g(ID,1)
      current%pardata%yp=         Particle_Data_g(ID,2)
      current%pardata%zp=         Particle_Data_g(ID,3)
      current%pardata%up=         Particle_Data_g(ID,4)
      current%pardata%vp=         Particle_Data_g(ID,5)
      current%pardata%wp=         Particle_Data_g(ID,6)
      current%pardata%rp=         Particle_Data_g(ID,7)
      current%pardata%sh_conf=    Particle_Data_g(ID,8)
      current%pardata%sh_shear=   Particle_Data_g(ID,9)
      current%pardata%sh_slip=    Particle_Data_g(ID,10)
      current%pardata%xpold=      Particle_Data_g(ID,11)
      current%pardata%ypold=      Particle_Data_g(ID,12)
      current%pardata%zpold=      Particle_Data_g(ID,13)
      current%pardata%upold=      Particle_Data_g(ID,14)
      current%pardata%vpold=      Particle_Data_g(ID,15)
      current%pardata%wpold=      Particle_Data_g(ID,16)
      current%pardata%rpold=      Particle_Data_g(ID,17)
      current%pardata%delNBbyCV=  Particle_Data_g(ID,18)
      current%pardata%par_conc=   Particle_Data_g(ID,19)
      current%pardata%bulk_conc=  Particle_Data_g(ID,20)
      current%pardata%S=          Particle_Data_g(ID,21)
      current%pardata%Sst=        Particle_Data_g(ID,22)
	    current%pardata%cur_part=   Particle_Data_g(ID,23)
      current%pardata%new_part=   Particle_Data_g(ID,23)
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
            current%pardata%sh_shear, current%pardata%rp, current%pardata%bulk_conc, current%pardata%delNBbyCV
            close(72)
         END IF

         CASE(2_lng)
         IF (mySub .EQ.current%pardata%cur_part) THEN
            open(73,file='History-Particle-2.dat',position='append')
            write(73,101) iter,iter*tcf,current%pardata%xp, current%pardata%yp, current%pardata%zp, &
            current%pardata%up, current%pardata%vp, current%pardata%wp, &
            current%pardata%sh_shear, current%pardata%rp, current%pardata%bulk_conc, current%pardata%delNBbyCV
            close(73)
         END IF

         IF (mySub .EQ.current%pardata%cur_part) THEN
            open(74,file='History-Particle-3.dat',position='append')
            write(74,101) iter,iter*tcf,current%pardata%xp, current%pardata%yp, current%pardata%zp, &
            current%pardata%up, current%pardata%vp, current%pardata%wp, &
            current%pardata%sh_shear, current%pardata%rp, current%pardata%bulk_conc, current%pardata%delNBbyCV
            close(74)
          END IF
                 
         CASE(4_lng)
         IF (mySub .EQ.current%pardata%cur_part) THEN
            open(75,file='History-Particle-4.dat',position='append')
            write(75,101) iter,iter*tcf,current%pardata%xp, current%pardata%yp, current%pardata%zp, &
            current%pardata%up, current%pardata%vp, current%pardata%wp, &
            current%pardata%sh_shear, current%pardata%rp, current%pardata%bulk_conc, current%pardata%delNBbyCV
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
