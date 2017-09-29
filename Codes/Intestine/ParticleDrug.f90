!==================================================================================================
MODULE ParticleDrug				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================
USE SetPrecision
USE Setup
USE IC
USE MPI
USE Output

IMPLICIT NONE

CONTAINS


!===================================================================================================
SUBROUTINE Compute_C_bulk				  ! Computes the mesh-independent bulk concentration
!===================================================================================================
IMPLICIT NONE

INTEGER(lng)  		 :: i,j,k, mpierr,ID
INTEGER(lng)  		 :: ix0,ix1,iy0,iy1,iz0,iz00,iz1,iz11		! Trilinear interpolation parameters
REAL(dbl)     		 :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd		! Trilinear interpolation parameters
REAL(dbl)  	   	 :: xp,yp,zp
REAL(dbl)		 :: delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
REAL(dbl)       	 :: n_b
REAL(dbl)    	         :: R_P, Sh_P, delta_P
REAl(dbl)                :: R_influence_p, L_influence_p		! Parameters related to particle's volume of influence
REAl(dbl)                :: V_influence_P	 			! Parameters related to particle's volume of influence
REAL(dbl)                :: Delta_L
REAL(dbl)                :: x_DP, y_DP, z_DP				! Coordinates of "Discretized Point" (DP)
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

delta_mesh = 1.0_dbl
zcf3 = 1.0_dbl
n_b = 2.0

Cb_Total_Veff   = 0_dbl
Cb_Total_Veff_l = 0_dbl
NumFluids_Veff_l= 0
NumFluids_Veff  = 0


current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next
   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN	
      ID=current%pardata%parid
      R_P = current%pardata%rp   !Particle length scale: delta= R/Sh & effective radius: R_influence_P= R+(n_b*delta)
      delta_P= R_P 
      R_influence_P= (R_P+n_b*delta_P)/xcf
!---- Computing equivalent cubic mesh length scale
      V_influence_P= (4.0_dbl/3.0_dbl)*PI* R_influence_P**3.0_dbl
      L_influence_P= V_influence_P **(1.0_dbl/3.0_dbl)
      V_eff_Ratio  = V_influence_P/zcf3 					! Ratio of the effective volume to cell size 

!----------------------------------------------------------------------------------------------------------------------
!--Veff is smaller than the mesh volume --> Cb = Trilinear interpolation of the concentration at particle location
!--No communication is necessary between processors
!----------------------------------------------------------------------------------------------------------------------
      IF (V_eff_Ratio .LE. 3.375) THEN 					
         CaseNo= 1
         IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------- Finding local particle location (at current processor)
            xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
            yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
            zp= current%pardata%zp - REAL(kMin-1_lng,dbl)
            ix0 =FLOOR(xp)
            ix1 =CEILING(xp)
            iy0 =FLOOR(yp)
            iy1 =CEILING(yp)
            iz0 =FLOOR(zp)
            iz1 =CEILING(zp)
!---------- TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
            IF (ix1 /= ix0) THEN
               xd= (xp-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
            ELSE
               xd= 0.0_dbl
            END IF

            IF (iy1 /= iy0) THEN
               yd= (yp-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
            ELSE
               yd= 0.0_dbl
            END IF

            IF (iz1 /= iz0) THEN
               zd= (zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
            ELSE
               zd= 0.0_dbl
            END IF
!---------- Interpolation in x-direction
            c00 = phi(ix0,iy0,iz0) * (1.0_dbl-xd) + phi(ix1,iy0,iz0) * xd
            c01 = phi(ix0,iy0,iz1) * (1.0_dbl-xd) + phi(ix1,iy0,iz1) * xd
            c10 = phi(ix0,iy1,iz0) * (1.0_dbl-xd) + phi(ix1,iy1,iz0) * xd
            c11 = phi(ix0,iy1,iz1) * (1.0_dbl-xd) + phi(ix1,iy1,iz1) * xd
!---------- Interpolation in y-direction
            c0  = c00 * (1.0_dbl-yd) + c10 * yd
            c1  = c01 * (1.0_dbl-yd) + c11 * yd
!---------- Interpolation in z-direction
            c   = c0 * (1.0_dbl-zd) + c1 * zd
            Cb_Hybrid= c 
            current%pardata%bulk_conc = Cb_Hybrid
         END IF !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END IF 





!-----------------------------------------------------------------------------------------------------------------------
!-----Veff is larger than 1 which means parallel communication might be necessary 
!-----------------------------------------------------------------------------------------------------------------------
      IF (V_eff_Ratio .GT. 3.375) THEN
!------- Finding particle global location (in the whole domain not this particular processor)
         xp= current%pardata%xp 
         yp= current%pardata%yp 
         zp= current%pardata%zp 
!------- Global Volume of Influence Border (GVIB) for this particle (in the whole domain not this particular processor)
         GVIB_x(1)= xp - 0.5_dbl * L_influence_P
         GVIB_x(2)= xp + 0.5_dbl * L_influence_P
         GVIB_y(1)= yp - 0.5_dbl * L_influence_P
         GVIB_y(2)= yp + 0.5_dbl * L_influence_P
         GVIB_z(1)= zp - 0.5_dbl * L_influence_P
         GVIB_z(2)= zp + 0.5_dbl * L_influence_P

!------- TAKING CARE OF THE PERIODIC BC
         GVIB_z_Per(1) = GVIB_z(1)
         GVIB_z_Per(2) = GVIB_z(2)

         IF (GVIB_z(1) .LT. 1) THEN
            GVIB_z_Per(1) = GVIB_z(1) + nz
            GVIB_z_Per(2) = GVIB_z(2) + nz
         ENDIF
         IF (GVIB_z(2) .GT. nz) THEN
            GVIB_z_Per(1) = GVIB_z(1) - nz
            GVIB_z_Per(2) = GVIB_z(2) - nz
         ENDIF

!------- Finding processor that have overlap with effective volume around the particle       
300      IF ((((GVIB_x(1) .GE. (iMin-1_lng)) .AND. (GVIB_x(1) .LT. iMax)) .OR. ((GVIB_x(2) .GE. (iMin-1_lng)) .AND. (GVIB_x(2) .LT. iMax))) .AND. &
             (((GVIB_y(1) .GE. (jMin-1_lng)) .AND. (GVIB_y(1) .LT. jMax)) .OR. ((GVIB_y(2) .GE. (jMin-1_lng)) .AND. (GVIB_y(2) .LT. jMax))) .AND. &
             (((GVIB_z(1) .GE. (kMin-1_lng)) .AND. (GVIB_z(1) .LT. kMax)) .OR. ((GVIB_z(2) .GE. (kMin-1_lng)) .AND. (GVIB_z(2) .LT. kMax)))  )THEN

!---------- Veff is slightly larger than lattice cell volume --> Volume of influence is discretized to provide 27 points
!---------- Cb= Average of concentration interpolated on each of the descritized nodes inside volume of influence
            IF ((V_eff_Ratio .GT. 3.375) .AND. (V_eff_Ratio .LT. 8.0)) THEN		
               CaseNo= 2
!------------- Discretizing the volume of influence to  make sure at least 27 points are available
               Delta_L = (GVIB_x(2)-GVIB_x(1)) / 2.0 
!------------- Loop over discretized points and averaging the concentration
               DO i= 0, 2
                  DO j= 0, 2
                     DO k= 0, 2
                        x_DP= GVIB_x(1) + (i * Delta_L) 
                        y_DP= GVIB_y(1) + (j * Delta_L)
                        z_DP= GVIB_z(1) + (k * Delta_L)
                        IF((x_DP .GE. (REAL(iMin,dbl)-1.0_dbl)) .AND. &               
                           (x_DP .LT.  REAL(iMax,dbl)         ) .AND. &
                           (y_DP .GE. (REAL(jMin,dbl)-1.0_dbl)) .AND. & 
                           (y_DP .LT.  REAL(jMax,dbl)         ) .AND. &
                           (z_DP .GE. (REAL(kMin,dbl)-1.0_dbl)) .AND. & 
                           (z_DP .LT.  REAL(kMax,dbl)         )) THEN
!------------------------- Finding Local lattice nodes surrounding this point (This point is discretized and is not a lattice node))
                           ix0 = FLOOR(x_DP)   - (REAL(iMin,dbl)-1.0_dbl)
                           ix1 = CEILING(x_DP) - (REAL(iMin,dbl)-1.0_dbl)
                           iy0 = FLOOR(y_DP)   - (REAL(jMin,dbl)-1.0_dbl)
                           iy1 = CEILING(y_DP) - (REAL(jMin,dbl)-1.0_dbl) 
                           iz0 = FLOOR(z_DP)   - (REAL(kMin,dbl)-1.0_dbl)
                           iz1 = CEILING(z_DP) - (REAL(kMin,dbl)-1.0_dbl)
                              
                           x_DP = x_DP - REAL(iMin-1_lng,dbl)
                           y_DP = y_DP - REAL(jMin-1_lng,dbl)
                           z_DP = z_DP - REAL(kMin-1_lng,dbl)
 
!------------------------- TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
                           IF (ix1 /= ix0) THEN
                              xd=(x_DP-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
                           ELSE
                              xd=0.0_dbl
                           END IF
                           IF (iy1 /= iy0) THEN
                              yd=(y_DP-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
                           ELSE
                              yd=0.0_dbl
                           END IF
                           IF (iz1 /= iz0) THEN
                              zd=(z_DP-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
                           ELSE
                              zd=0.0_dbl
                           END IF

!------------------------- Interpolation in x-direction
                           c00 = phi(ix0,iy0,iz0) * (1.0_dbl-xd) + phi(ix1,iy0,iz0) * xd
                           c01 = phi(ix0,iy0,iz1) * (1.0_dbl-xd) + phi(ix1,iy0,iz1) * xd
                           c10 = phi(ix0,iy1,iz0) * (1.0_dbl-xd) + phi(ix1,iy1,iz0) * xd
                           c11 = phi(ix0,iy1,iz1) * (1.0_dbl-xd) + phi(ix1,iy1,iz1) * xd
!------------------------- Interpolation in y-direction
                           c0  = c00 * (1.0_dbl-yd) + c10 * yd
                           c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------------------------- Interpolation in z-direction
                           c   = c0 * (1.0_dbl-zd) + c1 * zd

                           Cb_Total_Veff_l(ID)  = Cb_Total_Veff_l(ID)  + c
                           NumFluids_Veff_l(ID) = NumFluids_Veff_l(ID) + 1_lng
                        END IF
                     END DO
                  END DO
               END DO
          
!----------------------------------------------------------------------------------------------------------------------
!---------- Veff is much larger than mesh volume --> Cb= total number of moles in volume of influence / volume of influence 
!----------------------------------------------------------------------------------------------------------------------
            ELSE IF (V_eff_Ratio .GE. 8.0) THEN                             
               CaseNo = 3
!------------- Finding the lattice "Nodes Effected by Particle"
               GNEP_x(1)= CEILING(GVIB_x(1))
               GNEP_y(1)= CEILING(GVIB_y(1))
               GNEP_z(1)= CEILING(GVIB_z(1))
               GNEP_x(2)= FLOOR  (GVIB_x(2))
               GNEP_y(2)= FLOOR  (GVIB_y(2))
               GNEP_z(2)= FLOOR  (GVIB_z(2))
!------------- Finding the lattice "Nodes Effected by Particle"
               NEP_x(1)= Max(GNEP_x(1),iMin)- (iMin-1)
               NEP_y(1)= Max(GNEP_y(1),jMin)- (jMin-1)
               NEP_z(1)= Max(GNEP_z(1),kMin)- (kMin-1)
               NEP_x(2)= Min(GNEP_x(2),iMax)- (iMin-1)
               NEP_y(2)= Min(GNEP_y(2),jMax)- (jMin-1)
               NEP_z(2)= Min(GNEP_z(2),kMax)- (kMin-1)

               DO i= NEP_x(1),NEP_x(2) 
                  DO j= NEP_y(1),NEP_y(2)
                     DO k= NEP_z(1),NEP_z(2)
                        IF (node(i,j,k) .EQ. FLUID) THEN
                           Cb_Total_Veff_l(ID)  = Cb_Total_Veff_l(ID)  + phi(i,j,k)
                           NumFluids_Veff_l(ID) = NumFluids_Veff_l(ID) + 1_lng
                        END IF
                     END DO
                  END DO
               END DO
            END IF                   ! Conditional for cases 2 and 3
         END IF 								      ! Conditional for the processor which has overlap with effective volume 

!------- TAKING CARE OF THE PERIODIC BC
         IF (GVIB_z_Per(1) .NE. GVIB_z(1)) THEN
            GVIB_z(1)= GVIB_z_Per(1) 
            GVIB_z(2)= GVIB_z_Per(2) 
            GO TO 300
         ENDIF
      ENDIF
   ENDIF   
   current => next
END DO

!----Communication with other processors for V_eff greater than 1
CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
CALL MPI_ALLREDUCE(Cb_Total_Veff_l , Cb_Total_Veff , np, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
CALL MPI_ALLREDUCE(NumFluids_Veff_l, NumFluids_Veff, np, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, mpierr)

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next
   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN	
      ID= current%pardata%parid
      IF (NumFluids_Veff(ID) .GE. 1) THEN 
         Cb_Hybrid= Cb_Total_Veff(ID) / NumFluids_Veff(ID)
         current%pardata%bulk_conc = Cb_Hybrid
      END IF   
   END IF       			                                    		!End of conditional for V_eff greater than 1 
   current => next
END DO
!===================================================================================================
END SUBROUTINE Compute_C_bulk
!===================================================================================================







!===================================================================================================
SUBROUTINE Particle_Drug_Release                     ! Calculates drug release at every time step  
!===================================================================================================
IMPLICIT NONE
INTEGER(lng)  :: numFluids,i,j,k,RANK,mpierr
REAL(dbl)     :: deltaR,temp,cbt,zcf3,bulkconc 
REAL(dbl)     :: Sherwood
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

zcf3=xcf*ycf*zcf

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (mySub .EQ.current%pardata%cur_part) THEN
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                    !only calculate the drug release when particle radius is larger than 0.1 micron
         current%pardata%rpold = current%pardata%rp
         bulkconc = current%pardata%bulk_conc
         Sherwood = 1.0_dbl + current%pardata%sh_conf + current%pardata%sh_shear + current%pardata%sh_slip 
         temp = current%pardata%rpold**2.0_dbl-4.0_dbl*tcf*molarvol*diffm*Sherwood*max((current%pardata%par_conc-bulkconc),0.0_dbl)
         IF (temp.GE.0.0_dbl) THEN
            current%pardata%rp= 0.5_dbl*(current%pardata%rpold+sqrt(temp))
         ELSE
            current%pardata%rp= 0.5_dbl*(current%pardata%rpold)
         END IF
         deltaR=current%pardata%rpold-current%pardata%rp
         current%pardata%delNBbyCV = (4.0_dbl/3.0_dbl) * PI*(current%pardata%rpold**3.0_dbl - current%pardata%rp**3.0_dbl) /(molarvol*zcf3)

      ELSE IF ((current%pardata%rp .LE. Min_R_Acceptable) .AND. (current%pardata%rp .NE. 0.0)) THEN
         current%pardata%xp        = 0.0_dbl
         current%pardata%yp        = 0.0_dbl
         current%pardata%zp        = 0.0_dbl
         current%pardata%up        = 0.0_dbl
         current%pardata%vp        = 0.0_dbl
         current%pardata%wp        = 0.0_dbl
         current%pardata%rp        = 0.0_dbl
         current%pardata%delNBbyCV = 0.0_dbl
         current%pardata%bulk_conc = 0.0_dbl
         current%pardata%S		     = 0.0_Dbl
         current%pardata%Sst       = 0.0_Dbl
         current%pardata%sh_conf   = 0.0_Dbl
         current%pardata%sh_shear  = 0.0_Dbl
         current%pardata%sh_slip   = 0.0_Dbl
      END IF
   END IF 
   current => next
ENDDO

!===================================================================================================
END SUBROUTINE Particle_Drug_Release
!===================================================================================================







!===================================================================================================
SUBROUTINE Compute_Sherwood 
!===================================================================================================
! Called by Particle_Track (LBM.f90) used in Particle_Drug_Release 
! Incporates hierarchical mdoel to Sh(t) to include effect of shear/hydrodynamics and container effect

IMPLICIT NONE
INTEGER(lng)  :: mpierr, RANK
REAL(dbl)     :: S,Sst,Sh0 
REAL(dbl)     :: xp, yp, zp
REAL(dbl)     :: Re_P_du, Pe_P_du   !Reynolds and Peclet numbers for each particle based on the slip velocity
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (mySub .EQ.current%pardata%cur_part) THEN
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                     ! only calculate the drug release when particle radius is larger than 0.1 micron
        current%pardata%sh_conf = 0.0_dbl
        current%pardata%sh_shear= 0.0_dbl
        current%pardata%sh_slip = 0.0_dbl
         !----- If including the confinement effects ----------------------------------------------- 
         IF (Flag_Confinement_Effects) THEN                                 
            current%pardata%sh_conf= (current%pardata%gamma_cont / (1.0_dbl-current%pardata%gamma_cont)) 
         ENDIF

         !----- If including hydrodynamic shear -----------------------------------------------------          
         IF (Flag_Shear_Effects) THEN            
            S= current%pardata%S
            Sst= S* (current%pardata%rp**2.0) / diffm
            current%pardata%Sst= Sst
            IF (Sst .LE. 5.0_dbl) THEN
               current%pardata%sh_shear = 1.0_dbl + 0.281_dbl*(Sst**0.5_dbl)          -1.0_dbl
            ELSE IF ((Sst .GT. 5.0_dbl).AND.(Sst .LE. 80.0)) THEN     
               current%pardata%sh_shear = 1.181_dbl*(Sst**0.2_dbl)                    -1.0_dbl
            ELSE IF (Sst.GT.80.0) THEN
               current%pardata%sh_shear = 4.5_dbl - (7.389/(Sst**(1.0_dbl/3.0_dbl)) ) -1.0_dbl 
            ENDIF
         ENDIF 

         !----- If including hydrodynamic convection ------------------------------------------------          
         IF (Flag_Convection_Effects) THEN            
            Re_P_du= current%pardata%U_slip* current%pardata%rp / nu                          !Computing Reynolds number
            Pe_P_du= current%pardata%U_slip* current%pardata%rp / diffm                         !Computing Peclet number
            current%pardata%sh_slip = 0.424* (Pe_P_du**0.33) * (Re_P_du**0.17)     
            !write(*,*) 'iter,ID,nu,Dm,Re,Pe,Sh',iter,current%pardata%parid,nu,diffm,Re_P_du,Pe_P_du,current%pardata%sh_slip  
         ENDIF
      ENDIF
   ENDIF
   current => next
ENDDO
!===================================================================================================
END SUBROUTINE Compute_Sherwood
!===================================================================================================




!===================================================================================================
SUBROUTINE Compute_shear                
!===================================================================================================
! Computing components of velocity gradietns
! and strain rate tensor using central difference 

IMPLICIT NONE

INTEGER(lng)  :: i,j,k
INTEGER(lng)  :: ii,jj,kk
INTEGER(lng)  :: ix0,iy0,iz0
INTEGER(lng)  :: ix1,iy1,iz1

REAL(dbl)     :: xaxis,yaxis
REAL(dbl)     :: X_s,Y_s,R_s
REAL(dbl)     :: CosTheta_s,SinTheta_s
REAL(dbl)     :: Vel_s

REAL(dbl)     :: xp,yp,zp
REAL(dbl)     :: xd,yd,zd

REAL(dbl)     :: dudx00, dudx01, dudx10, dudx11
REAL(dbl)     :: dudx0, dudx1, dudx_P
REAL(dbl)     :: dudy00, dudy01, dudy10, dudy11
REAL(dbl)     :: dudy0, dudy1, dudy_P
REAL(dbl)     :: dudz00, dudz01, dudz10, dudz11
REAL(dbl)     :: dudz0, dudz1, dudz_P

REAL(dbl)     :: dvdx00, dvdx01, dvdx10, dvdx11
REAL(dbl)     :: dvdx0, dvdx1, dvdx_P
REAL(dbl)     :: dvdy00, dvdy01, dvdy10, dvdy11
REAL(dbl)     :: dvdy0, dvdy1, dvdy_P
REAL(dbl)     :: dvdz00, dvdz01, dvdz10, dvdz11
REAL(dbl)     :: dvdz0, dvdz1, dvdz_P

REAL(dbl)     :: dwdx00, dwdx01, dwdx10, dwdx11
REAL(dbl)     :: dwdx0, dwdx1, dwdx_P
REAL(dbl)     :: dwdy00, dwdy01, dwdy10, dwdy11
REAL(dbl)     :: dwdy0, dwdy1, dwdy_P
REAL(dbl)     :: dwdz00, dwdz01, dwdz10, dwdz11
REAL(dbl)     :: dwdz0, dwdz1, dwdz_P

REAL(dbl)     :: E11, E12, E13, E21, E22, E23, E31, E32, E33
REAL(dbl)     :: S

TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next

   IF (mySub .EQ. current%pardata%cur_part) THEN
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN						!only when particle radius is larger than 0.1 micron				
         xp = current%pardata%xp - REAL(iMin-1_lng,dbl)
         yp = current%pardata%yp - REAL(jMin-1_lng,dbl)
         zp = current%pardata%zp - REAL(kMin-1_lng,dbl)

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
            u_m= u 
            v_m= v 
            w_m= w
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
                        u_m(ii,jj,kk) = Vel_s * CosTheta_s
                        v_m(ii,jj,kk) = Vel_s * SinTheta_s
                        w_m(ii,jj,kk) = 0.0_dbl				
                     ELSE 									! Fluid nodes in the lattice cell encompassing the particle
                        u_m(ii,jj,kk) = u(ii,jj,kk) 
                        v_m(ii,jj,kk) = v(ii,jj,kk) 
                        w_m(ii,jj,kk) = w(ii,jj,kk)
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

!------- dudx -------------------------------------
         dudx00= u_m(ix1,iy0,iz0)- u_m(ix0,iy0,iz0)	
         dudx01= u_m(ix1,iy0,iz1)- u_m(ix0,iy0,iz1)	
         dudx10= u_m(ix1,iy1,iz0)- u_m(ix0,iy1,iz0)	
         dudx11= u_m(ix1,iy1,iz1)- u_m(ix0,iy1,iz1)
!------- y-dir interpolation ----------------------
         dudx0 = dudx00*(1.0_dbl-yd) + dudx10* yd
         dudx1 = dudx01*(1.0_dbl-yd) + dudx11* yd
!------- z-dir interpolation ----------------------
         dudx_P  = dudx0*(1.0_dbl-zd)+dudx1*zd
!--------------------------------------------------

!------- dvdx -------------------------------------
         dvdx00= v_m(ix1,iy0,iz0)- v_m(ix0,iy0,iz0)	
         dvdx01= v_m(ix1,iy0,iz1)- v_m(ix0,iy0,iz1)	
         dvdx10= v_m(ix1,iy1,iz0)- v_m(ix0,iy1,iz0)	
         dvdx11= v_m(ix1,iy1,iz1)- v_m(ix0,iy1,iz1)
!------- y-dir interpolation ----------------------
         dvdx0 = dvdx00*(1.0_dbl-yd) + dvdx10* yd
         dvdx1 = dvdx01*(1.0_dbl-yd) + dvdx11* yd
!------- z-dir interpolation ----------------------
         dvdx_P  = dvdx0*(1.0_dbl-zd)+dvdx1*zd
!--------------------------------------------------

!------- dwdx -------------------------------------
         dwdx00= w_m(ix1,iy0,iz0)- w_m(ix0,iy0,iz0)	
         dwdx01= w_m(ix1,iy0,iz1)- w_m(ix0,iy0,iz1)	
         dwdx10= w_m(ix1,iy1,iz0)- w_m(ix0,iy1,iz0)	
         dwdx11= w_m(ix1,iy1,iz1)- w_m(ix0,iy1,iz1)
!------- y-dir interpolation ----------------------
         dwdx0 = dwdx00*(1.0_dbl-yd) + dwdx10* yd
         dwdx1 = dwdx01*(1.0_dbl-yd) + dwdx11* yd
!------- z-dir interpolation ----------------------
         dwdx_P  = dwdx0*(1.0_dbl-zd)+dwdx1*zd
!--------------------------------------------------

!------- dudy -------------------------------------
         dudy00= u_m(ix0,iy1,iz0)- u_m(ix0,iy0,iz0)	
         dudy01= u_m(ix0,iy1,iz1)- u_m(ix0,iy0,iz1)	
         dudy10= u_m(ix1,iy1,iz0)- u_m(ix1,iy0,iz0)	
         dudy11= u_m(ix1,iy1,iz1)- u_m(ix1,iy0,iz1)
!------- x-dir interpolation ----------------------
         dudy0 = dudy00*(1.0_dbl-xd) + dudy10* xd
         dudy1 = dudy01*(1.0_dbl-xd) + dudy11* xd
!------- z-dir interpolation ----------------------
         dudy_P= dudy0*(1.0_dbl-zd)+dudy1*zd
!--------------------------------------------------

!------- dvdy -------------------------------------
         dvdy00= v_m(ix0,iy1,iz0)- v_m(ix0,iy0,iz0)	
         dvdy01= v_m(ix0,iy1,iz1)- v_m(ix0,iy0,iz1)	
         dvdy10= v_m(ix1,iy1,iz0)- v_m(ix1,iy0,iz0)	
         dvdy11= v_m(ix1,iy1,iz1)- v_m(ix1,iy0,iz1)
!------- x-dir interpolation ----------------------
         dvdy0 = dvdy00*(1.0_dbl-xd) + dvdy10* xd
         dvdy1 = dvdy01*(1.0_dbl-xd) + dvdy11* xd
!------- z-dir interpolation ----------------------
         dvdy_P= dvdy0*(1.0_dbl-zd)+dvdy1*zd
!--------------------------------------------------

!------- dwdy -------------------------------------
         dwdy00= w_m(ix0,iy1,iz0)- w_m(ix0,iy0,iz0)	
         dwdy01= w_m(ix0,iy1,iz1)- w_m(ix0,iy0,iz1)	
         dwdy10= w_m(ix1,iy1,iz0)- w_m(ix1,iy0,iz0)	
         dwdy11= w_m(ix1,iy1,iz1)- w_m(ix1,iy0,iz1)
!------- x-dir interpolation ----------------------
         dwdy0 = dwdy00*(1.0_dbl-xd) + dwdy10* xd
         dwdy1 = dwdy01*(1.0_dbl-xd) + dwdy11* xd
!------- z-dir interpolation ----------------------
         dwdy_P= dwdy0*(1.0_dbl-zd)+dwdy1*zd
!--------------------------------------------------

!------- dudz -------------------------------------
         dudz00= u_m(ix0,iy0,iz1)- u_m(ix0,iy0,iz0)	
         dudz01= u_m(ix0,iy1,iz1)- u_m(ix0,iy1,iz0)	
         dudz10= u_m(ix1,iy0,iz1)- u_m(ix1,iy0,iz0)	
         dudz11= u_m(ix1,iy1,iz1)- u_m(ix1,iy1,iz0)
!------- x-dir interpolation ----------------------
         dudz0 = dudz00*(1.0_dbl-xd) + dudz10* xd
         dudz1 = dudz01*(1.0_dbl-xd) + dudz11* xd
!------- y-dir interpolation ----------------------
         dudz_P= dudz0*(1.0_dbl-yd)+dudz1*yd
!--------------------------------------------------

!------- dvdz -------------------------------------
         dvdz00= v_m(ix0,iy0,iz1)- v_m(ix0,iy0,iz0)	
         dvdz01= v_m(ix0,iy1,iz1)- v_m(ix0,iy1,iz0)	
         dvdz10= v_m(ix1,iy0,iz1)- v_m(ix1,iy0,iz0)	
         dvdz11= v_m(ix1,iy1,iz1)- v_m(ix1,iy1,iz0)
!------- x-dir interpolation ----------------------
         dvdz0 = dvdz00*(1.0_dbl-xd) + dvdz10* xd
         dvdz1 = dvdz01*(1.0_dbl-xd) + dvdz11* xd
!------- y-dir interpolation ----------------------
         dvdz_P= dvdz0*(1.0_dbl-yd)+dvdz1*yd
!--------------------------------------------------

!------- dwdz -------------------------------------
         dwdz00= w_m(ix0,iy0,iz1)- w_m(ix0,iy0,iz0)	
         dwdz01= w_m(ix0,iy1,iz1)- w_m(ix0,iy1,iz0)	
         dwdz10= w_m(ix1,iy0,iz1)- w_m(ix1,iy0,iz0)	
         dwdz11= w_m(ix1,iy1,iz1)- w_m(ix1,iy1,iz0)
!------- x-dir interpolation ----------------------
         dwdz0 = dwdz00*(1.0_dbl-xd) + dwdz10* xd
         dwdz1 = dwdz01*(1.0_dbl-xd) + dwdz11* xd
!------- y-dir interpolation ----------------------
         dwdz_P= dwdz0*(1.0_dbl-yd)+dwdz1*yd
!--------------------------------------------------

!======= Computing 9 componenets of the strain rate tensor
!         E11= dudx_P
!         E12= 0.5_dbl*(dudy_P+dvdx_P)
!         E13= 0.5_dbl*(dudz_P+dwdx_P)
!
!         E21= 0.5_dbl*(dvdx_P+dudy_P)
!         E22= dvdy_P
!         E23= 0.5_dbl*(dvdz_P+dwdy_P)
!
!         E31= 0.5_dbl*(dwdx_P+dudz_P)
!         E32= 0.5_dbl*(dwdy_P+dvdz_P)
!         E33= dwdz_P

!======= Computing the strain rate magnitude
         S = 2.0_dbl*dudx_P**2.0_dbl  + &
             2.0_dbl*dvdy_P**2.0_dbl  + &
             2.0_dbl*dwdz_P**2.0_dbl  + &
             (dudy_P+dvdx_P)**2.0_dbl + &
             (dudz_P+dwdx_P)**2.0_dbl + &
             (dwdy_P+dvdz_P)**2.0_dbl
         S = sqrt(S) 
         S = S * vcf / zcf 
         current%pardata%S = S
       END IF 
    END IF
    current => next
ENDDO

!===================================================================================================
END SUBROUTINE Compute_shear
!===================================================================================================





!===================================================================================================
SUBROUTINE Compute_U_slip
!===================================================================================================
IMPLICIT NONE

TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next
INTEGER(lng):: ix0,ix1,iy0,iy1,iz0,iz00,iz1,iz11		! Trilinear interpolation parameters
REAL(dbl)  	:: xp,yp,zp
REAL(dbl)   :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd		! Trilinear interpolation parameters
REAL(dbl)   :: U_slip_x,U_slip_y, U_slip_z
REAL(dbl)   :: alpha_P
REAL(dbl)   :: R_P
REAL(dbl)   :: Laplacian_x_P, Laplacian_y_P, Laplacian_z_P
REAL(dbl)   :: DLaplacianDt_x_P, DLaplacianDt_y_P, DLaplacianDt_z_P
REAL(dbl)   :: DUDt_x_P, DUDt_y_P, DUDt_z_P
REAL(dbl)   :: A_P, B_P, C_P

alpha_P= den/(den_P+0.5_dbl*den)

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (mySub .EQ.current%pardata%cur_part) THEN
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN 
         R_P= current%pardata%rp 
!------- Finding local particle location (at current processor)
         xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
         yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
         zp= current%pardata%zp - REAL(kMin-1_lng,dbl)
         ix0 =FLOOR(xp)
         ix1 =CEILING(xp)
         iy0 =FLOOR(yp)
         iy1 =CEILING(yp)
         iz0 =FLOOR(zp)
         iz1 =CEILING(zp)
!------- TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES !!!!!!!!!!!!!!!!!!!!!!!11
         IF (ix1 /= ix0) THEN
            xd= (xp-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
         ELSE
            xd= 0.0_dbl
         END IF

         IF (iy1 /= iy0) THEN
            yd= (yp-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
         ELSE
            yd= 0.0_dbl
         END IF

         IF (iz1 /= iz0) THEN
            zd= (zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
         ELSE
            zd= 0.0_dbl
         END IF

!======= Laplacian_x: interpolation to the location of the particle ====================
!------- x-direction -----------------------------------------------------------------
         c00 = Laplacian_x(ix0,iy0,iz0) * (1.0_dbl-xd) + Laplacian_x(ix1,iy0,iz0) * xd
         c01 = Laplacian_x(ix0,iy0,iz1) * (1.0_dbl-xd) + Laplacian_x(ix1,iy0,iz1) * xd
         c10 = Laplacian_x(ix0,iy1,iz0) * (1.0_dbl-xd) + Laplacian_x(ix1,iy1,iz0) * xd
         c11 = Laplacian_x(ix0,iy1,iz1) * (1.0_dbl-xd) + Laplacian_x(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         Laplacian_x_P = c0 * (1.0_dbl-zd) + c1 * zd

!------- x-direction -----------------------------------------------------------------
         c00 = Laplacian_y(ix0,iy0,iz0) * (1.0_dbl-xd) + Laplacian_y(ix1,iy0,iz0) * xd
         c01 = Laplacian_y(ix0,iy0,iz1) * (1.0_dbl-xd) + Laplacian_y(ix1,iy0,iz1) * xd
         c10 = Laplacian_y(ix0,iy1,iz0) * (1.0_dbl-xd) + Laplacian_y(ix1,iy1,iz0) * xd
         c11 = Laplacian_y(ix0,iy1,iz1) * (1.0_dbl-xd) + Laplacian_y(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         Laplacian_y_P = c0 * (1.0_dbl-zd) + c1 * zd

!------- x-direction -----------------------------------------------------------------
         c00 = Laplacian_z(ix0,iy0,iz0) * (1.0_dbl-xd) + Laplacian_z(ix1,iy0,iz0) * xd
         c01 = Laplacian_z(ix0,iy0,iz1) * (1.0_dbl-xd) + Laplacian_z(ix1,iy0,iz1) * xd
         c10 = Laplacian_z(ix0,iy1,iz0) * (1.0_dbl-xd) + Laplacian_z(ix1,iy1,iz0) * xd
         c11 = Laplacian_z(ix0,iy1,iz1) * (1.0_dbl-xd) + Laplacian_z(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         Laplacian_z_P = c0 * (1.0_dbl-zd) + c1 * zd


!======= DUDt: interpolation to the location of the particle ===========================
!------- x-direction
         c00 = DUDt_x(ix0,iy0,iz0) * (1.0_dbl-xd) + DUDt_x(ix1,iy0,iz0) * xd
         c01 = DUDt_x(ix0,iy0,iz1) * (1.0_dbl-xd) + DUDt_x(ix1,iy0,iz1) * xd
         c10 = DUDt_x(ix0,iy1,iz0) * (1.0_dbl-xd) + DUDt_x(ix1,iy1,iz0) * xd
         c11 = DUDt_x(ix0,iy1,iz1) * (1.0_dbl-xd) + DUDt_x(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         DUDt_x_P = c0 * (1.0_dbl-zd) + c1 * zd

!------- x-direction
         c00 = DUDt_y(ix0,iy0,iz0) * (1.0_dbl-xd) + DUDt_y(ix1,iy0,iz0) * xd
         c01 = DUDt_y(ix0,iy0,iz1) * (1.0_dbl-xd) + DUDt_y(ix1,iy0,iz1) * xd
         c10 = DUDt_y(ix0,iy1,iz0) * (1.0_dbl-xd) + DUDt_y(ix1,iy1,iz0) * xd
         c11 = DUDt_y(ix0,iy1,iz1) * (1.0_dbl-xd) + DUDt_y(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         DUDt_y_P = c0 * (1.0_dbl-zd) + c1 * zd

!------- x-direction
         c00 = DUDt_z(ix0,iy0,iz0) * (1.0_dbl-xd) + DUDt_z(ix1,iy0,iz0) * xd
         c01 = DUDt_z(ix0,iy0,iz1) * (1.0_dbl-xd) + DUDt_z(ix1,iy0,iz1) * xd
         c10 = DUDt_z(ix0,iy1,iz0) * (1.0_dbl-xd) + DUDt_z(ix1,iy1,iz0) * xd
         c11 = DUDt_z(ix0,iy1,iz1) * (1.0_dbl-xd) + DUDt_z(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         DUDt_z_P = c0 * (1.0_dbl-zd) + c1 * zd



!======= DLaplacianDt: interpolation to the location of the particle ====================
!------- x-direction
         c00 = DLaplacianDt_x(ix0,iy0,iz0) * (1.0_dbl-xd) + DLaplacianDt_x(ix1,iy0,iz0) * xd
         c01 = DLaplacianDt_x(ix0,iy0,iz1) * (1.0_dbl-xd) + DLaplacianDt_x(ix1,iy0,iz1) * xd
         c10 = DLaplacianDt_x(ix0,iy1,iz0) * (1.0_dbl-xd) + DLaplacianDt_x(ix1,iy1,iz0) * xd
         c11 = DLaplacianDt_x(ix0,iy1,iz1) * (1.0_dbl-xd) + DLaplacianDt_x(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         DLaplacianDt_x_P = c0 * (1.0_dbl-zd) + c1 * zd

!------- x-direction
         c00 = DLaplacianDt_y(ix0,iy0,iz0) * (1.0_dbl-xd) + DLaplacianDt_y(ix1,iy0,iz0) * xd
         c01 = DLaplacianDt_y(ix0,iy0,iz1) * (1.0_dbl-xd) + DLaplacianDt_y(ix1,iy0,iz1) * xd
         c10 = DLaplacianDt_y(ix0,iy1,iz0) * (1.0_dbl-xd) + DLaplacianDt_y(ix1,iy1,iz0) * xd
         c11 = DLaplacianDt_y(ix0,iy1,iz1) * (1.0_dbl-xd) + DLaplacianDt_y(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         DLaplacianDt_y_P = c0 * (1.0_dbl-zd) + c1 * zd

!------- x-direction
         c00 = DLaplacianDt_z(ix0,iy0,iz0) * (1.0_dbl-xd) + DLaplacianDt_z(ix1,iy0,iz0) * xd
         c01 = DLaplacianDt_z(ix0,iy0,iz1) * (1.0_dbl-xd) + DLaplacianDt_z(ix1,iy0,iz1) * xd
         c10 = DLaplacianDt_z(ix0,iy1,iz0) * (1.0_dbl-xd) + DLaplacianDt_z(ix1,iy1,iz0) * xd
         c11 = DLaplacianDt_z(ix0,iy1,iz1) * (1.0_dbl-xd) + DLaplacianDt_z(ix1,iy1,iz1) * xd
!------- y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- z-direction
         DLaplacianDt_z_P = c0 * (1.0_dbl-zd) + c1 * zd


!======= Computing slip velocity =======================================================
         A_P= (R_P**2.0_dbl) / 6.0_dbl
         B_P= ((R_P**4.0_dbl)/(27.0_dbl*alpha_P*nu)) *(1.0_dbl-(3.0_dbl*alpha_P)/10.0_dbl)
         C_P= (2.0_dbl*(R_P**2.0_dbl)/(9.0_dbl*alpha_P*nu))*(1.0_dbl-(3.0_dbl*alpha_P/2.0_dbl))
         U_slip_x= A_P*Laplacian_x_P + B_P*DLaplacianDt_x_P - C_P*Dudt_x_P
         U_slip_y= A_P*Laplacian_y_P + B_P*DLaplacianDt_y_P - C_P*Dudt_y_P
         U_slip_z= A_P*Laplacian_z_P + B_P*DLaplacianDt_z_P - C_P*Dudt_z_P
         current%pardata%U_Slip_x= U_slip_x
         current%pardata%U_Slip_y= U_slip_y
         current%pardata%U_Slip_z= U_slip_z
         current%pardata%U_Slip = sqrt(U_slip_x**2.0_dbl + U_slip_y**2.0_dbl + U_slip_z**2.0_dbl)
        !IF(current%pardata%parid .EQ.1) THEN 
        !   write(*,*) iter,current%pardata%parid  
        !   write(*,*) iter, alpha_P,A_P,B_P,C_P
        !   write(*,*) iter, Laplacian_x_P,Laplacian_y_P,Laplacian_z_P   
        !   write(*,*) iter, DLaplacianDt_x_P,DLaplacianDt_y_P,DLaplacianDt_z_P   
        !   write(*,*) iter, DUDt_x_P,DUDt_y_P,DUDt_z_P  
        !   write(*,*) '------------------------------------------------------'
        !ENDIF
      ENDIF
   END IF  
   current => next
END DO   
!===================================================================================================
END SUBROUTINE Compute_U_slip
!===================================================================================================






!===============================================================================================================
SUBROUTINE Compute_C_surface_new !Calculating surface pH of drug with bicarbonate buffer using the IRR approach
!===============================================================================================================
IMPLICIT NONE

REAL(dbl)                :: c0,c1,c2,c3,c4,c5,c6
REAL(dbl)                :: S_ratio
REAL(dbl)                :: R_P,Sh_P,delta_P
REAL(dbl)                :: HA_o,pKaa,Kaa,Dha,Da,pKab,Kab,pKa1,Ka1
REAL(dbl)                :: Hh,Kw,OHh,BHh,Btotal,Dh,Doh,Dbh,Db,kd
REAL(dbl)                :: h,p,q,r,s,t,yy,yy_min,pH_s,proton,proton_bulk
INTEGER(lng)             :: ii
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next
proton_bulk=10.0**(-pH_bulk)
HA_o= S_intrinsic /1000.0_dbl     ! Intrinsic solubility of ibuprofen in mole/L (=0.00033)
pKaa= 4.43                        ! pKaa = the pka of the drug (ibuprofen)
Kaa=  10**(-pKaa)

!--- Dha and Da are the diffusion coefficients of the unionized and ionized form of the drug---------
Dha=7.93e-6
Da= 7.93e-6

!--- Bicarbonate Ka properties ---------------------------------------------------------------------
pKab= 3.55                        ! pka of bicarbonate used in the diffusion layer
Kab=  10**(-pKab) 
pKa1= 6.04                        ! pKa1 = the pka of of bicarbonate in the bulk solution
Ka1=  10**(-pKa1)

!--- Bicarbonate Buffer properties -----------------------------------------------------------------
Hh=10**(-pH_bulk)
Kw=2.57e-14
OHh=Kw/Hh
BHh=(Hh*Bh)/(10**(-3.55))          ! BHh = the concentration of H2CO3 in the bulk solution
Btotal=Bh+BHh
Dh=104.9e-6                       ! Dh = diffusion coefficient of H+
Doh=63.0e-6                       ! Doh= diffusion coefficient of OH-
Dbh=14.6e-6                       ! Dbh= diffusion coefficient of H2CO3
Db= 14.6e-6                       ! Db = diffusion coefficient of HCO3-
kd=50                             ! kd = dehydration rate constant for the IRR model

OPEN(131,FILE='Farhad_Debug.dat')
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (mySub .EQ.current%pardata%cur_part) THEN
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN 
         IF (Flag_Buffer) THEN 
            R_P  = 1000000.0* current%pardata%rp                        ! particle radius (um)
            Sh_P = 1.0_dbl+ current%pardata%sh_conf + current%pardata%sh_shear + current%pardata%sh_slip
            delta_P = (R_P / Sh_P)                                      ! diffusion layer thickness (um)
            !--- find roots of pH^3+qH^2+rH+s=0 to get the concentration of H+ ions at the particle's surface ------------------------
            h= delta_P/10000.0_dbl                                      ! diffusion layer thickness (cm)  
            p=1.0e16*(Dh*h*(sqrt(Dbh*kd)))
            q=1.0e16*((Kab*Db*Dh)-(h*(sqrt(Dbh*kd))*Dh*Hh)+(h*(sqrt(Dbh*kd))*Doh*OHh) +(h*(sqrt(Db*kd))*Db*Bh))
            r=1.0e16*((-Kab*Db*Dh*Hh)+(Kab*Db*Doh*OHh)-(h*(sqrt(Dbh*kd))*Doh*Kw) -(h*(sqrt(Dbh*kd))*Da*HA_o*Kaa) -(Kab*(h*((sqrt(Db*kd))*Db*BHh)))-(Db**(2)*Kab*Bh)+(Db**(2)*Kab*Bh))
            s=1.0e16*((-Kab*Db*Doh*Kw)-(Kab*Db*Da*HA_o*Kaa))
            yy_min=1.0e16
!           Do ii=1,10000000
!              proton= proton_bulk + ii*1e-8 ! /10000000.0
!              pH_s=-log10(proton) 
            Do ii=1,6000
               pH_s=pH_bulk- ii*0.001
               proton= 10.0_dbl**(-pH_s)
               yy= p*proton**3.0 + q*proton**2 + r*proton + s
               IF (yy.GE.0) THEN
                 GOTO 100
               ENDIF
            ENDDO
100         WRITE(*,*) 'pH_s',pH_s 
         ENDIF   
      END IF  
   END IF  
   current => next
END DO   

!===================================================================================================
END SUBROUTINE Compute_C_surface_new
!===================================================================================================





!===================================================================================================
SUBROUTINE Compute_C_surface
!===================================================================================================
IMPLICIT NONE

REAL(dbl)                :: c0,c1,c2,c3,c4,c5,c6
REAL(dbl)                :: S_ratio
REAL(dbl)                :: R_P,Sh_P,delta_P
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (mySub .EQ.current%pardata%cur_part) THEN
      IF (current%pardata%rp .GT. Min_R_Acceptable) THEN 

         IF (Flag_Buffer) THEN !--- Buffer Capacity =10.5 mM ------------------------------------------------
            R_P  = 1000000.0* current%pardata%rp   !units in  micron
            Sh_P = 1.0_dbl+ current%pardata%sh_conf + current%pardata%sh_shear + current%pardata%sh_slip
            delta_P = (R_P / Sh_P)                 !units in microns
            IF (delta_P .LE. 50.0) THEN
               c6= -0.000000002910474984
               c5=  0.000000504714268975
               c4= -0.000035024962744930
               c3=  0.001258840857793950
               c2= -0.026259477974747900
               c1=  0.435387995939737000
               c0=  2.398778083334100000
            ELSE IF ((delta_P .GT. 50.0) .AND. (delta_P.LE. 1000.0)) THEN
              c6= -0.000000000000000126
              c5=  0.000000000000461219
              c4= -0.000000000686316847
              c3=  0.000000540218636146
              c2= -0.000251691712619862
              c1=  0.083217749423939100
              c0=  5.755370460706180000
            ELSE                           !No correlations for larger diffusion layer thicknesses (YET)
               write(*,*) 'ERROR: delta_P > 1000 micron in Compute_Surface_Solubility'
               STOP
            END IF  
            S_ratio= (c6*delta_P**6) + (c5*delta_P**5) + (c4*delta_P**4) + (c3*delta_P**3) + (c2*delta_P**2) + (c1*delta_P) +(c0)

         ELSE !--- Buffer Capacity= 0.0 mM -----------------------------------------------------------------
            S_ratio =2.30196707
         END IF
         current%pardata%par_conc = S_ratio * S_intrinsic
         !write(*,*) 'iter,ID,delta,Cs',iter,current%pardata%parid,delta_P,current%pardata%par_conc 
      END IF  
   END IF  
   current => next

END DO   

!===================================================================================================
END SUBROUTINE Compute_C_surface
!===================================================================================================








!===================================================================================================
SUBROUTINE Particle_Drug_To_Nodes 
!===================================================================================================
!----- Interpolate Particle concentration release to node locations 
!----- Called by Particle_Track (LBM.f90) to get delphi_particle

IMPLICIT NONE
INTEGER(lng)  :: ID,i,j,k,mpierr
REAL(dbl)     :: xp,yp,zp
REAL(dbl)		  :: delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
REAL(dbl)     :: n_d         				! Modeling parameter to extend the volume of influence around 
REAL(dbl)     :: R_P, Sh_P, delta_P
REAL(dbl)     :: R_influence_P, L_influence_P
REAL(dbl)		  :: tmp
REAL(dbl)   	::  Overlap_test,Overlap_test_Global
TYPE(ParRecord), POINTER  :: current
TYPE(ParRecord), POINTER  :: next

delta_mesh = 1.0_dbl
zcf3 = xcf*ycf*zcf
n_d = 2.5
Overlap_sum_l = 0

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                   !only calculate the drug release when particle radius is larger than 0.1 micron
!-----Calculate length scale for jth particle:  delta = R / Sh
!-----Calculate effective radius: R_influence_P = R + (n_d *delta)
!-----Note: need to convert this into Lattice units and not use the physical length units
!-----Then compute equivalent cubic mesh length scale
      ID= current%pardata%parid 
      R_P  = current%pardata%rp
!     Sh_P = current%pardata%sh
!     delta_P = R_P / Sh_P
      delta_P = R_P
      R_influence_P = (R_P + n_d * delta_P) / xcf
  
!-----iomputing equivalent cubic mesh length scale
      L_influence_P = ( (4.0_dbl*PI/3.0_dbl) * R_influence_P**3.0_dbl)**(1.0_dbl/3.0_dbl)

!-----Global particle location (in whole domain and not in current processor)
      xp= current%pardata%xp 
      yp= current%pardata%yp
      zp= current%pardata%zp

!-----Global Volume of Influence Border (GVIB) for this particle
      GVIB_x(1)= xp - 0.5_dbl * L_influence_P
      GVIB_x(2)= xp + 0.5_dbl * L_influence_P
      GVIB_y(1)= yp - 0.5_dbl * L_influence_P
      GVIB_y(2)= yp + 0.5_dbl * L_influence_P
      GVIB_z(1)= zp - 0.5_dbl * L_influence_P
      GVIB_z(2)= zp + 0.5_dbl * L_influence_P

!-----Global Nodes Effected by Particle 
      GNEP_x(1)= FLOOR(GVIB_x(1))
      GNEP_y(1)= FLOOR(GVIB_y(1))
      GNEP_z(1)= FLOOR(GVIB_z(1))
      GNEP_x(2)= CEILING(GVIB_x(2))
      GNEP_y(2)= CEILING(GVIB_y(2))
      GNEP_z(2)= CEILING(GVIB_z(2))

!-----Taking care of the Z-dir Periodic BC
      GNEP_z_Per(1) = GNEP_z(1)
      GNEP_z_Per(2) = GNEP_z(2)

      IF (GNEP_z(1) .LT. 1) THEN
         GNEP_z_Per(1) = GNEP_z(1) + nz
         GNEP_z_Per(2) = GNEP_z(2) + nz
         GVIB_z_Per(1) = GVIB_z(1) + nz
         GVIB_z_Per(2) = GVIB_z(2) + nz
      ENDIF

      IF (GNEP_z(2) .GT. nz) THEN
          GNEP_z_Per(1) = GNEP_z(1) - nz
          GNEP_z_Per(2) = GNEP_z(2) - nz
          GVIB_z_Per(1) = GVIB_z(1) - nz
          GVIB_z_Per(2) = GVIB_z(2) - nz
      ENDIF

!     Overlap       = 0.0_dbl

!-----Finding processors with overlap with effective volume around the particle       
100   IF ((((GNEP_x(1) .GT. (iMin-1_lng)) .AND. (GNEP_x(1) .LE. iMax)) .OR. ((GNEP_x(2) .GT. (iMin-1_lng)) .AND. (GNEP_x(2) .LE. iMax))) .AND. &
         (((GNEP_y(1) .GT. (jMin-1_lng)) .AND. (GNEP_y(1) .LE. jMax)) .OR. ((GNEP_y(2) .GT. (jMin-1_lng)) .AND. (GNEP_y(2) .LE. jMax))) .AND. &
         (((GNEP_z(1) .GT. (kMin-1_lng)) .AND. (GNEP_z(1) .LE. kMax)) .OR. ((GNEP_z(2) .GT. (kMin-1_lng)) .AND. (GNEP_z(2) .LE. kMax)))  )THEN
    
         NEP_x(1) = Max(GNEP_x(1), iMin) - (iMin-1)          
         NEP_y(1) = Max(GNEP_y(1), jMin) - (jMin-1)
         NEP_z(1) = Max(GNEP_z(1), kMin) - (kMin-1)
         NEP_x(2) = Min(GNEP_x(2), iMax) - (iMin-1)
         NEP_y(2) = Min(GNEP_y(2), jMax) - (jMin-1)
         NEP_z(2) = Min(GNEP_z(2), kMax) - (kMin-1) 

         LVIB_x(1) = GVIB_x(1)- REAL(iMin-1.0_dbl , dbl)
         LVIB_x(2) = GVIB_x(2)- REAL(iMin-1.0_dbl , dbl)
         LVIB_y(1) = GVIB_y(1)- REAL(jMin-1.0_dbl , dbl)
         LVIB_y(2) = GVIB_y(2)- REAL(jMin-1.0_dbl , dbl)
         LVIB_z(1) = GVIB_z(1)- REAL(kMin-1.0_dbl , dbl)
         LVIB_z(2) = GVIB_z(2)- REAL(kMin-1.0_dbl , dbl)

!------- NEW: Finding the volume overlapping between particle-effetive-volume and the volume around each lattice node
         DO i= NEP_x(1),NEP_x(2) 
            DO j= NEP_y(1),NEP_y(2)
               DO k= NEP_z(1),NEP_z(2)
                  NVB_x(1) = REAL(i,dbl) - 0.5_dbl*delta_mesh
                  NVB_x(2) = REAL(i,dbl) + 0.5_dbl*delta_mesh
                  NVB_y(1) = REAL(j,dbl) - 0.5_dbl*delta_mesh
                  NVB_y(2) = REAL(j,dbl) + 0.5_dbl*delta_mesh
                  NVB_z(1) = REAL(k,dbl) - 0.5_dbl*delta_mesh
	                NVB_z(2) = REAL(k,dbl) + 0.5_dbl*delta_mesh
                  IF (node(i,j,k) .EQ. FLUID) THEN
                     Overlap(i,j,k)= MAX ( MIN(LVIB_x(2),NVB_x(2)) - MAX(LVIB_x(1),NVB_x(1)), 0.0_dbl) * & 
                                  MAX ( MIN(LVIB_y(2),NVB_y(2)) - MAX(LVIB_y(1),NVB_y(1)), 0.0_dbl) * &
                                  MAX ( MIN(LVIB_z(2),NVB_z(2)) - MAX(LVIB_z(1),NVB_z(1)), 0.0_dbl)
                     Overlap(i,j,k) = Overlap(i,j,k) * (max((S_bulk-phi(i,j,k) ),0.0_dbl) / S_bulk)
            		     Overlap_sum_l(ID)= Overlap_sum_l(ID) + Overlap(i,j,k)
                  END IF
               END DO
            END DO
         END DO
      END IF
!---- Taking care of the Z-dir Periodic BC
      IF (GNEP_z_Per(1) .ne. GNEP_z(1)) THEN
         GNEP_z(1) = GNEP_z_Per(1)
         GNEP_z(2) = GNEP_z_Per(2)
         GVIB_z(1) = GVIB_z_Per(1)
         GVIB_z(2) = GVIB_z_Per(2) 
         GOTO 100
      ENDIF
    END IF 						! Condition to check if R > R_min_acceptable
    current => next
ENDDO

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
CALL MPI_ALLREDUCE(Overlap_sum_l,Overlap_sum,np+1,MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)


current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN     
      ID= current%pardata%parid 
      IF (Overlap_sum(ID) .LT. 1.0e-9) THEN                        ! Detecting the case of overlap = 0.0
         current%pardata%rp =  (current%pardata%rp**3 + current%pardata%delNBbyCV * (molarvol*zcf3) * (3/(4*PI)) )**(1.0_dbl/3.0_dbl)
         current%pardata%delNBbyCV = 0.0_dbl
      ELSE
         R_P  = current%pardata%rp
         delta_P = R_P
         R_influence_P = (R_P + n_d * delta_P) / xcf
         L_influence_P = ( (4.0_dbl*PI/3.0_dbl) * R_influence_P**3.0_dbl)**(1.0_dbl/3.0_dbl)
         xp= current%pardata%xp 
         yp= current%pardata%yp
         zp= current%pardata%zp
!------- Global Volume of Influence Border (GVIB) for this particle
         GVIB_x(1)= xp - 0.5_dbl * L_influence_P
         GVIB_x(2)= xp + 0.5_dbl * L_influence_P
         GVIB_y(1)= yp - 0.5_dbl * L_influence_P
         GVIB_y(2)= yp + 0.5_dbl * L_influence_P
         GVIB_z(1)= zp - 0.5_dbl * L_influence_P
         GVIB_z(2)= zp + 0.5_dbl * L_influence_P
!------- Global Nodes Effected by Particle 
         GNEP_x(1)= FLOOR(GVIB_x(1))
         GNEP_y(1)= FLOOR(GVIB_y(1))
         GNEP_z(1)= FLOOR(GVIB_z(1))
         GNEP_x(2)= CEILING(GVIB_x(2))
         GNEP_y(2)= CEILING(GVIB_y(2))
         GNEP_z(2)= CEILING(GVIB_z(2))

!------- Taking care of the Z-dir Periodic BC
         GNEP_z_Per(1) = GNEP_z(1)
         GNEP_z_Per(2) = GNEP_z(2)
         IF (GNEP_z(1) .LT. 1) THEN
            GNEP_z_Per(1) = GNEP_z(1) + nz
            GNEP_z_Per(2) = GNEP_z(2) + nz
            GVIB_z_Per(1) = GVIB_z(1) + nz
            GVIB_z_Per(2) = GVIB_z(2) + nz
         ENDIF
         IF (GNEP_z(2) .GT. nz) THEN
            GNEP_z_Per(1) = GNEP_z(1) - nz
            GNEP_z_Per(2) = GNEP_z(2) - nz
            GVIB_z_Per(1) = GVIB_z(1) - nz
            GVIB_z_Per(2) = GVIB_z(2) - nz
         ENDIF
!------- Finding processor that have overlap with effective volume around the particle
200      IF ((((GNEP_x(1) .GT. (iMin-1_lng)) .AND. (GNEP_x(1) .LE. iMax)) .OR. ((GNEP_x(2) .GT. (iMin-1_lng)) .AND. (GNEP_x(2) .LE. iMax)))  .AND. &
             (((GNEP_y(1) .GT. (jMin-1_lng)) .AND. (GNEP_y(1) .LE. jMax)) .OR. ((GNEP_y(2) .GT. (jMin-1_lng)) .AND. (GNEP_y(2) .LE. jMax)))  .AND. &
             (((GNEP_z(1) .GT. (kMin-1_lng)) .AND. (GNEP_z(1) .LE. kMax)) .OR. ((GNEP_z(2) .GT. (kMin-1_lng)) .AND. (GNEP_z(2) .LE. kMax))) ) THEN
            NEP_x(1) = Max(GNEP_x(1), iMin) - (iMin-1)          
            NEP_y(1) = Max(GNEP_y(1), jMin) - (jMin-1)
            NEP_z(1) = Max(GNEP_z(1), kMin) - (kMin-1)
            NEP_x(2) = Min(GNEP_x(2), iMax) - (iMin-1)
            NEP_y(2) = Min(GNEP_y(2), jMax) - (jMin-1)
            NEP_z(2) = Min(GNEP_z(2), kMax) - (kMin-1)

            
            LVIB_x(1) = GVIB_x(1)- REAL(iMin-1.0_dbl , dbl)
            LVIB_x(2) = GVIB_x(2)- REAL(iMin-1.0_dbl , dbl)
            LVIB_y(1) = GVIB_y(1)- REAL(jMin-1.0_dbl , dbl)
            LVIB_y(2) = GVIB_y(2)- REAL(jMin-1.0_dbl , dbl)
            LVIB_z(1) = GVIB_z(1)- REAL(kMin-1.0_dbl , dbl)
            LVIB_z(2) = GVIB_z(2)- REAL(kMin-1.0_dbl , dbl)
!---------- Computing particle release contribution to scalar field at each lattice node
            DO i= NEP_x(1),NEP_x(2)
               DO j= NEP_y(1),NEP_y(2)
                  DO k= NEP_z(1),NEP_z(2)
                     NVB_x(1) = REAL(i,dbl) - 0.5_dbl*delta_mesh
                     NVB_x(2) = REAL(i,dbl) + 0.5_dbl*delta_mesh
                     NVB_y(1) = REAL(j,dbl) - 0.5_dbl*delta_mesh
                     NVB_y(2) = REAL(j,dbl) + 0.5_dbl*delta_mesh
                     NVB_z(1) = REAL(k,dbl) - 0.5_dbl*delta_mesh
	                   NVB_z(2) = REAL(k,dbl) + 0.5_dbl*delta_mesh

                     IF (node(i,j,k) .EQ. FLUID) THEN                 
                        Overlap(i,j,k)= MAX ( MIN(LVIB_x(2),NVB_x(2)) - MAX(LVIB_x(1),NVB_x(1)), 0.0_dbl) * & 
                                        MAX ( MIN(LVIB_y(2),NVB_y(2)) - MAX(LVIB_y(1),NVB_y(1)), 0.0_dbl) * &
                                        MAX ( MIN(LVIB_z(2),NVB_z(2)) - MAX(LVIB_z(1),NVB_z(1)), 0.0_dbl)
                        Overlap(i,j,k) = Overlap(i,j,k) * (max((S_bulk-phi(i,j,k) ),0.0_dbl) / S_bulk)
                        Overlap(i,j,k) = Overlap(i,j,k) / Overlap_sum(ID)
          	            delphi_particle(i,j,k)  = delphi_particle(i,j,k)  + current%pardata%delNBbyCV * Overlap(i,j,k) 
                     END IF 
                  END DO
               END DO
            END DO
         END IF

         IF (GNEP_z_Per(1) .ne. GNEP_z(1)) THEN
            GNEP_z(1) = GNEP_z_Per(1)
            GNEP_z(2) = GNEP_z_Per(2)
            GVIB_z(1) = GVIB_z_Per(1)
            GVIB_z(2) = GVIB_z_Per(2)
            GOTO 200
         ENDIF
      ENDIF            ! Condition to check if Overlap_Sum > 1e-9
   END IF             ! Condition to check if R > R_min_acceptable
   current => next
ENDDO
!===================================================================================================
END SUBROUTINE Particle_Drug_To_Nodes  		 
!===================================================================================================





!================================================
END MODULE ParticleDrug 
!================================================
