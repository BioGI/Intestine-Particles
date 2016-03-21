!==================================================================================================
MODULE ParticleDrug				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================
USE SetPrecision
USE Setup
USE ICBC
USE MPI

IMPLICIT NONE

CONTAINS

!!===================================================================================================
!SUBROUTINE Interp_bulkconc(Cb_Local)                                 ! Using Trilinear interpolation
!!===================================================================================================
!IMPLICIT NONE
!INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
!REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd
!REAL(dbl)     :: xp,yp,zp
!REAL(dbl)     :: Cb_Local
!TYPE(ParRecord), POINTER :: current
!TYPE(ParRecord), POINTER :: next
!
!current => ParListHead%next
!DO WHILE (ASSOCIATED(current))
!   next => current%next ! copy pointer of next node
!   IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      xp = current%pardata%xp - REAL(iMin-1_lng,dbl)
!      yp = current%pardata%yp - REAL(jMin-1_lng,dbl)
!      zp = current%pardata%zp - REAL(kMin-1_lng,dbl)
!
!      ix0= FLOOR(xp)
!      ix1= CEILING(xp)
!      iy0= FLOOR(yp)
!      iy1= CEILING(yp)
!      iz0= FLOOR(zp)
!      iz1= CEILING(zp)
!!!!!! MAKE SURE THE ABOVE NODES ARE FLUID NODES
!
!      IF (ix1 /= ix0) THEN 
!         xd= (xp-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
!      ELSE
!         xd= 0.0_dbl
!      END IF
!      IF (iy1 /= iy0) THEN 
!         yd=(yp-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
!      ELSE
!         yd = 0.0_dbl
!      END IF
!      IF (iz1 /= iz0) THEN 
!         zd=(zp-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
!      ELSE
!         zd = 0.0_dbl
!      END IF
!
!!-----phi-interpolation
!!-----1st level linear interpolation in x-direction
!      c00 = phi(ix0,iy0,iz0)*(1.0_dbl-xd)+phi(ix1,iy0,iz0)*xd	
!      c01 = phi(ix0,iy0,iz1)*(1.0_dbl-xd)+phi(ix1,iy0,iz1)*xd	
!      c10 = phi(ix0,iy1,iz0)*(1.0_dbl-xd)+phi(ix1,iy1,iz0)*xd	
!      c11 = phi(ix0,iy1,iz1)*(1.0_dbl-xd)+phi(ix1,iy1,iz1)*xd	
!
!!-----2nd level linear interpolation in y-direction
!      c0  = c00*(1.0_dbl-yd)+c10*yd
!      c1  = c01*(1.0_dbl-yd)+c11*yd
!
!!-----3rd level linear interpolation in z-direction
!	c   = c0*(1.0_dbl-zd)+c1*zd
!        Cb_Local= c        
!
!   END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   current => next
!ENDDO
!!===================================================================================================
!END SUBROUTINE Interp_bulkconc  
!!===================================================================================================
!
!
!
!
!!===================================================================================================
!SUBROUTINE Calc_Global_Bulk_Scalar_Conc(Cb_Domain)         !Bulk Conc= total moles/total domain size
!!===================================================================================================
!IMPLICIT NONE
!INTEGER(lng)  :: i,j,k
!REAL(dbl)     :: Cb_Domain
!
!Cb_global = 0.0_dbl
!Cb_numFluids = 0_lng
!
!DO k=1,nzSub
!   DO j=1,nySub
!      DO i=1,nxSub
!         IF (node(i,j,k) .EQ. FLUID) THEN
!            Cb_global = Cb_global + phi(i,j,k)
!            Cb_numFluids = Cb_numFluids + 1_lng
!         END IF
!      END DO
!   END DO
!END DO
!
!Cb_Domain = Cb_global/ Cb_numFluids
!!===================================================================================================
!END SUBROUTINE Calc_Global_Bulk_Scalar_Conc
!!===================================================================================================







!===================================================================================================
SUBROUTINE Compute_Cb				  ! Computes the mesh-independent bulk concentration
!===================================================================================================
IMPLICIT NONE

INTEGER(lng)  		 :: i,j,k, mpierr
INTEGER(lng)  		 :: ix0,ix1,iy0,iy1,iz0,iz00,iz1,iz11		! Trilinear interpolation parameters
INTEGER(dbl)		 :: NumFluids_Veff_l, NumFluids_Veff
INTEGER,DIMENSION(2)   	 :: LN_x,  LN_y,  LN_z				! Lattice Nodes Surronding the particle
INTEGER,DIMENSION(2)     :: GNEP_x, GNEP_y, GNEP_z                      ! Lattice Nodes Surronding the particle (Global: not considering the partitioning for parallel processing)
INTEGER,DIMENSION(2)     :: NEP_x,   NEP_y,  NEP_z                      ! Lattice Nodes Surronding the particle (Local: in current processor)
REAL(dbl)     		 :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd		! Trilinear interpolation parameters
REAL(dbl)  	   	 :: xp,yp,zp
REAL(dbl)		 :: delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
REAL(dbl)       	 :: N_b         				! Modeling parameter to extend the volume of influence  
REAL(dbl)    	         :: R_P, Sh_P, delta_P
REAl(dbl)                :: R_influence_p, L_influence_p		! Parameters related to particle's volume of influence
REAl(dbl)                :: V_influence_P	 			! Parameters related to particle's volume of influence
REAL(dbl)		 :: Cb_Total_Veff_l, Cb_Total_Veff
REAL(dbl),DIMENSION(2)   :: VIB_x, VIB_y, VIB_z	, VIB_z_Per 			! Volume of Influence's Borders
REAL(dbl),DIMENSION(2)   :: NVB_x, NVB_y, NVB_z				! Node Volume's Borders
REAL(dbl)                :: Delta_L
REAL(dbl)                :: x_DP, y_DP, z_DP				! Coordinates of "Discretized Point" (DP)
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

delta_mesh = 1.0_dbl
zcf3 = 1.0_dbl

current => ParListHead%next
DO WHILE (ASSOCIATED(current))

!------ Copy pointer of next node
	next => current%next

!------ Particle length scale: delta= R/Sh & effective radius: R_influence_P= R+(N_b*delta)
	N_b = 1.0
        R_P = current%pardata%rp
	Sh_P= current%pardata%sh
        delta_P= R_P/Sh_P
        R_influence_P= (R_P+N_b*delta_P)/xcf

!------ Computing equivalent cubic mesh length scale
        V_influence_P= (4.0_dbl/3.0_dbl)*PI* R_influence_P**3.0_dbl
        L_influence_P= V_influence_P **(1.0_dbl/3.0_dbl)
        V_eff_Ratio  = V_influence_P/zcf3 					! Ratio of the effective volume to cell size 

        Cb_Total_Veff_l  = 0.0_lng
        Cb_Total_Veff    = 0.0_lng
        NumFluids_Veff_l = 0.0_lng
        NumFluids_Veff   = 0.0_lng

!------ Veff is smaller than the mesh volume --> Cb = Trilinear interpolation of the concentration at particle location
        IF (V_eff_Ratio .LE. 1.0) THEN 					
           IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              CaseNo = 1
!------------ Finding local particle location (at current processor)
              xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
              yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
              zp= current%pardata%zp - REAL(kMin-1_lng,dbl)
              ix0 =FLOOR(xp)
              ix1 =CEILING(xp)
              iy0 =FLOOR(yp)
              iy1 =CEILING(yp)
              iz0 =FLOOR(zp)
              iz1 =CEILING(zp)
!------------ TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
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
!------------ Concentration Trilinear Iinterpolation
!------------ Interpolation in x-direction
              c00 = phi(ix0,iy0,iz0) * (1.0_dbl-xd) + phi(ix1,iy0,iz0) * xd
              c01 = phi(ix0,iy0,iz1) * (1.0_dbl-xd) + phi(ix1,iy0,iz1) * xd
              c10 = phi(ix0,iy1,iz0) * (1.0_dbl-xd) + phi(ix1,iy1,iz0) * xd
              c11 = phi(ix0,iy1,iz1) * (1.0_dbl-xd) + phi(ix1,iy1,iz1) * xd
!------------ Interpolation in y-direction
              c0  = c00 * (1.0_dbl-yd) + c10 * yd
              c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------------ Interpolation in z-direction
              c   = c0 * (1.0_dbl-zd) + c1 * zd
              Cb_Hybrid= c 
           END IF !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
!-----------------------------------------------------------------------------------------------------------------------------
        ELSE  !------------------------------------Veff is larger than 1 which means parallel communication might be necessary 
!-----------------------------------------------------------------------------------------------------------------------------   

!--------- Finding particle global location (in the whole domain not this particular processor)
           xp= current%pardata%xp 
           yp= current%pardata%yp 
           zp= current%pardata%zp 

!--------- Global Volume of Influence Border (VIB) for this particle (in the whole domain not this particular processor)
           VIB_x(1)= xp - 0.5_dbl * L_influence_P
           VIB_x(2)= xp + 0.5_dbl * L_influence_P
           VIB_y(1)= yp - 0.5_dbl * L_influence_P
           VIB_y(2)= yp + 0.5_dbl * L_influence_P
           VIB_z(1)= zp - 0.5_dbl * L_influence_P
           VIB_z(2)= zp + 0.5_dbl * L_influence_P

!--------- Finding processor that have overlap with effective volume around the particle       
           IF( (((VIB_x(1) .GE. (iMin-1_lng)) .AND. (VIB_x(1) .LT. iMax)) .OR. ((VIB_x(2) .GE. (iMin-1_lng)) .AND. (VIB_x(2) .LT. iMax))) .AND. &
               (((VIB_y(1) .GE. (jMin-1_lng)) .AND. (VIB_y(1) .LT. jMax)) .OR. ((VIB_y(2) .GE. (jMin-1_lng)) .AND. (VIB_y(2) .LT. jMax))) .AND. &
               (((VIB_z(1) .GE. (kMin-1_lng)) .AND. (VIB_z(1) .LT. kMax)) .OR. ((VIB_z(2) .GE. (kMin-1_lng)) .AND. (VIB_z(2) .LT. kMax)))  )THEN

!-------------------------------------------------------------------------------------------------------------------------
!------------- Veff is slightly larger than mesh volume --> Volume of influence is discretized
!------------- Cb= Average of concentration interpolated on each of the descritized nodes inside volume of influence
!-------------------------------------------------------------------------------------------------------------------------
               IF ( (V_eff_Ratio .GT. 1.0) .AND. (V_eff_Ratio .LT. 27.0 ) ) THEN		
                  CaseNo = 2

!---------------- Discretizing the volume of influence to  make sure at least 27 points are available
                  Delta_L = (VIB_x(2)-VIB_x(1)) / 2.0 

!---------------- Loop over discretized points and averaging the concentration
                  DO i= 0, 2
                     DO j= 0, 2
                        DO k= 0, 2
                           x_DP = VIB_x(1) + (i * Delta_L) 
                           y_DP = VIB_y(1) + (j * Delta_L)
                           z_DP = VIB_z(1) + (k * Delta_L)
                           IF( (x_DP .GE. (REAL(iMin,dbl)-1.0_dbl)) .AND. &               
			       (x_DP .LT.  REAL(iMax,dbl)         ) .AND. &
                               (y_DP .GE. (REAL(jMin,dbl)-1.0_dbl)) .AND. & 
			       (y_DP .LT.  REAL(jMax,dbl)         ) .AND. &
                               (z_DP .GE. (REAL(kMin,dbl)-1.0_dbl)) .AND. & 
			       (z_DP .LT.  REAL(kMax,dbl)         ) ) THEN

!----------------------------- Finding Local lattice nodes surrounding this point (This point is discretized and is not a lattice node))
                               ix0 = FLOOR(x_DP)   - (REAL(iMin,dbl)-1.0_dbl)
                               ix1 = CEILING(x_DP) - (REAL(iMin,dbl)-1.0_dbl)
                               iy0 = FLOOR(y_DP)   - (REAL(jMin,dbl)-1.0_dbl)
                               iy1 = CEILING(y_DP) - (REAL(jMin,dbl)-1.0_dbl) 
                               iz0 = FLOOR(z_DP)   - (REAL(kMin,dbl)-1.0_dbl)
                               iz1 = CEILING(z_DP) - (REAL(kMin,dbl)-1.0_dbl)
                              
            		       x_DP = x_DP - REAL(iMin-1_lng,dbl)
                               y_DP = y_DP - REAL(jMin-1_lng,dbl)
                               z_DP = z_DP - REAL(kMin-1_lng,dbl)
 
!----------------------------- TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
                               IF (ix1 /= ix0) THEN
                                  xd=(x_DP-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
                               ELSE
                                  xd = 0.0_dbl
                               END IF
 
                               IF (iy1 /= iy0) THEN
                                  yd=(y_DP-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
                               ELSE
                                  yd = 0.0_dbl
                               END IF
        
                               IF (iz1 /= iz0) THEN
                                  zd=(z_DP-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
                               ELSE
                                  zd = 0.0_dbl
                               END IF

!----------------------------- Concentration Trilinear Iinterpolation
!----------------------------- Interpolation in x-direction
                               c00 = phi(ix0,iy0,iz0) * (1.0_dbl-xd) + phi(ix1,iy0,iz0) * xd
                               c01 = phi(ix0,iy0,iz1) * (1.0_dbl-xd) + phi(ix1,iy0,iz1) * xd
                               c10 = phi(ix0,iy1,iz0) * (1.0_dbl-xd) + phi(ix1,iy1,iz0) * xd
                               c11 = phi(ix0,iy1,iz1) * (1.0_dbl-xd) + phi(ix1,iy1,iz1) * xd
!----------------------------- Interpolation in y-direction
                               c0  = c00 * (1.0_dbl-yd) + c10 * yd
                               c1  = c01 * (1.0_dbl-yd) + c11 * yd
!----------------------------- Interpolation in z-direction
                               c   = c0 * (1.0_dbl-zd) + c1 * zd
 
                               Cb_Total_Veff_l  = Cb_Total_Veff_l  + c
                               NumFluids_Veff_l = NumFluids_Veff_l + 1_lng
                           END IF
                        END DO
                    END DO
                 END DO
          
!----------------------------------------------------------------------------------------------------------------------
!------------- Veff is much larger than mesh volume --> Cb= total number of moles in volume of influence / volume of influence 
!----------------------------------------------------------------------------------------------------------------------
               ELSE IF (V_eff_Ratio .GE. 27.0) THEN                             
                  CaseNo = 3


!---------------- Finding the lattice "Nodes Effected by Particle"
                  GNEP_x(1)= CEILING(VIB_x(1))
                  GNEP_x(2)= FLOOR  (VIB_x(2))
                  GNEP_y(1)= CEILING(VIB_y(1))
                  GNEP_y(2)= FLOOR  (VIB_y(2))
                  GNEP_z(1)= CEILING(VIB_z(1))
                  GNEP_z(2)= FLOOR  (VIB_z(2))

!---------------- Finding the lattice "Nodes Effected by Particle"
                  NEP_x(1) = Max(GNEP_x(1) , iMin)
                  NEP_y(1) = Max(GNEP_y(1) , jMin)
                  NEP_z(1) = Max(GNEP_z(1) , kMin)

                  NEP_x(2) = Min(GNEP_x(2) , iMax)
                  NEP_y(2) = Min(GNEP_y(2) , jMax)
                  NEP_z(2) = Min(GNEP_z(2) , kMax)

                  NEP_x(1) = NEP_x(1) - (iMin-1)
                  NEP_x(2) = NEP_x(2) - (iMin-1)
                  NEP_y(1) = NEP_y(1) - (jMin-1)
                  NEP_y(2) = NEP_y(2) - (jMin-1)
                  NEP_z(1) = NEP_z(1) - (kMin-1)
                  NEP_z(2) = NEP_z(2) - (kMin-1)

                  !write(*,*) iter,mySub,'A:NEP',NEP_x(1),NEP_x(2),NEP_y(1),NEP_y(2),NEP_z(1),NEP_z(2)

                  DO i= NEP_x(1),NEP_x(2) 
                     DO j= NEP_y(1),NEP_y(2)
                        DO k= NEP_z(1),NEP_z(2)
                           IF (node(i,j,k) .EQ. FLUID) THEN
                              Cb_Total_Veff_l  = Cb_Total_Veff_l  + phi(i,j,k)
                              NumFluids_Veff_l = NumFluids_Veff_l + 1_lng
                           END IF
                        END DO
                     END DO
                  END DO

                  !write(*,*) iter,mySub, ' A-Cb1',Cb_Total_Veff_l, NumFluids_Veff_l

               END IF  									! Conditional for cases 2 and 3
         END IF 									! Conditional for the processor which has overlap with effective volume 




!--------------------------------------------------------------
!------- TAKING CARE OF THE PERIODIC BC
!--------------------------------------------------------------
         VIB_z_Per(1) = VIB_z(1) 
         VIB_z_Per(2) = VIB_z(2)  

         IF (VIB_z(1) .LT. 1) THEN
            VIB_z_Per(1) = VIB_z(1) + nz
            VIB_z_Per(2) = VIB_z(2) + nz
         ENDIF

         IF (VIB_z(2) .GT. nz) THEN
            VIB_z_Per(1) = VIB_z(1) - nz
            VIB_z_Per(2) = VIB_z(2) - nz
         ENDIF

        IF (VIB_z_Per(1) .NE. VIB_z(1)) THEN
!--------- Finding processor that have overlap with effective volume around the particle       
           IF( (((VIB_x(1) .GE. (iMin-1_lng)) .AND. (VIB_x(1) .LT. iMax)) .OR. ((VIB_x(2) .GE. (iMin-1_lng)) .AND. (VIB_x(2) .LT. iMax))) .AND. &
               (((VIB_y(1) .GE. (jMin-1_lng)) .AND. (VIB_y(1) .LT. jMax)) .OR. ((VIB_y(2) .GE. (jMin-1_lng)) .AND. (VIB_y(2) .LT. jMax))) .AND. &
               (((VIB_z_Per(1) .GE. (kMin-1_lng)) .AND. (VIB_z_Per(1) .LT. kMax)) .OR. ((VIB_z_Per(2) .GE. (kMin-1_lng)) .AND. (VIB_z_Per(2) .LT. kMax)))  )THEN

!-------------------------------------------------------------------------------------------------------------------------
!------------- Veff is slightly larger than mesh volume --> Volume of influence is discretized
!------------- Cb= Average of concentration interpolated on each of the descritized nodes inside volume of influence
!-------------------------------------------------------------------------------------------------------------------------
               IF ( (V_eff_Ratio .GT. 1.0) .AND. (V_eff_Ratio .LT. 27.0 ) ) THEN		
                  CaseNo = 2

!---------------- Discretizing the volume of influence to  make sure at least 27 points are available
                  Delta_L = (VIB_x(2)-VIB_x(1)) / 2.0 

!---------------- Loop over discretized points and averaging the concentration
                  DO i= 0, 2
                     DO j= 0, 2
                        DO k= 0, 2
                           x_DP = VIB_x(1) + (i * Delta_L) 
                           y_DP = VIB_y(1) + (j * Delta_L)
                           z_DP = VIB_z_Per(1) + (k * Delta_L)
                           IF( (x_DP .GE. (REAL(iMin,dbl)-1.0_dbl)) .AND. &               
			       (x_DP .LT.  REAL(iMax,dbl)         ) .AND. &
                               (y_DP .GE. (REAL(jMin,dbl)-1.0_dbl)) .AND. & 
			       (y_DP .LT.  REAL(jMax,dbl)         ) .AND. &
                               (z_DP .GE. (REAL(kMin,dbl)-1.0_dbl)) .AND. & 
			       (z_DP .LT.  REAL(kMax,dbl)         ) ) THEN

!----------------------------- Finding Local lattice nodes surrounding this point (This point is discretized and is not a lattice node))
                               ix0 = FLOOR(x_DP)   - (REAL(iMin,dbl)-1.0_dbl)
                               ix1 = CEILING(x_DP) - (REAL(iMin,dbl)-1.0_dbl)
                               iy0 = FLOOR(y_DP)   - (REAL(jMin,dbl)-1.0_dbl)
                               iy1 = CEILING(y_DP) - (REAL(jMin,dbl)-1.0_dbl) 
                               iz0 = FLOOR(z_DP)   - (REAL(kMin,dbl)-1.0_dbl)
                               iz1 = CEILING(z_DP) - (REAL(kMin,dbl)-1.0_dbl)
                              
            		       x_DP = x_DP - REAL(iMin-1_lng,dbl)
                               y_DP = y_DP - REAL(jMin-1_lng,dbl)
                               z_DP = z_DP - REAL(kMin-1_lng,dbl)
 
!----------------------------- TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
                               IF (ix1 /= ix0) THEN
                                  xd=(x_DP-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
                               ELSE
                                  xd = 0.0_dbl
                               END IF
 
                               IF (iy1 /= iy0) THEN
                                  yd=(y_DP-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
                               ELSE
                                  yd = 0.0_dbl
                               END IF
        
                               IF (iz1 /= iz0) THEN
                                  zd=(z_DP-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
                               ELSE
                                  zd = 0.0_dbl
                               END IF

!----------------------------- Concentration Trilinear Iinterpolation
!----------------------------- Interpolation in x-direction
                               c00 = phi(ix0,iy0,iz0) * (1.0_dbl-xd) + phi(ix1,iy0,iz0) * xd
                               c01 = phi(ix0,iy0,iz1) * (1.0_dbl-xd) + phi(ix1,iy0,iz1) * xd
                               c10 = phi(ix0,iy1,iz0) * (1.0_dbl-xd) + phi(ix1,iy1,iz0) * xd
                               c11 = phi(ix0,iy1,iz1) * (1.0_dbl-xd) + phi(ix1,iy1,iz1) * xd
!----------------------------- Interpolation in y-direction
                               c0  = c00 * (1.0_dbl-yd) + c10 * yd
                               c1  = c01 * (1.0_dbl-yd) + c11 * yd
!----------------------------- Interpolation in z-direction
                               c   = c0 * (1.0_dbl-zd) + c1 * zd
 
                               Cb_Total_Veff_l  = Cb_Total_Veff_l  + c
                               NumFluids_Veff_l = NumFluids_Veff_l + 1_lng
                           END IF
                        END DO
                    END DO
                 END DO
          
!----------------------------------------------------------------------------------------------------------------------
!------------- Veff is much larger than mesh volume --> Cb= total number of moles in volume of influence / volume of influence 
!----------------------------------------------------------------------------------------------------------------------
               ELSE IF (V_eff_Ratio .GE. 27.0) THEN                             
                  CaseNo = 3

!---------------- Finding the lattice "Nodes Effected by Particle"
                  GNEP_x(1)= CEILING(VIB_x(1))
                  GNEP_x(2)= FLOOR  (VIB_x(2))
                  GNEP_y(1)= CEILING(VIB_y(1))
                  GNEP_y(2)= FLOOR  (VIB_y(2))
                  GNEP_z(1)= CEILING(VIB_z_Per(1))
                  GNEP_z(2)= FLOOR  (VIB_z_Per(2))

!---------------- Finding the lattice "Nodes Effected by Particle"
                  NEP_x(1) = Max(GNEP_x(1) , iMin)
                  NEP_y(1) = Max(GNEP_y(1) , jMin)
                  NEP_z(1) = Max(GNEP_z(1) , kMin)

                  NEP_x(2) = Min(GNEP_x(2) , iMax)
                  NEP_y(2) = Min(GNEP_y(2) , jMax)
                  NEP_z(2) = Min(GNEP_z(2) , kMax)

                  NEP_x(1) = NEP_x(1) - (iMin-1)
                  NEP_x(2) = NEP_x(2) - (iMin-1)
                  NEP_y(1) = NEP_y(1) - (jMin-1)
                  NEP_y(2) = NEP_y(2) - (jMin-1)
                  NEP_z(1) = NEP_z(1) - (kMin-1)
                  NEP_z(2) = NEP_z(2) - (kMin-1)

                  !write(*,*) iter,mySub,'B:NEP',NEP_x(1),NEP_x(2),NEP_y(1),NEP_y(2),NEP_z(1),NEP_z(2)

                  DO i= NEP_x(1),NEP_x(2) 
                     DO j= NEP_y(1),NEP_y(2)
                        DO k= NEP_z(1),NEP_z(2)
                           IF (node(i,j,k) .EQ. FLUID) THEN
                              Cb_Total_Veff_l  = Cb_Total_Veff_l  + phi(i,j,k)
                              NumFluids_Veff_l = NumFluids_Veff_l + 1_lng
                           END IF
                        END DO
                     END DO
                  END DO

                  !write(*,*) iter,mySub, ' B-Cb1',Cb_Total_Veff_l, NumFluids_Veff_l

               END IF  									! Conditional for cases 2 and 3
         END IF 									! Conditional for the processor which has overlap with effective volume 

       ENDIF 

         CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
         CALL MPI_ALLREDUCE(Cb_Total_Veff_l , Cb_Total_Veff , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
         CALL MPI_ALLREDUCE(NumFluids_Veff_l, NumFluids_Veff, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, mpierr)

         Cb_Hybrid= Cb_Total_Veff / NumFluids_Veff
         current%pardata%bulk_conc = Cb_Hybrid
	
      END IF       			                                    		!Conditional for V_eff 
      
      open(172,file='Cb-'//sub//'.dat', position='append')

      current => next
END DO
!===================================================================================================
END SUBROUTINE Compute_Cb
!===================================================================================================







!===================================================================================================
SUBROUTINE Particle_Drug_Release                     ! Calculates drug release at every time step  
!===================================================================================================
IMPLICIT NONE
INTEGER(lng)  :: numFluids,i,j,k,RANK,mpierr
REAL(dbl)     :: deltaR,temp,cbt,zcf3,bulkconc 
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

zcf3=xcf*ycf*zcf

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (mySub .EQ.current%pardata%cur_part) THEN !+++++++++++++++++++++++++++++++++++++++++++++++++++
      current%pardata%rpold = current%pardata%rp
      bulkconc = current%pardata%bulk_conc
      temp = current%pardata%rpold**2.0_dbl-4.0_dbl*tcf*molarvol*diffm*current%pardata%sh*max((current%pardata%par_conc-bulkconc),0.0_dbl)
 
      IF (temp.GE.0.0_dbl) THEN
          current%pardata%rp=0.5_dbl*(current%pardata%rpold+sqrt(temp))
      ELSE
          temp = 0.0_dbl
          current%pardata%rp=0.5_dbl*(current%pardata%rpold+sqrt(temp))
      END IF

      deltaR=current%pardata%rpold-current%pardata%rp
      current%pardata%delNBbyCV = (4.0_dbl/3.0_dbl) * PI*(current%pardata%rpold**3.0_dbl - current%pardata%rp**3.0_dbl) /(molarvol*zcf3)
   END IF !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   RANK= current%pardata%cur_part - 1
   CALL MPI_BCast(current%pardata%delNBbyCV, 1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD, mpierr)
   CALL MPI_BCast(current%pardata%rp,        1, MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD, mpierr)

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
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
   IF (mySub .EQ.current%pardata%cur_part) THEN !+++++++++++++++++++++++++++++++++++++++++++++++++++
      current%pardata%sh= 1.0_dbl + (current%pardata%gamma_cont / (1.0_dbl-current%pardata%gamma_cont)) 	! Add container effect
      S= current%pardata%S
      Sst= S* (current%pardata%rp**2.0) / diffm
      current%pardata%Sst= Sst

      IF (Sst.LT.5.0_dbl) THEN
         current%pardata%sh = current%pardata%sh + 0.296_dbl*(Sst**0.5_dbl)
      ELSE
         Sh0 = exp(0.162_dbl + 0.202_dbl*log(Sst) - 7.5e-6_dbl*(log(Sst)**5.4_dbl)) 
         current%pardata%sh = current%pardata%sh + Sh0-1.0_dbl
      END IF
   END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
   RANK= current%pardata%cur_part - 1
   CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   CALL MPI_BCast(current%pardata%sh,1,MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)

   current => next
ENDDO
!===================================================================================================
END SUBROUTINE Compute_Sherwood
!===================================================================================================






!===================================================================================================
SUBROUTINE Compute_shear                
!===================================================================================================
! Computing components of strain rate tensor using central difference everywhere except near the
! processor boundaries where sided difference is used. 

IMPLICIT NONE
INTEGER(lng)  :: i,j,k
INTEGER(lng)  :: it,jt,kt,ib,jb,kb
REAL(dbl)     :: xp,yp,zp,xd,yd,zd
INTEGER(lng)  :: ix0,iy0,iz0,ix1,iy1,iz1
REAL(dbl)     :: temp,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
REAL(dbl)     :: S
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next
   IF (mySub .EQ.current%pardata%cur_part) THEN !+++++++++++++++++++++++++++++++++++++++++++++++++++
      xp = current%pardata%xp - REAL(iMin-1_lng,dbl)
      yp = current%pardata%yp - REAL(jMin-1_lng,dbl)
      zp = current%pardata%zp - REAL(kMin-1_lng,dbl)

      ix0= FLOOR(xp)
      ix1= CEILING(xp)
      iy0= FLOOR(yp)
      iy1= CEILING(yp)
      iz0= FLOOR(zp)
      iz1= CEILING(zp)

      ib = ix0
      jb = iy0
      kb = iz0
      it = ix0 + 1_lng
      jt = iy0
      kt = iz0
!!!!! MAKE SURE THE ABOVE NODES ARE FLUID NODES

      dwdz = w(it,jt,kt) - w(ib,jb,kb)
      S = abs(dwdz*vcf/zcf)
      current%pardata%S = S
    END IF !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    current => next
ENDDO

!===================================================================================================
END SUBROUTINE Compute_shear
!===================================================================================================







!===================================================================================================
SUBROUTINE Particle_Drug_To_Nodes 
!===================================================================================================
!----- Interpolate Particle concentration release to node locations 
!----- Called by Particle_Track (LBM.f90) to get delphi_particle

IMPLICIT NONE
INTEGER(lng)  		  :: i,j,k,mpierr
REAL(dbl)     		  :: xp,yp,zp
REAL(dbl)		  :: delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
REAL(dbl)                 :: N_d         				! Modeling parameter to extend the volume of influence around 
REAL(dbl)                 :: R_P, Sh_P, delta_P
REAL(dbl)                 :: R_influence_P, L_influence_P
REAL(dbl),DIMENSION(2)    :: VIB_x, VIB_y, VIB_z, VIB_z_Per 		! Volume of Influence's Borders
REAL(dbl),DIMENSION(2)    :: NVB_x, NVB_y, NVB_z			! Node Volume's Borders
INTEGER  ,DIMENSION(2)    :: LN_x,  LN_y,  LN_z				! Lattice Nodes Surronding the particle
INTEGER  ,DIMENSION(2)    :: GNEP_x, GNEP_y, GNEP_z, GNEP_z_Per         ! Lattice Nodes Surronding the particle (Global: not considering the partitioning for parallel processing)
INTEGER  ,DIMENSION(2)    :: NEP_x,   NEP_y,  NEP_z                     ! Lattice Nodes Surronding the particle (Local: in current processor)
REAL(dbl)		  :: tmp, Overlap_sum_l, Overlap_sum
TYPE(ParRecord), POINTER  :: current
TYPE(ParRecord), POINTER  :: next

delta_mesh = 1.0_dbl
zcf3 = xcf*ycf*zcf
current => ParListHead%next

DO WHILE (ASSOCIATED(current))
   next => current%next 

!--Calculate length scale for jth particle:  delta = R / Sh
!--Calculate effective radius: R_influence_P = R + (N_d *delta)
!--Note: need to convert this into Lattice units and not use the physical length units
!--Then compute equivalent cubic mesh length scale
   N_d = 1.0
   R_P  = current%pardata%rp
   Sh_P = current%pardata%sh
   delta_P = R_P / Sh_P
   R_influence_P = (R_P + N_d * delta_P) / xcf

!--iomputing equivalent cubic mesh length scale
   L_influence_P = ( (4.0_dbl*PI/3.0_dbl) * R_influence_P**3.0_dbl)**(1.0_dbl/3.0_dbl)

!--Global particle location (in whole domain and not in current processor)
   xp= current%pardata%xp 
   yp= current%pardata%yp
   zp= current%pardata%zp

!--Global Volume of Influence Border (VIB) for this particle
   VIB_x(1)= xp - 0.5_dbl * L_influence_P
   VIB_x(2)= xp + 0.5_dbl * L_influence_P
   VIB_y(1)= yp - 0.5_dbl * L_influence_P
   VIB_y(2)= yp + 0.5_dbl * L_influence_P
   VIB_z(1)= zp - 0.5_dbl * L_influence_P
   VIB_z(2)= zp + 0.5_dbl * L_influence_P

!--Global Nodes Effected by Particle 
   GNEP_x(1)= FLOOR(VIB_x(1))
   GNEP_x(2)= CEILING(VIB_x(2))
   GNEP_y(1)= FLOOR(VIB_y(1))
   GNEP_y(2)= CEILING(VIB_y(2))
   GNEP_z(1)= FLOOR(VIB_z(1))
   GNEP_z(2)= CEILING(VIB_z(2))

!--Taking care of the Z-dir Periodic BC
   GNEP_z_Per(1) = GNEP_z(1)
   GNEP_z_Per(2) = GNEP_z(2)

   IF (GNEP_z(1) .LT. 1) THEN
       GNEP_z_Per(1) = GNEP_z(1) + nz
       GNEP_z_Per(2) = GNEP_z(2) + nz
       VIB_z_Per(1)  = VIB_z(1)  + nz
       VIB_z_Per(2)  = VIB_z(2)  + nz
   ENDIF

   IF (GNEP_z(2) .GT. nz) THEN
       GNEP_z_Per(1) = GNEP_z(1) - nz
       GNEP_z_Per(2) = GNEP_z(2) - nz
       VIB_z_Per(1)  = VIB_z(1)  - nz
       VIB_z_Per(2)  = VIB_z(2)  - nz
   ENDIF

   Overlap_sum_l = 0.0_dbl
   Overlap= 0.0

!--Finding processors with overlap with effective volume around the particle       
100 IF((((GNEP_x(1) .GT. (iMin-1_lng)) .AND. (GNEP_x(1) .LE. iMax)) .OR. ((GNEP_x(2) .GT. (iMin-1_lng)) .AND. (GNEP_x(2) .LE. iMax))) .AND. &
      (((GNEP_y(1) .GT. (jMin-1_lng)) .AND. (GNEP_y(1) .LE. jMax)) .OR. ((GNEP_y(2) .GT. (jMin-1_lng)) .AND. (GNEP_y(2) .LE. jMax))) .AND. &
      (((GNEP_z(1) .GT. (kMin-1_lng)) .AND. (GNEP_z(1) .LE. kMax)) .OR. ((GNEP_z(2) .GT. (kMin-1_lng)) .AND. (GNEP_z(2) .LE. kMax)))  )THEN
    
      NEP_x(1) = Max(GNEP_x(1), iMin)            
      NEP_y(1) = Max(GNEP_y(1), jMin)
      NEP_z(1) = Max(GNEP_z(1), kMin)

      NEP_x(2) = Min(GNEP_x(2), iMax)
      NEP_y(2) = Min(GNEP_y(2), jMax)
      NEP_z(2) = Min(GNEP_z(2), kMax)

      NEP_x(1) = NEP_x(1)- (iMin-1)
      NEP_x(2) = NEP_x(2)- (iMin-1)
      NEP_y(1) = NEP_y(1)- (jMin-1)
      NEP_y(2) = NEP_y(2)- (jMin-1)
      NEP_z(1) = NEP_z(1)- (kMin-1)
      NEP_z(2) = NEP_z(2)- (kMin-1)

      VIB_x(1) = VIB_x(1)- REAL(iMin-1.0_dbl , dbl)
      VIB_x(2) = VIB_x(2)- REAL(iMin-1.0_dbl , dbl)
      VIB_y(1) = VIB_y(1)- REAL(jMin-1.0_dbl , dbl)
      VIB_y(2) = VIB_y(2)- REAL(jMin-1.0_dbl , dbl)
      VIB_z(1) = VIB_z(1)- REAL(kMin-1.0_dbl , dbl)
      VIB_z(2) = VIB_z(2)- REAL(kMin-1.0_dbl , dbl)

!---- NEW: Finding the volume overlapping between particle-effetive-volume and the volume around each lattice node
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
                  Overlap(i,j,k)= MAX ( MIN(VIB_x(2),NVB_x(2)) - MAX(VIB_x(1),NVB_x(1)), 0.0_dbl) * & 
                                  MAX ( MIN(VIB_y(2),NVB_y(2)) - MAX(VIB_y(1),NVB_y(1)), 0.0_dbl) * &
                                  MAX ( MIN(VIB_z(2),NVB_z(2)) - MAX(VIB_z(1),NVB_z(1)), 0.0_dbl)
                  Overlap_sum_l= Overlap_sum_l + Overlap(i,j,k)
               END IF
            END DO
         END DO
      END DO
   END IF

!--Taking care of the Z-dir Periodic BC
   IF (GNEP_z_Per(1) .ne. GNEP_z(1)) THEN
       GNEP_z(1) = GNEP_z_Per(1)
       GNEP_z(2) = GNEP_z_Per(2)
       VIB_z(1)  = VIB_z_Per(1)
       VIB_z(2)  = VIB_z_Per(2) 
       GOTO 100
   ENDIF

   CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   CALL MPI_ALLREDUCE(Overlap_sum_l, Overlap_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)








!--Computing NB_j and Veff for each particle
!  Nbj = 0.0_dbl                                                           ! initialize Nbj - the number of moles of drug in the effective volume surrounding the particle
!  Veff = 0.0_dbl                                                          ! initialize Veff - the eff. volume of each particle
!--Solving an equation for Rj/Reff in order to estimate Veff and Nbj (see notes form July 2015)
!  CALL Find_Root(current%pardata%parid, current%pardata%bulk_conc, current%pardata%par_conc &
!                ,current%pardata%gamma_cont,current%pardata%rp,Nbj,Veff)
!  current%pardata%Veff = Veff                                             ! store Veff in particle record
!  current%pardata%Nbj = Nbj                                               ! store Nbj in particle record
!  Nbj = Nbj/zcf3  


!--Global Volume of Influence Border (VIB) for this particle
   VIB_x(1)= xp - 0.5_dbl* L_influence_P
   VIB_x(2)= xp + 0.5_dbl* L_influence_P
   VIB_y(1)= yp - 0.5_dbl* L_influence_P
   VIB_y(2)= yp + 0.5_dbl* L_influence_P
   VIB_z(1)= zp - 0.5_dbl* L_influence_P
   VIB_z(2)= zp + 0.5_dbl* L_influence_P

!--Global Nodes Effected by Particle
   GNEP_x(1)= FLOOR(VIB_x(1))
   GNEP_x(2)= CEILING(VIB_x(2))
   GNEP_y(1)= FLOOR(VIB_y(1))
   GNEP_y(2)= CEILING(VIB_y(2))
   GNEP_z(1)= FLOOR(VIB_z(1))
   GNEP_z(2)= CEILING(VIB_z(2))

!--Taking care of the Z-dir Periodic BC
   GNEP_z_Per(1) = GNEP_z(1)
   GNEP_z_Per(2) = GNEP_z(2)

   IF (GNEP_z(1) .LT. 1) THEN
       GNEP_z_Per(1) = GNEP_z(1) + nz
       GNEP_z_Per(2) = GNEP_z(2) + nz
       VIB_z_Per(1)  = VIB_z(1)  + nz
       VIB_z_Per(2)  = VIB_z(2)  + nz
   ENDIF

   IF (GNEP_z(2) .GT. nz) THEN
       GNEP_z_Per(1) = GNEP_z(1) - nz
       GNEP_z_Per(2) = GNEP_z(2) - nz
       VIB_z_Per(1)  = VIB_z(1)  - nz
       VIB_z_Per(2)  = VIB_z(2)  - nz
   ENDIF

!--Finding processor that have overlap with effective volume around the particle
200 IF((((GNEP_x(1) .GT. (iMin-1_lng)) .AND. (GNEP_x(1) .LE. iMax)) .OR. ((GNEP_x(2) .GT. (iMin-1_lng)) .AND. (GNEP_x(2) .LE. iMax))) .AND. &
      (((GNEP_y(1) .GT. (jMin-1_lng)) .AND. (GNEP_y(1) .LE. jMax)) .OR. ((GNEP_y(2) .GT. (jMin-1_lng)) .AND. (GNEP_y(2) .LE. jMax))) .AND. &
      (((GNEP_z(1) .GT. (kMin-1_lng)) .AND. (GNEP_z(1) .LE. kMax)) .OR. ((GNEP_z(2) .GT. (kMin-1_lng)) .AND. (GNEP_z(2) .LE. kMax)))  )THEN

      NEP_x(1) = Max(GNEP_x(1) , 1)
      NEP_y(1) = Max(GNEP_y(1) , 1)
      NEP_z(1) = Max(GNEP_z(1) , 1)

      NEP_x(2) = Min(GNEP_x(2) , iMax )
      NEP_y(2) = Min(GNEP_y(2) , jMax )
      NEP_z(2) = Min(GNEP_z(2) , kMax )

      NEP_x(1) = NEP_x(1) - (iMin-1)
      NEP_x(2) = NEP_x(2) - (iMin-1)
      NEP_y(1) = NEP_y(1) - (jMin-1)
      NEP_y(2) = NEP_y(2) - (jMin-1)
      NEP_z(1) = NEP_z(1) - (kMin-1)
      NEP_z(2) = NEP_z(2) - (kMin-1)

!-----Computing particle release contribution to scalar field at each lattice node
      DO i= NEP_x(1),NEP_x(2)
         DO j= NEP_y(1),NEP_y(2)
            DO k= NEP_z(1),NEP_z(2)
               IF (node(i,j,k) .EQ. FLUID) THEN                 
                  IF (Overlap_sum .GT. 1e-40) THEN 			! Overlap_sum going to zero when the particle is disapearing
                     Overlap(i,j,k) = Overlap(i,j,k) / Overlap_sum
                  ELSE
                     Overlap(i,j,k) = 0.0
                  END IF
                   
	          delphi_particle(i,j,k)  = delphi_particle(i,j,k)  + current%pardata%delNBbyCV * Overlap(i,j,k) 
!                 tausgs_particle_x(i,j,k)= tausgs_particle_x(i,j,k)- current%pardata%up*Nbj   * (Overlap(i,j,k)/Overlap_sum)
!                 tausgs_particle_y(i,j,k)= tausgs_particle_y(i,j,k)- current%pardata%up*Nbj   * (Overlap(i,j,k)/Overlap_sum)
!                 tausgs_particle_z(i,j,k)= tausgs_particle_z(i,j,k)- current%pardata%up*Nbj   * (Overlap(i,j,k)/Overlap_sum)
               END IF 
            END DO
         END DO
      END DO
   END IF

   IF (GNEP_z_Per(1) .ne. GNEP_z(1)) THEN
       GNEP_z(1) = GNEP_z_Per(1)
       GNEP_z(2) = GNEP_z_Per(2)
       VIB_z(1)  = VIB_z_Per(1)
       VIB_z(2)  = VIB_z_Per(2)
       GOTO 200
   ENDIF

   current => next
ENDDO

!===================================================================================================
END SUBROUTINE Particle_Drug_To_Nodes  		 
!===================================================================================================





!===================================================================================================
SUBROUTINE Find_Root(parid,conc,cs,gammaj,Rj,Nbj,Veff)
!===================================================================================================
IMPLICIT NONE
INTEGER(lng)   :: iter,nmax
REAL(dbl),intent(in) :: conc,cs,gammaj,Rj
INTEGER(lng),intent(in) :: parid
REAL(dbl),intent(out) :: Nbj,Veff
REAL(dbl)      :: xnew,xold,f,fprime,error,Reff,parcb,parcs
REAL(dbl) :: Aj,Bj,a,b,c,AjBj
REAL(dbl),parameter :: eps = 1.0e-12_dbl, tol = 1.0e-8_dbl, largenum = 1.0e8_dbl

parcb = conc
parcs = cs
nmax = 100_lng
iter = 0_lng
xnew = 0.0_dbl
xold = 0.5_dbl

!---The conc values are quite small. So we will make them larger so that the
!--- coeffs in the eq for Rj/Reff are not small.
parcb = parcb*largenum
parcs = parcs*largenum

Aj = (parcb-gammaj*parcs)/(1.0_dbl-gammaj)
Bj = (parcs - parcb)/max((parcb-gammaj*cs),eps)
AjBj = (parcs-parcb)/(1.0_dbl-gammaj)

a = Aj + 1.5_dbl*AjBj!*Aj*Bj
b = -1.5_dbl*AjBj!Aj*Bj
c = parcb - Aj
error = abs(xold-xnew)

DO WHILE ((error.gt.tol).AND.(iter.LE.nmax))
   f = a*(xold**3.0_dbl)+b*xold+c
   fprime = 3.0_dbl*a*(xold**2.0_dbl)+b
   if (fprime.GE.0.0_dbl) then
      xnew = xold - (f/max(fprime,1.0_dbl*eps))
   else
      !xnew = xold - (f/min(fprime,-1.0_dbl*eps))
       xnew = xold - (f/fprime)
   endif
   error = abs(xnew-xold)
   iter = iter+1_lng
   xold = xnew
END DO

xnew= max(min(xnew,1.0_dbl),0.01) 			! Limit xnew (radius ratio) to values that are meaningful and not too small or large. 
Reff= Rj/xnew
Veff= (88.0_dbl/21.0_dbl)*(Reff**3.0_dbl)
Nbj = conc*Veff
!===================================================================================================
END SUBROUTINE Find_Root
!===================================================================================================




!================================================
END MODULE ParticleDrug 
!================================================
