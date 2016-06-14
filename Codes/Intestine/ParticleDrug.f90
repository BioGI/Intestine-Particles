!==================================================================================================
MODULE ParticleDrug				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================
USE SetPrecision
USE Setup
USE IC
USE MPI

IMPLICIT NONE

CONTAINS


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
REAL(dbl)       	 :: N_b
REAL(dbl)    	         :: R_P, Sh_P, delta_P
REAl(dbl)                :: R_influence_p, L_influence_p		! Parameters related to particle's volume of influence
REAl(dbl)                :: V_influence_P	 			! Parameters related to particle's volume of influence
REAL(dbl)		 :: Cb_Total_Veff_l, Cb_Total_Veff
REAL(dbl),DIMENSION(2)   :: GVIB_x, GVIB_y, GVIB_z, GVIB_z_Per 		! Volume of Influence's Borders
REAL(dbl),DIMENSION(2)   :: NVB_x, NVB_y, NVB_z				! Node Volume's Borders
REAL(dbl)                :: Delta_L
REAL(dbl)                :: x_DP, y_DP, z_DP				! Coordinates of "Discretized Point" (DP)
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

delta_mesh = 1.0_dbl
zcf3 = 1.0_dbl

current => ParListHead%next

DO WHILE (ASSOCIATED(current))

next => current%next
IF (current%pardata%rp .GT. Min_R_Acceptable) THEN	

!--Particle length scale: delta= R/Sh & effective radius: R_influence_P= R+(N_b*delta)
   N_b = 2.0
   R_P = current%pardata%rp
   Sh_P= current%pardata%sh
   delta_P= R_P/Sh_P
   R_influence_P= (R_P+N_b*delta_P)/xcf

!--Computing equivalent cubic mesh length scale
   V_influence_P= (4.0_dbl/3.0_dbl)*PI* R_influence_P**3.0_dbl
   L_influence_P= V_influence_P **(1.0_dbl/3.0_dbl)
   V_eff_Ratio  = V_influence_P/zcf3 					! Ratio of the effective volume to cell size 

   Cb_Total_Veff_l  = 0.0_lng
   Cb_Total_Veff    = 0.0_lng
   NumFluids_Veff_l = 0.0_lng
   NumFluids_Veff   = 0.0_lng

!----------------------------------------------------------------------------------------------------------------------
!--Veff is smaller than the mesh volume --> Cb = Trilinear interpolation of the concentration at particle location
!--No communication is necessary between processors
!----------------------------------------------------------------------------------------------------------------------
   IF (V_eff_Ratio .LE. 1.0) THEN 					
      CaseNo= 1
      IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
!------- TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
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
!------- Interpolation in x-direction
         c00 = phi(ix0,iy0,iz0) * (1.0_dbl-xd) + phi(ix1,iy0,iz0) * xd
         c01 = phi(ix0,iy0,iz1) * (1.0_dbl-xd) + phi(ix1,iy0,iz1) * xd
         c10 = phi(ix0,iy1,iz0) * (1.0_dbl-xd) + phi(ix1,iy1,iz0) * xd
         c11 = phi(ix0,iy1,iz1) * (1.0_dbl-xd) + phi(ix1,iy1,iz1) * xd
!------- Interpolation in y-direction
         c0  = c00 * (1.0_dbl-yd) + c10 * yd
         c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------- Interpolation in z-direction
         c   = c0 * (1.0_dbl-zd) + c1 * zd
         Cb_Hybrid= c 
         current%pardata%bulk_conc = Cb_Hybrid
      END IF !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   END IF 





!-----------------------------------------------------------------------------------------------------------------------
!--Veff is larger than 1 which means parallel communication might be necessary 
!-----------------------------------------------------------------------------------------------------------------------
   IF (V_eff_Ratio .GT. 1.0) THEN

!-----Finding particle global location (in the whole domain not this particular processor)
      xp= current%pardata%xp 
      yp= current%pardata%yp 
      zp= current%pardata%zp 

!-----Global Volume of Influence Border (GVIB) for this particle (in the whole domain not this particular processor)
      GVIB_x(1)= xp - 0.5_dbl * L_influence_P
      GVIB_x(2)= xp + 0.5_dbl * L_influence_P
      GVIB_y(1)= yp - 0.5_dbl * L_influence_P
      GVIB_y(2)= yp + 0.5_dbl * L_influence_P
      GVIB_z(1)= zp - 0.5_dbl * L_influence_P
      GVIB_z(2)= zp + 0.5_dbl * L_influence_P

!-----TAKING CARE OF THE PERIODIC BC
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

!-----Finding processor that have overlap with effective volume around the particle       
300   IF ((((GVIB_x(1) .GE. (iMin-1_lng)) .AND. (GVIB_x(1) .LT. iMax)) .OR. ((GVIB_x(2) .GE. (iMin-1_lng)) .AND. (GVIB_x(2) .LT. iMax))) .AND. &
          (((GVIB_y(1) .GE. (jMin-1_lng)) .AND. (GVIB_y(1) .LT. jMax)) .OR. ((GVIB_y(2) .GE. (jMin-1_lng)) .AND. (GVIB_y(2) .LT. jMax))) .AND. &
          (((GVIB_z(1) .GE. (kMin-1_lng)) .AND. (GVIB_z(1) .LT. kMax)) .OR. ((GVIB_z(2) .GE. (kMin-1_lng)) .AND. (GVIB_z(2) .LT. kMax)))  )THEN

!-------------------------------------------------------------------------------------------------------------------------
!------- Veff is slightly larger than lattice cell volume --> Volume of influence is discretized to provide 27 points
!------- Cb= Average of concentration interpolated on each of the descritized nodes inside volume of influence
!-------------------------------------------------------------------------------------------------------------------------
         IF ((V_eff_Ratio .GT. 1.0) .AND. (V_eff_Ratio .LT. 27.0)) THEN		
             CaseNo= 2

!------------Discretizing the volume of influence to  make sure at least 27 points are available
             Delta_L = (GVIB_x(2)-GVIB_x(1)) / 2.0 

!------------Loop over discretized points and averaging the concentration
             DO i= 0, 2
                DO j= 0, 2
                   DO k= 0, 2
                      x_DP= GVIB_x(1) + (i * Delta_L) 
                      y_DP= GVIB_y(1) + (j * Delta_L)
                      z_DP= GVIB_z(1) + (k * Delta_L)
                      IF ((x_DP .GE. (REAL(iMin,dbl)-1.0_dbl)) .AND. &               
		          (x_DP .LT.  REAL(iMax,dbl)         ) .AND. &
                          (y_DP .GE. (REAL(jMin,dbl)-1.0_dbl)) .AND. & 
			  (y_DP .LT.  REAL(jMax,dbl)         ) .AND. &
                          (z_DP .GE. (REAL(kMin,dbl)-1.0_dbl)) .AND. & 
			  (z_DP .LT.  REAL(kMax,dbl)         ) ) THEN

!-------------------------Finding Local lattice nodes surrounding this point (This point is discretized and is not a lattice node))
                          ix0 = FLOOR(x_DP)   - (REAL(iMin,dbl)-1.0_dbl)
                          ix1 = CEILING(x_DP) - (REAL(iMin,dbl)-1.0_dbl)
                          iy0 = FLOOR(y_DP)   - (REAL(jMin,dbl)-1.0_dbl)
                          iy1 = CEILING(y_DP) - (REAL(jMin,dbl)-1.0_dbl) 
                          iz0 = FLOOR(z_DP)   - (REAL(kMin,dbl)-1.0_dbl)
                          iz1 = CEILING(z_DP) - (REAL(kMin,dbl)-1.0_dbl)
                              
            		  x_DP = x_DP - REAL(iMin-1_lng,dbl)
                          y_DP = y_DP - REAL(jMin-1_lng,dbl)
                          z_DP = z_DP - REAL(kMin-1_lng,dbl)
 
!-------------------------TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
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

!-------------------------Interpolation in x-direction
                          c00 = phi(ix0,iy0,iz0) * (1.0_dbl-xd) + phi(ix1,iy0,iz0) * xd
                          c01 = phi(ix0,iy0,iz1) * (1.0_dbl-xd) + phi(ix1,iy0,iz1) * xd
                          c10 = phi(ix0,iy1,iz0) * (1.0_dbl-xd) + phi(ix1,iy1,iz0) * xd
                          c11 = phi(ix0,iy1,iz1) * (1.0_dbl-xd) + phi(ix1,iy1,iz1) * xd
!------------------------ Interpolation in y-direction
                          c0  = c00 * (1.0_dbl-yd) + c10 * yd
                          c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------------------------ Interpolation in z-direction
                          c   = c0 * (1.0_dbl-zd) + c1 * zd

                          Cb_Total_Veff_l  = Cb_Total_Veff_l  + c
                          NumFluids_Veff_l = NumFluids_Veff_l + 1_lng
                      END IF
                   END DO
               END DO
            END DO
          
!----------------------------------------------------------------------------------------------------------------------
!--------Veff is much larger than mesh volume --> Cb= total number of moles in volume of influence / volume of influence 
!----------------------------------------------------------------------------------------------------------------------
         ELSE IF (V_eff_Ratio .GE. 27.0) THEN                             
             CaseNo = 3
!------------Finding the lattice "Nodes Effected by Particle"
             GNEP_x(1)= CEILING(GVIB_x(1))
             GNEP_y(1)= CEILING(GVIB_y(1))
             GNEP_z(1)= CEILING(GVIB_z(1))
             GNEP_x(2)= FLOOR  (GVIB_x(2))
             GNEP_y(2)= FLOOR  (GVIB_y(2))
             GNEP_z(2)= FLOOR  (GVIB_z(2))

!------------Finding the lattice "Nodes Effected by Particle"
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
                         Cb_Total_Veff_l  = Cb_Total_Veff_l  + phi(i,j,k)
                         NumFluids_Veff_l = NumFluids_Veff_l + 1_lng
                      END IF
                   END DO
                END DO
             END DO

         END IF  									! Conditional for cases 2 and 3
     END IF 										! Conditional for the processor which has overlap with effective volume 




!----TAKING CARE OF THE PERIODIC BC
     IF (GVIB_z_Per(1) .NE. GVIB_z(1)) THEN
        GVIB_z(1)= GVIB_z_Per(1) 
        GVIB_z(2)= GVIB_z_Per(2) 
        GO TO 300
     ENDIF

!----Communication with other processors for V_eff greater than 1
     CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
     CALL MPI_ALLREDUCE(Cb_Total_Veff_l , Cb_Total_Veff , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
     CALL MPI_ALLREDUCE(NumFluids_Veff_l, NumFluids_Veff, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, mpierr)
     Cb_Hybrid= Cb_Total_Veff / NumFluids_Veff
     current%pardata%bulk_conc = Cb_Hybrid
	
   END IF       			                                    		!End of conditional for V_eff greater than 1 
      
ENDIF  

!open(172,file='Cb-'//sub//'.dat', position='append')
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

   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                                           !only calculate the drug release when particle radius is larger than 0.1 micron
      current%pardata%rpold = current%pardata%rp
      bulkconc = current%pardata%bulk_conc
      temp = current%pardata%rpold**2.0_dbl-4.0_dbl*tcf*molarvol*diffm*current%pardata%sh*max((current%pardata%par_conc-bulkconc),0.0_dbl)
      IF (temp.GE.0.0_dbl) THEN
         current%pardata%rp= 0.5_dbl*(current%pardata%rpold+sqrt(temp))
      ELSE
         temp = 0.0_dbl
         current%pardata%rp= 0.5_dbl*(current%pardata%rpold+sqrt(temp))
      END IF
      deltaR=current%pardata%rpold-current%pardata%rp
      current%pardata%delNBbyCV = (4.0_dbl/3.0_dbl) * PI*(current%pardata%rpold**3.0_dbl - current%pardata%rp**3.0_dbl) /(molarvol*zcf3)
   ELSE IF ((current%pardata%rp .LT. Min_R_Acceptable) .AND. (current%pardata%rp .NE. 0.0)) THEN
      current%pardata%rp= 0.0_dbl
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
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 

   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                                           !only calculate the drug release when particle radius is larger than 0.1 micron
      IF (mySub .EQ.current%pardata%cur_part) THEN
         current%pardata%sh= 1.0_dbl + (current%pardata%gamma_cont / (1.0_dbl-current%pardata%gamma_cont)) 
         S= current%pardata%S
         Sst= S* (current%pardata%rp**2.0) / diffm
         current%pardata%Sst= Sst
         IF (Sst.LT.5.0_dbl) THEN
            current%pardata%sh = current%pardata%sh + 0.296_dbl*(Sst**0.5_dbl)
         ELSE
            Sh0 = exp(0.162_dbl + 0.202_dbl*log(Sst) - 7.5e-6_dbl*(log(Sst)**5.4_dbl)) 
            current%pardata%sh = current%pardata%sh + Sh0-1.0_dbl
         END IF
      END IF

      RANK= current%pardata%cur_part - 1
      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      CALL MPI_BCast(current%pardata%sh,1,MPI_DOUBLE_PRECISION, RANK, MPI_COMM_WORLD,mpierr)
   END IF 

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

   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN						!only calculate the drug release when particle radius is larger than 0.1 micron				
      IF (mySub .EQ.current%pardata%cur_part) THEN
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
!!!!!!!! MAKE SURE THE ABOVE NODES ARE FLUID NODES

         dwdz = w(it,jt,kt) - w(ib,jb,kb)
         S = abs(dwdz*vcf/zcf)
         current%pardata%S = S
       END IF 
    END IF

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
REAL(dbl),DIMENSION(2)    :: GVIB_x, GVIB_y, GVIB_z, GVIB_z_Per 	! Global Volume of Influence's Borders (in whole domain)
REAL(dbl),DIMENSION(2)    :: LVIB_x, LVIB_y, LVIB_z                     ! Local  Volume of Influence's Borders (in current procesor) 

REAL(dbl),DIMENSION(2)    :: NVB_x, NVB_y, NVB_z			! Node Volume's Borders
INTEGER  ,DIMENSION(2)    :: LN_x,  LN_y,  LN_z				! Lattice Nodes Surronding the particle
INTEGER  ,DIMENSION(2)    :: GNEP_x, GNEP_y, GNEP_z, GNEP_z_Per         ! Lattice Nodes Surronding the particle (Global: not considering the partitioning for parallel processing)
INTEGER  ,DIMENSION(2)    :: NEP_x,   NEP_y,  NEP_z                     ! Lattice Nodes Surronding the particle (Local: in current processor)
REAL(dbl)		  :: tmp, Overlap_sum_l, Overlap_sum
REAL(dbl)   		  ::  Overlap_test,Overlap_test_Global
TYPE(ParRecord), POINTER  :: current
TYPE(ParRecord), POINTER  :: next

delta_mesh = 1.0_dbl
zcf3 = xcf*ycf*zcf
current => ParListHead%next

DO WHILE (ASSOCIATED(current))
   next => current%next 

   IF (current%pardata%rp .GT. Min_R_Acceptable) THEN                   !only calculate the drug release when particle radius is larger than 0.1 micron

!--Calculate length scale for jth particle:  delta = R / Sh
!--Calculate effective radius: R_influence_P = R + (N_d *delta)
!--Note: need to convert this into Lattice units and not use the physical length units
!--Then compute equivalent cubic mesh length scale
   N_d = 3.0
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

!--Global Volume of Influence Border (GVIB) for this particle
   GVIB_x(1)= xp - 0.5_dbl * L_influence_P
   GVIB_x(2)= xp + 0.5_dbl * L_influence_P
   GVIB_y(1)= yp - 0.5_dbl * L_influence_P
   GVIB_y(2)= yp + 0.5_dbl * L_influence_P
   GVIB_z(1)= zp - 0.5_dbl * L_influence_P
   GVIB_z(2)= zp + 0.5_dbl * L_influence_P

!--Global Nodes Effected by Particle 
   GNEP_x(1)= FLOOR(GVIB_x(1))
   GNEP_y(1)= FLOOR(GVIB_y(1))
   GNEP_z(1)= FLOOR(GVIB_z(1))
   GNEP_x(2)= CEILING(GVIB_x(2))
   GNEP_y(2)= CEILING(GVIB_y(2))
   GNEP_z(2)= CEILING(GVIB_z(2))

!--Taking care of the Z-dir Periodic BC
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

   Overlap_sum_l = 0
   Overlap       = 0.0_dbl

!--Finding processors with overlap with effective volume around the particle       
100 IF((((GNEP_x(1) .GT. (iMin-1_lng)) .AND. (GNEP_x(1) .LE. iMax)) .OR. ((GNEP_x(2) .GT. (iMin-1_lng)) .AND. (GNEP_x(2) .LE. iMax))) .AND. &
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
                  Overlap(i,j,k)= MAX ( MIN(LVIB_x(2),NVB_x(2)) - MAX(LVIB_x(1),NVB_x(1)), 0.0_dbl) * & 
                                  MAX ( MIN(LVIB_y(2),NVB_y(2)) - MAX(LVIB_y(1),NVB_y(1)), 0.0_dbl) * &
                                  MAX ( MIN(LVIB_z(2),NVB_z(2)) - MAX(LVIB_z(1),NVB_z(1)), 0.0_dbl)
		  Overlap(i,j,k) = Overlap(i,j,k) * (max((Cs_mol-phi(i,j,k) ),0.0_dbl) / Cs_mol)
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
       GVIB_z(1) = GVIB_z_Per(1)
       GVIB_z(2) = GVIB_z_Per(2) 
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
   GVIB_x(1)= xp - 0.5_dbl* L_influence_P
   GVIB_x(2)= xp + 0.5_dbl* L_influence_P
   GVIB_y(1)= yp - 0.5_dbl* L_influence_P
   GVIB_y(2)= yp + 0.5_dbl* L_influence_P
   GVIB_z(1)= zp - 0.5_dbl* L_influence_P
   GVIB_z(2)= zp + 0.5_dbl* L_influence_P

!--Global Nodes Effected by Particle
   GNEP_x(1)= FLOOR(GVIB_x(1))
   GNEP_y(1)= FLOOR(GVIB_y(1))
   GNEP_z(1)= FLOOR(GVIB_z(1))
   GNEP_x(2)= CEILING(GVIB_x(2))
   GNEP_y(2)= CEILING(GVIB_y(2))
   GNEP_z(2)= CEILING(GVIB_z(2))

!--Taking care of the Z-dir Periodic BC
   GNEP_z_Per(1) = GNEP_z(1)
   GNEP_z_Per(2) = GNEP_z(2)

   IF (GNEP_z(1) .LT. 1) THEN
       GNEP_z_Per(1) = GNEP_z(1) + nz
       GNEP_z_Per(2) = GNEP_z(2) + nz
   ENDIF

   IF (GNEP_z(2) .GT. nz) THEN
       GNEP_z_Per(1) = GNEP_z(1) - nz
       GNEP_z_Per(2) = GNEP_z(2) - nz
   ENDIF

!--Finding processor that have overlap with effective volume around the particle

OVERLAP_TEST = 0.0_dbl 

200 IF((((GNEP_x(1) .GT. (iMin-1_lng)) .AND. (GNEP_x(1) .LE. iMax)) .OR. ((GNEP_x(2) .GT. (iMin-1_lng)) .AND. (GNEP_x(2) .LE. iMax))) .AND. &
      (((GNEP_y(1) .GT. (jMin-1_lng)) .AND. (GNEP_y(1) .LE. jMax)) .OR. ((GNEP_y(2) .GT. (jMin-1_lng)) .AND. (GNEP_y(2) .LE. jMax))) .AND. &
      (((GNEP_z(1) .GT. (kMin-1_lng)) .AND. (GNEP_z(1) .LE. kMax)) .OR. ((GNEP_z(2) .GT. (kMin-1_lng)) .AND. (GNEP_z(2) .LE. kMax)))  )THEN

      NEP_x(1) = Max(GNEP_x(1), 1)    - (iMin-1)
      NEP_y(1) = Max(GNEP_y(1), 1)    - (jMin-1)
      NEP_z(1) = Max(GNEP_z(1), 1)    - (kMin-1)
      NEP_x(2) = Min(GNEP_x(2), iMax) - (iMin-1)
      NEP_y(2) = Min(GNEP_y(2), jMax) - (jMin-1)
      NEP_z(2) = Min(GNEP_z(2), kMax) - (kMin-1)

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
                  OVERLAP_TEST= OVERLAP_TEST + Overlap(i,j,k)

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
       GOTO 200
   ENDIF

   CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   CALL MPI_ALLREDUCE(Overlap_test, Overlap_test_Global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)

   IF (abs(Overlap_test_Global - 1.0) .GT. 0.5) THEN                        ! Detecting the case of overlap = 0.0
      current%pardata%rp =  (current%pardata%rp**3 + current%pardata%delNBbyCV * (molarvol*zcf3) * (3/(4*PI)) )**(1.0_dbl/3.0_dbl)
      current%pardata%delNBbyCV = 0.0_dbl
   END IF

 END IF 						! Condition to check if R > R_min_acceptable
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
