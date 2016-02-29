!==================================================================================================
MODULE LBM				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================

USE SetPrecision
USE Setup
USE ICBC
USE MPI

IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------------------------------------
SUBROUTINE LBM_Setup	! setup the LBM simulation
!---------------------------------------------------------------------------------------------------

IMPLICIT NONE

!---- Initialize variables and arrays---------------------------------------------------------------
f     = 0.0_dbl			! distribution functions
fplus = 0.0_dbl			! post-collision distribution functions
u     = 0.0_dbl			! x-velocity
v     = 0.0_dbl			! y-velocity
w     = 0.0_dbl			! z-velocity
rho   = 0.0_lng			! density
ex    = 0.0_dbl			! lattice discretized velocity vector (x-component)
ey    = 0.0_dbl			! lattice discretized velocity vector (x-component)
ez    = 0.0_dbl			! lattice discretized velocity vector (x-component)
bb    = 0_lng			! bounceback directions
wt    = 0.0_dbl			! weighting coefficients

!---- Fill out weighting coefficient array ---------------------------------------------------------
wt(0)    = 2.0_dbl/9.0_dbl	
wt(1:6)  = 1.0_dbl/9.0_dbl
wt(7:14) = 1.0_dbl/72.0_dbl

!---- Fill out bounceback array --------------------------------------------------------------------
bb(0)  = 0_lng
bb(1)  = 2_lng
bb(2)  = 1_lng
bb(3)  = 4_lng
bb(4)  = 3_lng
bb(5)  = 6_lng
bb(6)  = 5_lng
bb(7)  = 8_lng
bb(8)  = 7_lng
bb(9)  = 10_lng
bb(10) = 9_lng
bb(11) = 12_lng
bb(12) = 11_lng
bb(13) = 14_lng
bb(14) = 13_lng

!----Fill out symmetry array -----------------------------------------------------------------------
! iComm=2, -ZY FACE
sym(0,2)  = 0_lng
sym(1,2)  = 2_lng
sym(2,2)  = 1_lng
sym(3,2)  = 4_lng
sym(4,2)  = 3_lng
sym(5,2)  = 6_lng
sym(6,2)  = 5_lng
sym(7,2)  = 11_lng
sym(8,2)  = 12_lng
sym(9,2)  = 14_lng
sym(10,2) = 13_lng
sym(11,2) = 7_lng
sym(12,2) = 8_lng
sym(13,2) = 10_lng
sym(14,2) = 9_lng

! iComm=3, -ZX FACE
sym(0,4)  = 0_lng
sym(1,4)  = 2_lng
sym(2,4)  = 1_lng
sym(3,4)  = 4_lng
sym(4,4)  = 3_lng
sym(5,4)  = 6_lng
sym(6,4)  = 5_lng
sym(7,4)  = 13_lng
sym(8,4)  = 14_lng
sym(9,4)  = 12_lng
sym(10,4) = 11_lng
sym(11,4) = 10_lng
sym(12,4) = 9_lng
sym(13,4) = 7_lng
sym(14,4) = 8_lng

! iComm=8, Z AXIS
sym(0,8)  = 0_lng
sym(1,8)  = 2_lng
sym(2,8)  = 1_lng
sym(3,8)  = 4_lng
sym(4,8)  = 3_lng
sym(5,8)  = 6_lng
sym(6,8)  = 5_lng
sym(7,8)  = 10_lng
sym(8,8)  = 9_lng
sym(9,8)  = 8_lng
sym(10,8) = 7_lng
sym(11,8) = 13_lng
sym(12,8) = 14_lng
sym(13,8) = 11_lng
sym(14,8) = 12_lng 

!---- Fill velocity direction vector arrays --------------------------------------------------------
ex(0) =	 0.0_dbl		! direction 0
ey(0) =	 0.0_dbl
ez(0) =	 0.0_dbl
ex(1) =	 1.0_dbl		! direction 1
ey(1) =	 0.0_dbl
ez(1) =	 0.0_dbl
ex(2) =  -1.0_dbl		! direction 2
ey(2) =   0.0_dbl
ez(2) =   0.0_dbl
ex(3) =   0.0_dbl		! direction 3
ey(3) =   1.0_dbl
ez(3) =   0.0_dbl
ex(4) =   0.0_dbl		! direction 4
ey(4) =  -1.0_dbl
ez(4) =   0.0_dbl
ex(5) =   0.0_dbl		! direction 5
ey(5) =   0.0_dbl
ez(5) =   1.0_dbl
ex(6) =   0.0_dbl		! direction 6
ey(6) =   0.0_dbl
ez(6) =  -1.0_dbl
ex(7) =   1.0_dbl		! direction 7
ey(7) =   1.0_dbl
ez(7) =   1.0_dbl
ex(8) =  -1.0_dbl		! direction 8
ey(8) =  -1.0_dbl
ez(8) =  -1.0_dbl
ex(9) =   1.0_dbl		! direction 9
ey(9) =   1.0_dbl
ez(9) =  -1.0_dbl
ex(10) = -1.0_dbl		! direction 10
ey(10) = -1.0_dbl
ez(10) =  1.0_dbl
ex(11) = -1.0_dbl		! direction 11
ey(11) =  1.0_dbl
ez(11) =  1.0_dbl
ex(12) =  1.0_dbl		! direction 12
ey(12) = -1.0_dbl
ez(12) = -1.0_dbl
ex(13) =  1.0_dbl		! direction 13
ey(13) = -1.0_dbl
ez(13) =  1.0_dbl
ex(14) = -1.0_dbl		! direction 14
ey(14) =  1.0_dbl
ez(14) = -1.0_dbl

!---- Define other simulation parameters -----------------------------------------------------------
nuL   		= (2.0_dbl*tau - 1.0_dbl)/6.0_dbl	! lattice kinematic viscosity
denL 		= 1.0_dbl				! arbitrary lattice density (1.0 for convenience)
oneOVERtau 	= 1.0_dbl/tau				! reciprical of tau
cs		= (1.0_dbl)/(SQRT(3.0_dbl))		! speed of sound on the lattice

!---- Initialize timestep --------------------------------------------------------------------------
iter = 0_lng						! intialize the starting timestep to 0 - will get reset in 'ICs' in ICBCM.f90

!---- Calculate feq for initial condition ----------------------------------------------------------
CALL Equilibrium

!---- Set f-from wall motion sums to zero at initial timestep --------------------------------------
fmovingsum = 0.0_dbl
fmovingrhosum = 0.0_dbl

!===================================================================================================
END SUBROUTINE LBM_Setup
!===================================================================================================







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
SUBROUTINE Interp_Parvel ! Using Trilinear interpolation
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
! Do first level linear interpolation in x-direction
	c00 = u(ix0,iy0,iz0)*(1.0_dbl-xd)+u(ix1,iy0,iz0)*xd	
	c01 = u(ix0,iy0,iz1)*(1.0_dbl-xd)+u(ix1,iy0,iz1)*xd	
	c10 = u(ix0,iy1,iz0)*(1.0_dbl-xd)+u(ix1,iy1,iz0)*xd	
	c11 = u(ix0,iy1,iz1)*(1.0_dbl-xd)+u(ix1,iy1,iz1)*xd	
	
! Do second level linear interpolation in y-direction
	c0  = c00*(1.0_dbl-yd)+c10*yd
	c1  = c01*(1.0_dbl-yd)+c11*yd

! Do third level linear interpolation in z-direction
	c   = c0*(1.0_dbl-zd)+c1*zd
        current%pardata%up=c


! v-interpolation
! Do first level linear interpolation in x-direction
	c00 = v(ix0,iy0,iz0)*(1.0_dbl-xd)+v(ix1,iy0,iz0)*xd
	c01 = v(ix0,iy0,iz1)*(1.0_dbl-xd)+v(ix1,iy0,iz1)*xd
	c10 = v(ix0,iy1,iz0)*(1.0_dbl-xd)+v(ix1,iy1,iz0)*xd
	c11 = v(ix0,iy1,iz1)*(1.0_dbl-xd)+v(ix1,iy1,iz1)*xd	

! Do second level linear interpolation in y-direction
	c0  = c00*(1.0_dbl-yd)+c10*yd
	c1  = c01*(1.0_dbl-yd)+c11*yd

! Do third level linear interpolation in z-direction
	c   = c0*(1.0_dbl-zd)+c1*zd
        current%pardata%vp=c

! w-interpolation
! Do first level linear interpolation in x-direction
	c00 = w(ix0,iy0,iz0)*(1.0_dbl-xd)+w(ix1,iy0,iz0)*xd	
	c01 = w(ix0,iy0,iz1)*(1.0_dbl-xd)+w(ix1,iy0,iz1)*xd	
	c10 = w(ix0,iy1,iz0)*(1.0_dbl-xd)+w(ix1,iy1,iz0)*xd	
	c11 = w(ix0,iy1,iz1)*(1.0_dbl-xd)+w(ix1,iy1,iz1)*xd	

! Do second level linear interpolation in y-direction
	c0  = c00*(1.0_dbl-yd)+c10*yd
	c1  = c01*(1.0_dbl-yd)+c11*yd

! Do third level linear interpolation in z-direction
	c   = c0*(1.0_dbl-zd)+c1*zd
        current%pardata%wp=c
      END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! point to next node in the list
	current => next

ENDDO

!===================================================================================================
END SUBROUTINE Interp_Parvel ! Using Trilinear interpolation
!===================================================================================================









!===================================================================================================
SUBROUTINE Interp_bulkconc(Cb_Local) ! Using Trilinear interpolation
!===================================================================================================

IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd
REAL(dbl)     :: xp,yp,zp
REAL(dbl)     :: Cb_Local
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

	!yd=0.0_dbl ! TEST: used to keep particle motion plainly 2-D ! Balaji added

! phi-interpolation
! Do first level linear interpolation in x-direction
	c00 = phi(ix0,iy0,iz0)*(1.0_dbl-xd)+phi(ix1,iy0,iz0)*xd	
	c01 = phi(ix0,iy0,iz1)*(1.0_dbl-xd)+phi(ix1,iy0,iz1)*xd	
	c10 = phi(ix0,iy1,iz0)*(1.0_dbl-xd)+phi(ix1,iy1,iz0)*xd	
	c11 = phi(ix0,iy1,iz1)*(1.0_dbl-xd)+phi(ix1,iy1,iz1)*xd	

! Do second level linear interpolation in y-direction
	c0  = c00*(1.0_dbl-yd)+c10*yd
	c1  = c01*(1.0_dbl-yd)+c11*yd

! Do third level linear interpolation in z-direction
	c   = c0*(1.0_dbl-zd)+c1*zd
!       current%pardata%bulk_conc=c
        Cb_Local= c        

      END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! point to next node in the list
	current => next
ENDDO

!===================================================================================================
END SUBROUTINE Interp_bulkconc  ! Using Trilinear interpolation
!===================================================================================================








!===================================================================================================
SUBROUTINE Calc_Global_Bulk_Scalar_Conc(Cb_Domain)  !Calculate Global Bulk SCalar COnc for use in the scalar drug relese model 
!===================================================================================================

IMPLICIT NONE
INTEGER(lng)  :: i,j,k
REAL(dbl)     :: Cb_Domain
! Calculate the bulk Conc = total number of moles/total domain size or it is the average conc in the domain
Cb_global = 0.0_dbl
Cb_numFluids = 0_lng
DO k=1,nzSub
  DO j=1,nySub
    DO i=1,nxSub
      IF(node(i,j,k) .EQ. FLUID) THEN
        !if (mySub.eq.22) then
        !        write(*,*) iter,i,j,k,Cb_global,phi(i,j,k),delphi_particle(i,j,k)
        !end if
        Cb_global = Cb_global + phi(i,j,k)
        Cb_numFluids = Cb_numFluids + 1_lng
      END IF
    END DO
  END DO
END DO

Cb_Domain = Cb_global/ Cb_numFluids

!===================================================================================================
END SUBROUTINE Calc_Global_Bulk_Scalar_Conc
!===================================================================================================








!===================================================================================================
SUBROUTINE Compute_Cb(V_eff_Ratio,CaseNo,Cb_Hybrid) ! Computes the mesh-independent bulk concentration
!===================================================================================================

IMPLICIT NONE

INTEGER(lng)  			:: i,j,k, kk, CaseNo
INTEGER(lng)  			:: ix0,ix1,iy0,iy1,iz0,iz00,iz1,iz11			! Trilinear interpolation parameters
INTEGER(lng)			:: NumFluids_Veff = 0_lng
INTEGER,DIMENSION(2)   		:: LN_x,  LN_y,  LN_z				! Lattice Nodes Surronding the particle
INTEGER,DIMENSION(2)   		:: NEP_x, NEP_y, NEP_z                  	! Lattice Nodes Surronding the particle

REAL(dbl)     			:: c00,c01,c10,c11,c0,c1,c,xd,yd,zd		! Trilinear interpolation parameters
REAL(dbl)  	   		:: xp,yp,zp
REAL(dbl)			:: delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
REAL(dbl)       	        :: N_b         					! Modeling parameter to extend the volume of influence  
REAL(dbl)    	        	:: R_P, Sh_P, delta_P
REAl(dbl)               	:: R_influence_p, L_influence_p			! Parameters related to particle's volume of influence
REAl(dbl)               	:: V_influence_P, V_eff_Ratio			! Parameters related to particle's volume of influence
REAL(dbl)			:: Cb_Total_Veff
REAL(dbl),DIMENSION(2)   	:: VIB_x, VIB_y, VIB_z	 			! Volume of Influence's Borders
REAL(dbl),DIMENSION(2)    	:: NVB_x, NVB_y, NVB_z				! Node Volume's Borders
REAL(dbl)                       :: Delta_X, Delta_Y, Delta_Z
REAL(dbl)                       :: x_DP, y_DP, z_DP				! Coordinates of "Discretized Point" (DP)
REAL(dbl)                       :: Cb_Hybrid

TYPE(ParRecord), POINTER  	:: current
TYPE(ParRecord), POINTER  	:: next

delta_mesh = 1.0_dbl

!zcf3 = xcf*ycf*zcf
zcf3 = 1.0

current => ParListHead%next
DO WHILE (ASSOCIATED(current))

!------ Copy pointer of next node
	next => current%next
        IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!------ Calculate length scale for jth particle:  delta = R / Sh
!------ Calculate effective radius: R_influence_P = R + (N_b *delta)
!------ Note: need to convert this into Lattice units and not use the physical length units
!------ Then compute equivalent cubic mesh length scale
	N_b = 1.0
        R_P  = current%pardata%rp
	Sh_P = current%pardata%sh
        delta_P = R_P / Sh_P
        R_influence_P = (R_P + N_b * delta_P) / xcf

!------ Computing equivalent cubic mesh length scale
        V_influence_P= (4.0_dbl/3.0_dbl) * PI * R_influence_P**3.0_dbl
        L_influence_P= V_influence_P **(1.0_dbl/3.0_dbl)

        V_eff_Ratio = V_influence_P / zcf3 					! Ratio of the effective volume to cell size 

!------ Finding particle location in this processor
	xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
	yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
	zp= current%pardata%zp - REAL(kMin-1_lng,dbl)


!----------------------------------------------------------------------------------------------------------------------
!------ Veff is smaller than the mesh volume --> Cb = Trilinear interpolation of the concentration at particle location
!----------------------------------------------------------------------------------------------------------------------
        IF (V_eff_Ratio .LE. 1.0) THEN 					
           CaseNo = 1
	   ix0 =FLOOR(xp)
           ix1 =CEILING(xp)
           iy0 =FLOOR(yp)
           iy1 =CEILING(yp)
           iz0 =FLOOR(zp)
           iz1 =CEILING(zp)

!----------TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES

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

!--------- Concentration Trilinear Iinterpolation
!--------- Interpolation in x-direction
           c00 = phi(ix0,iy0,iz0) * (1.0_dbl-xd) + phi(ix1,iy0,iz0) * xd
           c01 = phi(ix0,iy0,iz1) * (1.0_dbl-xd) + phi(ix1,iy0,iz1) * xd
           c10 = phi(ix0,iy1,iz0) * (1.0_dbl-xd) + phi(ix1,iy1,iz0) * xd
           c11 = phi(ix0,iy1,iz1) * (1.0_dbl-xd) + phi(ix1,iy1,iz1) * xd
!--------- Interpolation in y-direction
           c0  = c00 * (1.0_dbl-yd) + c10 * yd
           c1  = c01 * (1.0_dbl-yd) + c11 * yd
!--------- Interpolation in z-direction
           c   = c0 * (1.0_dbl-zd) + c1 * zd

           Cb_Hybrid= c 

!----------------------------------------------------------------------------------------------------------------------
!------ Veff is slightly larger than mesh volume --> Volume of influence is discretized 
!------ Cb= Average of concentration interpolated on each of the descritized nodes inside volume of influence 
!----------------------------------------------------------------------------------------------------------------------
 	ELSE IF ( (V_eff_Ratio .GT. 1.0) .AND. (V_eff_Ratio .LT. 27.0 ) ) THEN		
           CaseNo = 2
!--------- NEW: Volume of Influence Border (VIB) for this particle
           VIB_x(1)= xp - 0.5_dbl * L_influence_P
           VIB_x(2)= xp + 0.5_dbl * L_influence_P
           VIB_y(1)= yp - 0.5_dbl * L_influence_P
           VIB_y(2)= yp + 0.5_dbl * L_influence_P
           VIB_z(1)= zp - 0.5_dbl * L_influence_P
           VIB_z(2)= zp + 0.5_dbl * L_influence_P

!--------- Discretizing the volume of influence to  make sure at least 64 points are available
           Delta_X = (VIB_x(2)-VIB_x(1)) / 3.0 
	   Delta_Y = (VIB_y(2)-VIB_y(1)) / 3.0
	   Delta_Z = (VIB_z(2)-VIB_z(1)) / 3.0 

           Cb_Total_Veff = 0
           NumFluids_Veff = 0

!--------- Loop over discretized points and averaging the concentration
           DO i= 0, 3
              DO j= 0, 3
                 DO k= 0, 3
                    x_DP = VIB_x(1) + (i * Delta_X) 
                    y_DP = VIB_y(1) + (j * Delta_Y)
                    z_DP = VIB_z(1) + (k * Delta_Z)

!------------------ Finding Lattice nodes surrounding this point (This point is discretized and is not a lattice node))
                    ix0 = FLOOR(x_DP)
                    ix1 = CEILING(x_DP)
                    iy0 = FLOOR(y_DP)
                    iy1 = CEILING(y_DP)
                    iz0 = FLOOR(z_DP)
                    iz1 = CEILING(z_DP)

!------------------ TO BE DONE: MAKE SURE THE ABOVE NODES ARE FLUID NODES
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

!------------------ Taking care of the periodic BC in Z-dir
                    iz00 = iz0
                    IF (iz0 .gt. nz) THEN
                       iz00 = iz0 - (nz - 1)
                    ELSE IF (iz0 .lt. 1) THEN
                       iz00 = iz0 + (nz-1)
                    END IF

                    iz11 = iz1         
                    IF (iz1 .gt. nz) THEN
                       iz11 = iz1 - (nz-1)
                    ELSE IF (iz1 .lt. 1) THEN
                       iz11 = iz1 + (nz-1)
                    END IF

!------------------ Concentration Trilinear Iinterpolation
!------------------ Interpolation in x-direction
                    c00 = phi(ix0,iy0,iz00) * (1.0_dbl-xd) + phi(ix1,iy0,iz00) * xd
                    c01 = phi(ix0,iy0,iz11) * (1.0_dbl-xd) + phi(ix1,iy0,iz11) * xd
                    c10 = phi(ix0,iy1,iz00) * (1.0_dbl-xd) + phi(ix1,iy1,iz00) * xd
                    c11 = phi(ix0,iy1,iz11) * (1.0_dbl-xd) + phi(ix1,iy1,iz11) * xd
!------------------ Interpolation in y-direction
                    c0  = c00 * (1.0_dbl-yd) + c10 * yd
                    c1  = c01 * (1.0_dbl-yd) + c11 * yd
!------------------ Interpolation in z-direction
                    c   = c0 * (1.0_dbl-zd) + c1 * zd

                    Cb_Total_Veff  = Cb_Total_Veff  + c
                    NumFluids_Veff = NumFluids_Veff + 1_lng
                 END DO
             END DO
          END DO
          
          Cb_Hybrid= Cb_Total_Veff / NumFluids_Veff

!----------------------------------------------------------------------------------------------------------------------
!------ Veff is much larger than mesh volume --> Cb= total number of moles in volume of influence / volume of influence 
!----------------------------------------------------------------------------------------------------------------------
        ELSE IF (V_eff_Ratio .GE. 27.0) THEN                             
           CaseNo = 3
!--------- NEW: Volume of Influence Border (VIB) for this particle
           VIB_x(1)= xp - 0.5_dbl * L_influence_P
           VIB_x(2)= xp + 0.5_dbl * L_influence_P
           VIB_y(1)= yp - 0.5_dbl * L_influence_P
           VIB_y(2)= yp + 0.5_dbl * L_influence_P
           VIB_z(1)= zp - 0.5_dbl * L_influence_P
           VIB_z(2)= zp + 0.5_dbl * L_influence_P

!--------- NEW: Finding the lattice "Nodes Effected by Particle"
           NEP_x(1)= CEILING(VIB_x(1))
           NEP_x(2)= FLOOR  (VIB_x(2))
           NEP_y(1)= CEILING(VIB_y(1))
           NEP_y(2)= FLOOR  (VIB_y(2))
           NEP_z(1)= CEILING(VIB_z(1))
           NEP_z(2)= FLOOR  (VIB_z(2))

           Cb_Total_Veff = 0.0_dbl
           NumFluids_Veff = 0_lng

           DO i= NEP_x(1),NEP_x(2) 
              DO j= NEP_y(1),NEP_y(2)
                 DO k= NEP_z(1),NEP_z(2)

                    !---- Taking care of the periodic BC in Z-dir
                    kk = k
		    IF (k .gt. nz) THEN
                       kk = k - (nz - 1)  
           	    ELSE IF (k .lt. 1) THEN
                       kk = k + (nz-1)
		    END IF   

                    IF (node(i,j,kk) .EQ. FLUID) THEN
                       Cb_Total_Veff  = Cb_Total_Veff  + phi(i,j,kk)
                       NumFluids_Veff = NumFluids_Veff + 1_lng
                    END IF
                 END DO
             END DO
          END DO
          
          Cb_Hybrid= Cb_Total_Veff / NumFluids_Veff

       END IF
 
       current%pardata%bulk_conc = Cb_Hybrid

       END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       current => next

END DO

!===================================================================================================
END SUBROUTINE Compute_Cb
!===================================================================================================







!===================================================================================================
SUBROUTINE Calc_Scalar_Release! Calculate rate of scalar release at every time step  
!===================================================================================================

! Called by Particle_Track (LBM.f90) to get delNBbyCV, update particle radius,
! Sh(t)- sherwood number

IMPLICIT NONE
INTEGER(lng)  :: numFluids,i,j,k
REAL(dbl)     :: deltaR,bulkVolume,temp,cbt,zcf3,bulkconc
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

bulkVolume=xcf*ycf*zcf
zcf3=xcf*ycf*zcf

!------calculate delNBbyCV for each particle in the domain
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next ! copy pointer of next node
   IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	current%pardata%rpold=current%pardata%rp
        bulkconc = current%pardata%bulk_conc
	temp = current%pardata%rpold**2.0_dbl-4.0_dbl*tcf*molarvol*diffm*current%pardata%sh*max((current%pardata%par_conc-bulkconc),0.0_dbl)
 
        IF (temp.GE.0.0_dbl) THEN
		current%pardata%rp=0.5_dbl*(current%pardata%rpold+sqrt(temp))
	ELSE
          temp = 0.0_dbl
          current%pardata%rp=0.5_dbl*(current%pardata%rpold+sqrt(temp))
	END IF

	deltaR=current%pardata%rpold-current%pardata%rp
	current%pardata%delNBbyCV=(4.0_dbl/3.0_dbl)*PI*(current%pardata%rpold*current%pardata%rpold*current%pardata%rpold &
				-current%pardata%rp*current%pardata%rp*current%pardata%rp) &
				/(molarvol*bulkVolume)

	IF (associated(current,ParListHead%next)) THEN
           write(9,*) iter*tcf,current%pardata%parid,current%pardata%rp,current%pardata%Sh,Cb_global*zcf3*Cb_numFluids,current%pardata%delNBbyCV,Cb_global,Cb_numFluids
           CALL FLUSH(9)
	ENDIF

   END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   current => next
ENDDO


IF (myid .EQ. 0) THEN
   	open(799,file='testoutput.dat',position='append')
	write(799,*) iter*tcf,Cb_global*zcf3*Cb_numFluids,Cb_global,Cb_numFluids
        close(799)
ENDIF

!===================================================================================================
END SUBROUTINE Calc_Scalar_Release
!===================================================================================================







!------------------------------------------------
SUBROUTINE Update_Sh! Incporate hierarchical mdoel to Sh(t) to include effect of shear/hydrodynamics and container effect
! Called by Particle_Track (LBM.f90) to get Calc_SCalar_Release, delNBbyCV, update particle radius
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)  :: ix0,iy0,iz0,ix1,iy1,iz1
INTEGER(lng)  :: it,jt,kt,ib,jb,kb
REAL(dbl)     :: temp,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
REAL(dbl)     :: S,Sst,Sh0
REAL(dbl)     :: xp,yp,zp,xd,yd,zd
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next



!calculate delNBbyCV for each particle in the domain
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node
        IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	! Initialize Sh for this particle
	current%pardata%sh=1.0_dbl/(1.0_dbl-current%pardata%gamma_cont)
	! Add container effect
	current%pardata%sh=current%pardata%sh + (current%pardata%gamma_cont/(1.0_dbl-current%pardata%gamma_cont))

	! Add Shear effect from Yanxing's correlations from Wang et al., (2015)
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

	ib = ix0
	jb = iy0
	kb = iz0
	it = ix0+1_lng
	jt = iy0
	kt = iz0

	dwdz = (w(it,jt,kt)-w(ib,jb,kb))
	S = abs(dwdz*vcf/zcf)
	Sst = S*(current%pardata%rp**2.0)/diffm

	current%pardata%S = S
	current%pardata%Sst = Sst

	IF (Sst.LT.5.0_dbl) THEN
		current%pardata%sh=current%pardata%sh+0.296_dbl*(Sst**0.5_dbl)
	ELSE
		Sh0 = exp(0.162_dbl+0.202_dbl*log(Sst)-7.5e-6_dbl*(log(Sst)**5.4_dbl)) 
		current%pardata%sh=current%pardata%sh+Sh0-1.0_dbl
	END IF


        END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	current => next
ENDDO


!------------------------------------------------
END SUBROUTINE Update_Sh
!------------------------------------------------




!===================================================================================================
SUBROUTINE Compute_shear                
!===================================================================================================
! Computing components of strain rate tensor using central difference everywhere except near the
! processor boundaries where sided difference is used. 

IMPLICIT NONE
INTEGER(lng)  :: i,j,k
REAL(dbl)  :: temp

!===================================================================================================
END SUBROUTINE Compute_shear
!===================================================================================================








!===================================================================================================
SUBROUTINE Interp_ParToNodes_Conc_New 
!===================================================================================================

!--- Interpolate Particle concentration release to node locations 
!--- Called by Particle_Track (LBM.f90) to get delphi_particle

IMPLICIT NONE
INTEGER(lng)  		  :: i,j,k,kk
REAL(dbl)     		  :: xp,yp,zp
REAL(dbl)		  :: delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
REAL(dbl)                 :: N_d         				! Modeling parameter to extend the volume of influence around 
REAL(dbl)                 :: R_P, Sh_P, delta_P
REAL(dbl)                 :: R_influence_P, L_influence_P
REAL(dbl),DIMENSION(2)    :: VIB_x, VIB_y, VIB_z 			! Volume of Influence's Borders
REAL(dbl),DIMENSION(2)    :: NVB_x, NVB_y, NVB_z			! Node Volume's Borders
INTEGER  ,DIMENSION(2)    :: LN_x,  LN_y,  LN_z				! Lattice Nodes Surronding the particle
INTEGER  ,DIMENSION(2)    :: NEP_x, NEP_y, NEP_z                        ! Lattice Nodes Surronding the particle
REAL(dbl)		  :: tmp, Overlap_sum
TYPE(ParRecord), POINTER  :: current
TYPE(ParRecord), POINTER  :: next


delta_mesh = 1.0_dbl
zcf3 = xcf*ycf*zcf
current => ParListHead%next

DO WHILE (ASSOCIATED(current))

!------ Copy pointer of next node
	next => current%next 
        IF (mySub .EQ.current%pardata%cur_part) THEN !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!------ Calculate length scale for jth particle:  delta = R / Sh
!------ Calculate effective radius: R_influence_P = R + (N_d *delta)
!------ Note: need to convert this into Lattice units and not use the physical length units
!------ Then compute equivalent cubic mesh length scale
	N_d = 1.0
        R_P  = current%pardata%rp
	Sh_P = current%pardata%sh
        delta_P = R_P / Sh_P
        R_influence_P = (R_P + N_d * delta_P) / xcf

!------ Computing equivalent cubic mesh length scale
        L_influence_P = ( (4.0_dbl*PI/3.0_dbl) * R_influence_P**3.0_dbl)**(1.0_dbl/3.0_dbl)
 
!------ Finding particle location in this processor
	xp= current%pardata%xp - REAL(iMin-1_lng,dbl)
	yp= current%pardata%yp - REAL(jMin-1_lng,dbl)
	zp= current%pardata%zp - REAL(kMin-1_lng,dbl)

!------ NEW: Volume of Influence Border (VIB) for this particle
        VIB_x(1)= xp - 0.5_dbl * L_influence_P
        VIB_x(2)= xp + 0.5_dbl * L_influence_P
        VIB_y(1)= yp - 0.5_dbl * L_influence_P
        VIB_y(2)= yp + 0.5_dbl * L_influence_P
        VIB_z(1)= zp - 0.5_dbl * L_influence_P
        VIB_z(2)= zp + 0.5_dbl * L_influence_P

!------ NEW: Finding the lattice "Nodes Effected by Particle" 
        NEP_x(1)= FLOOR(VIB_x(1))
        NEP_x(2)= CEILING(VIB_x(2))
        NEP_y(1)= FLOOR(VIB_y(1))
        NEP_y(2)= CEILING(VIB_y(2))
        NEP_z(1)= FLOOR(VIB_z(1))
        NEP_z(2)= CEILING(VIB_z(2))

!------ NEW: Finding the volume overlapping between particle-effetive-volume and the volume around each lattice node
        Overlap_sum = 0.0_dbl
        Overlap= 0.0
 
        DO i= NEP_x(1),NEP_x(2) 
           DO j= NEP_y(1),NEP_y(2)
              DO k= NEP_z(1),NEP_z(2)
	         NVB_x(1) = REAL(i,dbl) - 0.5_dbl*delta_mesh
        	 NVB_x(2) = REAL(i,dbl) + 0.5_dbl*delta_mesh
		 NVB_y(1) = REAL(j,dbl) - 0.5_dbl*delta_mesh
		 NVB_y(2) = REAL(j,dbl) + 0.5_dbl*delta_mesh
		 NVB_z(1) = REAL(k,dbl) - 0.5_dbl*delta_mesh
		 NVB_z(2) = REAL(k,dbl) + 0.5_dbl*delta_mesh

		 tmp = MAX ( MIN(VIB_x(2),NVB_x(2)) - MAX(VIB_x(1),NVB_x(1)), 0.0_dbl) * &
                       MAX ( MIN(VIB_y(2),NVB_y(2)) - MAX(VIB_y(1),NVB_y(1)), 0.0_dbl) * &
                       MAX ( MIN(VIB_z(2),NVB_z(2)) - MAX(VIB_z(1),NVB_z(1)), 0.0_dbl)
                
		!---- Taking care of the periodic BC in Z-dir
                kk = k 
      		IF (k .gt. nz) THEN
                    kk = k - (nz - 1)
                ELSE IF (k .lt. 1) THEN
                    kk = k + (nz-1)
                END IF

                IF (node(i,j,kk) .EQ. FLUID) THEN
                    Overlap(i,j,kk)= tmp
                    Overlap_sum= Overlap_sum + Overlap(i,j,kk)
                 END IF
              END DO
           END DO
        END DO

!------ Computing NB_j and Veff for each particle
!        Nbj = 0.0_dbl                                                           ! initialize Nbj - the number of moles of drug in the effective volume surrounding the particle
!        Veff = 0.0_dbl                                                          ! initialize Veff - the eff. volume of each particle
!
!------ Solving an equation for Rj/Reff in order to estimate Veff and Nbj (see notes form July 2015)
!        CALL Find_Root(current%pardata%parid, current%pardata%bulk_conc, current%pardata%par_conc &
!                      ,current%pardata%gamma_cont,current%pardata%rp,Nbj,Veff)
!        current%pardata%Veff = Veff                                             ! store Veff in particle record
!        current%pardata%Nbj = Nbj                                               ! store Nbj in particle record
!        Nbj = Nbj/zcf3  
!
!------ Computing particle release contribution to scalar field at each lattice node
        DO i= NEP_x(1),NEP_x(2)
           DO j= NEP_y(1),NEP_y(2)
              DO k= NEP_z(1),NEP_z(2)
                 
		 !----- Taking care of the periodic BC in Z-dir
                 kk = k 
                 IF (k .gt. nz) THEN
                     kk = k - (nz-1)
                 ELSE IF (k .lt. 1) THEN
                     kk = k + (nz-1)
                 END IF
 
                 IF (node(i,j,kk) .EQ. FLUID) THEN                 
                 
		    !----- Overlap_sum going to zero when the particle is disapearing
		    IF (Overlap_sum .gt. 1e-40) THEN
                       Overlap(i,j,kk) = Overlap(i,j,kk) / Overlap_sum
                    ELSE
                       Overlap(i,j,kk) = 0.0
                    END IF

                    delphi_particle(i,j,kk)  = delphi_particle(i,j,kk)  + current%pardata%delNBbyCV * Overlap(i,j,kk) 

!                   tausgs_particle_x(i,j,k)= tausgs_particle_x(i,j,k)- current%pardata%up*Nbj   * (Overlap(i,j,k)/Overlap_sum)
!                   tausgs_particle_y(i,j,k)= tausgs_particle_y(i,j,k)- current%pardata%up*Nbj   * (Overlap(i,j,k)/Overlap_sum)
!                   tausgs_particle_z(i,j,k)= tausgs_particle_z(i,j,k)- current%pardata%up*Nbj   * (Overlap(i,j,k)/Overlap_sum)
                 END IF 
              END DO
           END DO
        END DO

END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
current => next
ENDDO

!===================================================================================================
END SUBROUTINE Interp_ParToNodes_Conc_New  		 
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
REAL(dbl),parameter :: eps = 1.0e-12_dbl,tol = 1.0e-8_dbl,largenum = 1.0e8_dbl

parcb = conc
parcs = cs

nmax = 100_lng
iter = 0_lng
xnew = 0.0_dbl
xold = 0.5_dbl

! The conc values are quite small. So we will make them larger so that the
! coeffs in the eq for Rj/Reff are not small.
parcb = parcb*largenum
parcs = parcs*largenum


Aj = (parcb-gammaj*parcs)/(1.0_dbl-gammaj)
Bj = (parcs - parcb)/max((parcb-gammaj*cs),eps)
AjBj = (parcs-parcb)/(1.0_dbl-gammaj)

a = Aj + 1.5_dbl*AjBj!*Aj*Bj
b = -1.5_dbl*AjBj!Aj*Bj
c = parcb - Aj
error = abs(xold-xnew)

!write(*,*) 'test code', Nbj,conc,Veff
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
!xnew = max(min(xnew,0.5_dbl),0.01)
xnew = max(min(xnew,1.0_dbl),0.01) ! Limit xnew (radius ratio) to values that are meaningful and not too small or large. 
Reff = Rj/xnew
Veff = (88.0_dbl/21.0_dbl)*(Reff**3.0_dbl)
Nbj = conc*Veff

!===================================================================================================
END SUBROUTINE Find_Root
!===================================================================================================






!===================================================================================================
SUBROUTINE Particle_Track
!===================================================================================================
IMPLICIT NONE

REAL(dbl)      		 :: xpold(1:np),ypold(1:np),zpold(1:np),z_tmp 		! old particle coordinates (working coordinates are stored in xp,yp,zp)
REAL(dbl)      		 :: upold(1:np),vpold(1:np),wpold(1:np) 		! old particle velocity components (new vales are stored in up, vp, wp)
REAL(dbl)                :: Cb_Domain, Cb_Local, Cb_Hybrid, V_eff_Ratio
REAL(dbl)                :: ZZZP
INTEGER(lng)   		 :: i,ipartition,ii,jj,kk, CaseNo
INTEGER(dbl)             :: RANKK
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
   
!-- Particle tracking is done, now time for drug relaes calculations------------------------------------------------------------------------------------------------------------

   CALL Interp_bulkconc(Cb_Local)  					! interpolate final bulk_concentration after the final position is ascertained.
   CALL Calc_Global_Bulk_Scalar_Conc(Cb_Domain)
   CALL Compute_Cb(V_eff_Ratio,CaseNo,Cb_Hybrid)  
   open(172,file='Cb-history.dat',position='append')
   write(172,*) iter, V_eff_Ratio, CaseNo, Cb_Local, Cb_Domain, Cb_Hybrid
   CALL Update_Sh 							! Update the Sherwood number for each particle depending on the shear rate at the particle location. 
   CALL Calc_Scalar_Release 						! Updates particle radius, calculates new drug conc release rate delNBbyCV. 
   CALL Interp_ParToNodes_Conc_New 					! distributes released drug concentration to neightbouring nodes 
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
	
	
	SELECT CASE(current%pardata%parid)
	CASE(1_lng)
      open(72,file='particle1-history.dat',position='append')
      write(72,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up*vcf,current%pardata%vp*vcf,current%pardata%wp*vcf,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(72)
	CASE(3_lng)
      open(73,file='particle3-history.dat',position='append')
      write(73,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(73) 
	CASE(5_lng)
      open(74,file='particle5-history.dat',position='append')
      write(74,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(74)
	CASE(7_lng)
      open(75,file='particle7-history.dat',position='append')
      write(75,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(75)
	CASE(9_lng)
      open(76,file='particle9-history.dat',position='append')
      write(76,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(76) 
	CASE(10_lng)
      open(77,file='particle10-history.dat',position='append')
      write(77,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(77)
	CASE(8_lng)
      open(78,file='particle8-history.dat',position='append')
      write(78,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(78)
	CASE(6_lng)
      open(79,file='particle6-history.dat',position='append')
      write(79,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(79)
	CASE(4_lng)
      open(80,file='particle4-history.dat',position='append')
      write(80,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(80)
	CASE(2_lng)
      open(81,file='particle2-history.dat',position='append')
      write(81,*) iter,iter*tcf,current%pardata%xp,current%pardata%yp,current%pardata%zp,current%pardata%up,current%pardata%vp,current%pardata%wp,current%pardata%sh,current%pardata%rp,current%pardata%bulk_conc,current%pardata%delNBbyCV,current%pardata%cur_part,current%pardata%new_part
      close(81)
     
      END SELECT
   END IF !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   current => next 
ENDDO


!---- Parallel communication between all processors
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
   next => current%next 
	RANKK= current%pardata%cur_part - 1
        current%pardata%cur_part = current%pardata%new_part 

	CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
	CALL MPI_BCast(current%pardata%xp,        1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%yp,        1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%zp,        1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%up,        1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%vp,        1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%wp,        1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%rp,        1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%sh,        1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%xpold,     1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%ypold,     1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%zpold,     1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%upold,     1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%vpold,     1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%wpold,     1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%rpold,     1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%delNBbyCV, 1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%par_conc,  1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%bulk_conc, 1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%gamma_cont,1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%Nbj,       1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%S,         1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%Sst,       1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
	CALL MPI_BCast(current%pardata%cur_part,  1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        CALL MPI_BCast(current%pardata%new_part,  1, MPI_DOUBLE_PRECISION, RANKK, MPI_COMM_WORLD,mpierr)
        
        write(*,*) iter,'CPU, cur_P, Rank, z', mySub, current%pardata%cur_part, RANKK, current%pardata%zp
   current => next  
ENDDO


!===================================================================================================
END SUBROUTINE Particle_Track
!===================================================================================================






!--------------------------------------------------------------------------------------------------
SUBROUTINE Equilibrium		! calculate the equilibrium distribution function and set f to feq (initial condition)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng)	:: i,j,k,m					! index variables
REAL(dbl)		:: uu,ue,ve,we,Usum		! precalculated quantities for use in the feq equation
REAL(dbl)		:: feq						! equilibrium distribution function

! Balaji modified to change indices form 0 to nzSub+1
DO k=1,nzSub+0
  DO j=1,nySub+0
    DO i=1,nxSub+0

      IF(node(i,j,k) .EQ. FLUID) THEN
      
        uu = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)						! u . u
      
        DO m=0,NumDistDirs
        
          ue	= u(i,j,k)*ex(m)																			! u . e
          ve	= v(i,j,k)*ey(m)																			! v . e
          we	= w(i,j,k)*ez(m)																			! w . e

          Usum	= ue + ve + we																				! U . e
        
          feq = (wt(m)*rho(i,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum - 1.5_dbl*uu)	! equilibrium distribution function

          f(m,i,j,k) = feq    

        END DO

      END IF
      
    END DO
  END DO
END DO

!------------------------------------------------
END SUBROUTINE Equilibrium
!------------------------------------------------



!--------------------------------------------------------------------------------------------------
SUBROUTINE Collision		! calculates equilibrium distribution function AND collision step for each node
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng)	:: i,j,k,m				! index variables
REAL(dbl)	:: UU,ue,ve,we,Usum			! precalculated quantities for use in the feq equation
REAL(dbl)	:: feq					! equilibrium distribution function

! Balaji modified to change indices form 0 to nzSub+1
DO k=1,nzSub+0
  DO j=1,nySub+0
    DO i=1,nxSub+0

      IF(node(i,j,k) .EQ. FLUID) THEN

        UU = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)						! U . U
      
        DO m=0,NumDistDirs
        
          ue	= u(i,j,k)*ex(m)										! u . e
          ve	= v(i,j,k)*ey(m)										! v . e
          we	= w(i,j,k)*ez(m)										! w . e

          Usum	= ue + ve + we											! U . e
        
          feq	= (wt(m)*rho(i,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5*Usum*Usum - 1.5*uu)	! equilibrium distribution function
          f(m,i,j,k)		= f(m,i,j,k) - oneOVERtau*(f(m,i,j,k) - feq)					! collision
	!write(*,*) 0.5_dbl*s1/vcf,s1,vcf,xcf,nu,nuL,tau!m,i,j,k,f(m,i,j,k),w(i,j,k)
        
        END DO 

      END IF
      
    END DO
  END DO
END DO

!------------------------------------------------
END SUBROUTINE Collision
!------------------------------------------------







!--------------------------------------------------------------------------------------------------
SUBROUTINE Stream	! stream the distribution functions between neighboring nodes (stream - using Lallemand 2nd order moving BB)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1		! index variables
REAL(dbl) :: fbb								! bounced back distribution function

fplus = f									! store the post-collision distribution function

! interior nodes (away from other subdomains)
DO k=2,nzSub-1
  DO j=2,nySub-1
    DO i=2,nxSub-1

      IF(node(i,j,k) .EQ. FLUID) THEN
      
        DO m=1,NumDistDirs
        
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)
    
          IF(node(im1,jm1,km1) .EQ. FLUID) THEN 
            f(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN					! macro- boundary
	    ! Balaji added after commenting out the earlier method
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)					  					! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi
            CALL BounceBackV2(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)							! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM.f90 at Line 1570: node(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter,'mySub',mySub
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
            WRITE(1000,*) "node(i,j,k)=",node(i,j,k)
            WRITE(1000,*) "node(im1,jm1,km1)=",node(im1,jm1,km1)
            CLOSE(1000)
            STOP
          END IF

        END DO    

      END IF
 
    END DO
  END DO
END DO

! XY faces
DO k=1,nzSub,(nzSub-1)
  DO j=1,nySub
    DO i=1,nxSub

      IF(node(i,j,k) .EQ. FLUID) THEN
      
        DO m=1,NumDistDirs
        
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)

          !IF(km1.eq.0) km1=nzSub ! Balaji added

          IF(node(im1,jm1,km1) .EQ. FLUID) THEN 
            f(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN				! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  						! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN				! villi
            CALL BounceBackVL(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)						! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM.f90 at Line 1618: node(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
            WRITE(1000,*) "node(i,j,k)=",node(i,j,k)
            WRITE(1000,*) "node(im1,jm1,km1)=",node(im1,jm1,km1)
            CLOSE(1000)
            STOP
          END IF

        END DO    

      END IF

    END DO
  END DO
END DO

! XZ faces
DO j=1,nySub,(nySub-1)
  DO k=1,nzSub
    DO i=1,nxSub

      IF(node(i,j,k) .EQ. FLUID) THEN
      
        DO m=1,NumDistDirs
        
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)

        
          IF(node(im1,jm1,km1) .EQ. FLUID) THEN
            f(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN		! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  				! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi
            CALL BounceBackVL(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)				! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM.f90 at Line :1664 node(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
            WRITE(1000,*) "node(i,j,k)=",node(i,j,k)
            WRITE(1000,*) "node(im1,jm1,km1)=",node(im1,jm1,km1)
            CLOSE(1000)
            STOP
          END IF

        END DO    

      END IF

    END DO
  END DO
END DO

! YZ faces
DO i=1,nxSub,(nxSub-1)
  DO k=1,nzSub
    DO j=1,nySub

      IF(node(i,j,k) .EQ. FLUID) THEN
      
        DO m=1,NumDistDirs
        
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)
        
          IF(node(im1,jm1,km1) .EQ. FLUID) THEN 
            f(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN		! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  				! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi
            CALL BounceBackVL(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)				! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
		    WRITE(1000,'(A75)') "error in LBM.f90 at Line 1709: node(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter", iter
            WRITE(1000,*) "m=",m
            WRITE(1000,*) "i=",i,"j=",j,"k=",k
            WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
            WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
            WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
            WRITE(1000,*) "node(i,j,k)=",node(i,j,k)
            WRITE(1000,*) "node(im1,jm1,km1)=",node(im1,jm1,km1)
            CLOSE(1000)
            STOP
          END IF

        END DO    

      END IF

    END DO
  END DO
END DO

!------------------------------------------------
END SUBROUTINE Stream
!------------------------------------------------





!--------------------------------------------------------------------------------------------------
SUBROUTINE Macro	! calculate the macroscopic quantities
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m						! index variables
INTEGER(lng) :: ii,jj,kk

LOGICAL :: nodebounflag

! Balaji modified to include 0 to nzSub+1
DO k=1,nzSub
  DO j=1,nySub
    DO i=1,nxSub
      
      IF(node(i,j,k) .EQ. FLUID) THEN

        ! initialize arrays
        rho(i,j,k)		= 0.0_dbl				! density
        u(i,j,k)		= 0.0_dbl				! x-velocity
        v(i,j,k)		= 0.0_dbl				! y-velocity
        w(i,j,k)		= 0.0_dbl				! z-velocity     

        DO m=0,NumDistDirs  
          rho(i,j,k)	= rho(i,j,k) + f(m,i,j,k)			! density
          u(i,j,k)	= u(i,j,k)   + f(m,i,j,k)*ex(m)			! x-velocity
          v(i,j,k)	= v(i,j,k)   + f(m,i,j,k)*ey(m)			! y-velocity
          w(i,j,k)	= w(i,j,k)   + f(m,i,j,k)*ez(m)			! z-velocity
        END DO

        IF(rho(i,j,k) .NE. 0) THEN
          u(i,j,k) = u(i,j,k)/rho(i,j,k)				! x-velocity
          v(i,j,k) = v(i,j,k)/rho(i,j,k)				! y-velocity
          w(i,j,k) = w(i,j,k)/rho(i,j,k)				! z-velocity
        ELSE          

          OPEN(6678,FILE='error.'//sub//'.txt')
          WRITE(6678,*) 'rho(i,j,k) = 0: Line 362 in Macro in LBM.f90'
          WRITE(6678,*) 'iter', iter
          WRITE(6678,*) 'i,j,k:', i,j,k
          WRITE(6678,*) 'node(i,j,k)', node(i,j,k)
          WRITE(6678,*) 'rho(i,j,k)', rho(i,j,k)
          WRITE(6678,*)
          WRITE(6678,*)
          DO m=1,NumDistDirs
            ii = i + ex(m)
            jj = j + ey(m)
            kk = k + ez(m)       
            WRITE(6678,*) 'ii,jj,kk:', ii,jj,kk
            WRITE(6678,*) 'node(ii,jj,kk)', node(ii,jj,kk)
            WRITE(6678,*) 'rho(ii,jj,kk)', rho(ii,jj,kk)
            WRITE(6678,*)
          END DO
          CLOSE(6678)

          OPEN(1001,FILE='rhoMacro.'//sub//'.dat')
          WRITE(1001,'(A60)') 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
          WRITE(1001,'(8E15.5,I6)') x(i), y(j), z(k), u(i,j,k), v(i,j,k), w(i,j,k), (rho(i,j,k)-denL)*dcf*pcf, phi(i,j,k), node(i,j,k)
          CLOSE(1001)

          CALL PrintFieldsTEST									! output the velocity, density, and scalar fields [MODULE: Output]

          STOP

        END IF

      ELSE
        rho(i,j,k)	= denL									! density (zero gauge pressure)
        u(i,j,k)		= 0.0_dbl							! x-velocity
        v(i,j,k)		= 0.0_dbl							! y-velocity
        w(i,j,k)		= 0.0_dbl							! z-velocity
        phi(i,j,k)	= phiWall								! scalar
      END IF

    END DO
  END DO
END DO   

!------------------------------------------------
END SUBROUTINE Macro
!------------------------------------------------



!================================================
END MODULE LBM
!================================================
