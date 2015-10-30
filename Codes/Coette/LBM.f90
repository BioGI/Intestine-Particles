!==================================================================================================
MODULE LBM				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================

USE SetPrecision
USE Setup
USE ICBC

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
SUBROUTINE Interp_Parvel_1 ! Using a crude interpolation approach
!===================================================================================================

IMPLICIT NONE
INTEGER(lng)  :: i
!REAL(dbl)     :: s1,s2,s3,s4,x1,y1,a,b,c,d
REAL(dbl)     :: xp,yp,zp
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node
	xp = current%pardata%xp - REAL(iMin-1_lng,dbl)
	yp = current%pardata%yp - REAL(jMin-1_lng,dbl)
	zp = current%pardata%zp - REAL(kMin-1_lng,dbl)
	current%pardata%up=0.5*(u(FLOOR(xp),FLOOR(yp),FLOOR(zp))+u(CEILING(xp),CEILING(yp),CEILING(zp)))
	current%pardata%vp=0.5*(v(FLOOR(xp),FLOOR(yp),FLOOR(zp))+v(CEILING(xp),CEILING(yp),CEILING(zp)))
	current%pardata%wp=0.5*(w(FLOOR(xp),FLOOR(yp),FLOOR(zp))+w(CEILING(xp),CEILING(yp),CEILING(zp)))
	! point to next node in the list
	current => next
	!write(*,*) i
ENDDO

!===================================================================================================
END SUBROUTINE Interp_Parvel_1 ! Using a crude interpolation approach
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

! point to next node in the list
	current => next

ENDDO

!===================================================================================================
END SUBROUTINE Interp_Parvel ! Using Trilinear interpolation
!===================================================================================================





!===================================================================================================
SUBROUTINE Interp_Parvel_Crude ! Using a crde interpolation approach
!===================================================================================================

IMPLICIT NONE
INTEGER(lng)  :: i
REAL(dbl)     :: s1,s2,s3,s4,x1,y1,a,b,c,d
REAL(dbl)     :: xp,yp,zp
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node
	xp = current%pardata%xp - REAL(iMin-1_lng,dbl)
	yp = current%pardata%yp - REAL(jMin-1_lng,dbl)
	zp = current%pardata%zp - REAL(kMin-1_lng,dbl)
	current%pardata%up=0.5*(u(FLOOR(xp),FLOOR(yp),FLOOR(zp))+u(CEILING(xp),CEILING(yp),CEILING(zp)))
	current%pardata%vp=0.5*(v(FLOOR(xp),FLOOR(yp),FLOOR(zp))+v(CEILING(xp),CEILING(yp),CEILING(zp)))
	current%pardata%wp=0.5*(w(FLOOR(xp),FLOOR(yp),FLOOR(zp))+w(CEILING(xp),CEILING(yp),CEILING(zp)))
	! point to next node in the list
	current => next
	!write(*,*) i
ENDDO


!===================================================================================================
END SUBROUTINE Interp_Parvel_Crude ! Using a crde interpolation approach
!===================================================================================================







!===================================================================================================
SUBROUTINE Interp_bulkconc ! Using Trilinear interpolation
!===================================================================================================

IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd
REAL(dbl)     :: xp,yp,zp
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node

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
        current%pardata%bulk_conc=c

! point to next node in the list
	current => next
ENDDO

!===================================================================================================
END SUBROUTINE Interp_bulkconc ! Using Trilinear interpolation
!===================================================================================================







!===================================================================================================
SUBROUTINE Calc_Global_Bulk_Scalar_Conc !Calculate Global Bulk SCalar COnc for use in the scalar drug relese model 
!===================================================================================================

IMPLICIT NONE
INTEGER(lng)  :: i,j,k

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

!===================================================================================================
END SUBROUTINE Calc_Global_Bulk_Scalar_Conc
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


!Cb_global = 0.0
!bulkVolume=xcf*ycf*zcf*Cb_numFluids/num_particles
bulkVolume=xcf*ycf*zcf
zcf3=xcf*ycf*zcf

!calculate delNBbyCV for each particle in the domain
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node

	current%pardata%rpold=current%pardata%rp

	bulkconc = Cb_global
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

! point to next node in the list
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

	!IF (associated(current,ParListHead%next)) THEN
	IF (current%pardata%parid.eq.1) THEN
		write(*,*) iter*tcf,S,Sst,current%pardata%sh,current%pardata%cur_part,dwdz,w(it,jt,kt),w(ib,jb,kb),ib,it,jt,kt,node(ib,jb,kb),node(it,jt,kt),zp,FLOOR(zp),CEILING(zp),w(it,jt,0),w(it,jt,nzSub+1_lng),nzSub
	ENDIF

	! point to next node in the list
	current => next
ENDDO


!------------------------------------------------
END SUBROUTINE Update_Sh
!------------------------------------------------






!------------------------------------------------
SUBROUTINE Compute_shear
!------------------------------------------------
! The goal of this is subroutine is to compute a shear field
IMPLICIT NONE
INTEGER(lng)  :: i,j,k
REAL(dbl)  :: temp

! Notes: We will compute the components of the strain rate tensor using central difference everywhere except near the processor boundaries where we will use one sided difference. 

!------------------------------------------------
END SUBROUTINE Compute_shear
!------------------------------------------------







!------------------------------------------------
SUBROUTINE Interp_ParToNodes_Conc ! Interpolate Particle concentration release to node locations 
! Called by Particle_Track (LBM.f90) to get delphi_particle
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
REAL(dbl)  :: ax0,ax1,ay0,ay1,az0,az1
REAL(dbl)  :: bx0,bx1,by0,by1,bz0,bz1
REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd
REAL(dbl)     :: c000,c010,c100,c110,c001,c011,c101,c111,csum
REAL(dbl)     :: xp,yp,zp,delta_par,delta_mesh,zcf3,Nbj,Veff,bulkconc
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next


delta_mesh = 1.0_dbl
zcf3=xcf*ycf*zcf

current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node

	!write(*,*) current%pardata%parid,iter,mySub,current%pardata%new_part,xp,iMin,yp,jMin,zp,kMin

	xp = current%pardata%xp - REAL(iMin-1_lng,dbl)
	yp = current%pardata%yp - REAL(jMin-1_lng,dbl)
	zp = current%pardata%zp - REAL(kMin-1_lng,dbl)

        delta_par = ((current%pardata%rp/xcf)/current%pardata%sh)+(current%pardata%rp/xcf) ! calculate length scale delta_j = R_j/Sh_j for jth particle and effective radius delta_par = delta_j + R_j (note: need to convert this into Lattice units and not use the physical length units)
        delta_par = ((88.0_dbl/21.0_dbl)*((delta_par)**3.0_dbl))**(1.0_dbl/3.0_dbl) ! compute equivalent cubic mesh length scale
        !delta_par = delta_mesh

	ix0=FLOOR(xp)
	ix1=CEILING(xp)
	iy0=FLOOR(yp)
	iy1=CEILING(yp)
	iz0=FLOOR(zp)
	iz1=CEILING(zp)

! Boundaries of the volume of influence of the particle
	ax0 = xp - 0.5_dbl*delta_par
	ax1 = xp + 0.5_dbl*delta_par
	ay0 = yp - 0.5_dbl*delta_par
	ay1 = yp + 0.5_dbl*delta_par
	az0 = zp - 0.5_dbl*delta_par
	az1 = zp + 0.5_dbl*delta_par

	if (iz0.GT.nzSub) then
		iz0 = iz0-1_lng
	elseif (iz0.LT.1_lng) then
		iz0 = iz0 +1_lng
	endif
	if (ix0.GT.nxSub) then
		ix0 = ix0-1_lng
	elseif (ix0.LT.1) then
		ix0 = ix0 +1_lng
	endif
	if (iy0.GT.nySub) then
		iy0 = iy0-1_lng
	elseif (iy0.LT.1) then
		iy0 = iy0 +1_lng
	endif
	if (iz1.GT.nzSub) then
		iz1 = iz1-1_lng
	elseif (iz1.LT.1_lng) then
		iz1 = iz1 +1_lng
	endif
	if (ix1.GT.nxSub) then
		ix1 = ix1-1_lng
	elseif (ix1.LT.1) then
		ix1 = ix1 +1_lng
	endif
	if (iy1.GT.nySub) then
		iy1 = iy1-1_lng
	elseif (iy1.LT.1) then
		iy1 = iy1 +1_lng
	endif

	bx0 = REAL(ix0,dbl) - 0.5_dbl*delta_mesh
	bx1 = REAL(ix0,dbl) + 0.5_dbl*delta_mesh
	by0 = REAL(iy0,dbl) - 0.5_dbl*delta_mesh
	by1 = REAL(iy0,dbl) + 0.5_dbl*delta_mesh
	bz0 = REAL(iz0,dbl) - 0.5_dbl*delta_mesh
	bz1 = REAL(iz0,dbl) + 0.5_dbl*delta_mesh
	c000 = MAX(MIN(ax1,bx1)-MAX(ax0,bx0),0.0_dbl)*MAX(MIN(ay1,by1)-MAX(ay0,by0),0.0_dbl)!*MAX(MIN(az1,bz1)-MAX(az0,bz0),0.0_dbl)
	bx0 = REAL(ix0,dbl) - 0.5_dbl*delta_mesh
	bx1 = REAL(ix0,dbl) + 0.5_dbl*delta_mesh
	by0 = REAL(iy1,dbl) - 0.5_dbl*delta_mesh
	by1 = REAL(iy1,dbl) + 0.5_dbl*delta_mesh
	bz0 = REAL(iz0,dbl) - 0.5_dbl*delta_mesh
	bz1 = REAL(iz0,dbl) + 0.5_dbl*delta_mesh
	c010 = MAX(MIN(ax1,bx1)-MAX(ax0,bx0),0.0_dbl)*MAX(MIN(ay1,by1)-MAX(ay0,by0),0.0_dbl)!*MAX(MIN(az1,bz1)-MAX(az0,bz0),0.0_dbl)
	bx0 = REAL(ix1,dbl) - 0.5_dbl*delta_mesh
	bx1 = REAL(ix1,dbl) + 0.5_dbl*delta_mesh
	by0 = REAL(iy0,dbl) - 0.5_dbl*delta_mesh
	by1 = REAL(iy0,dbl) + 0.5_dbl*delta_mesh
	bz0 = REAL(iz0,dbl) - 0.5_dbl*delta_mesh
	bz1 = REAL(iz0,dbl) + 0.5_dbl*delta_mesh
	c100 = MAX(MIN(ax1,bx1)-MAX(ax0,bx0),0.0_dbl)*MAX(MIN(ay1,by1)-MAX(ay0,by0),0.0_dbl)!*MAX(MIN(az1,bz1)-MAX(az0,bz0),0.0_dbl)
	bx0 = REAL(ix1,dbl) - 0.5_dbl*delta_mesh
	bx1 = REAL(ix1,dbl) + 0.5_dbl*delta_mesh
	by0 = REAL(iy1,dbl) - 0.5_dbl*delta_mesh
	by1 = REAL(iy1,dbl) + 0.5_dbl*delta_mesh
	bz0 = REAL(iz0,dbl) - 0.5_dbl*delta_mesh
	bz1 = REAL(iz0,dbl) + 0.5_dbl*delta_mesh
	c110 = MAX(MIN(ax1,bx1)-MAX(ax0,bx0),0.0_dbl)*MAX(MIN(ay1,by1)-MAX(ay0,by0),0.0_dbl)*MAX(MIN(az1,bz1)-MAX(az0,bz0),0.0_dbl)

	bx0 = REAL(ix0,dbl) - 0.5_dbl*delta_mesh
	bx1 = REAL(ix0,dbl) + 0.5_dbl*delta_mesh
	by0 = REAL(iy0,dbl) - 0.5_dbl*delta_mesh
	by1 = REAL(iy0,dbl) + 0.5_dbl*delta_mesh
	bz0 = REAL(iz1,dbl) - 0.5_dbl*delta_mesh
	bz1 = REAL(iz1,dbl) + 0.5_dbl*delta_mesh
	c001 = MAX(MIN(ax1,bx1)-MAX(ax0,bx0),0.0_dbl)*MAX(MIN(ay1,by1)-MAX(ay0,by0),0.0_dbl)*MAX(MIN(az1,bz1)-MAX(az0,bz0),0.0_dbl)
	bx0 = REAL(ix0,dbl) - 0.5_dbl*delta_mesh
	bx1 = REAL(ix0,dbl) + 0.5_dbl*delta_mesh
	by0 = REAL(iy1,dbl) - 0.5_dbl*delta_mesh
	by1 = REAL(iy1,dbl) + 0.5_dbl*delta_mesh
	bz0 = REAL(iz1,dbl) - 0.5_dbl*delta_mesh
	bz1 = REAL(iz1,dbl) + 0.5_dbl*delta_mesh
	c011 = MAX(MIN(ax1,bx1)-MAX(ax0,bx0),0.0_dbl)*MAX(MIN(ay1,by1)-MAX(ay0,by0),0.0_dbl)*MAX(MIN(az1,bz1)-MAX(az0,bz0),0.0_dbl)
	bx0 = REAL(ix1,dbl) - 0.5_dbl*delta_mesh
	bx1 = REAL(ix1,dbl) + 0.5_dbl*delta_mesh
	by0 = REAL(iy0,dbl) - 0.5_dbl*delta_mesh
	by1 = REAL(iy0,dbl) + 0.5_dbl*delta_mesh
	bz0 = REAL(iz1,dbl) - 0.5_dbl*delta_mesh
	bz1 = REAL(iz1,dbl) + 0.5_dbl*delta_mesh
	c101 = MAX(MIN(ax1,bx1)-MAX(ax0,bx0),0.0_dbl)*MAX(MIN(ay1,by1)-MAX(ay0,by0),0.0_dbl)*MAX(MIN(az1,bz1)-MAX(az0,bz0),0.0_dbl)
	bx0 = REAL(ix1,dbl) - 0.5_dbl*delta_mesh
	bx1 = REAL(ix1,dbl) + 0.5_dbl*delta_mesh
	by0 = REAL(iy1,dbl) - 0.5_dbl*delta_mesh
	by1 = REAL(iy1,dbl) + 0.5_dbl*delta_mesh
	bz0 = REAL(iz1,dbl) - 0.5_dbl*delta_mesh
	bz1 = REAL(iz1,dbl) + 0.5_dbl*delta_mesh
	c111 = MAX(MIN(ax1,bx1)-MAX(ax0,bx0),0.0_dbl)*MAX(MIN(ay1,by1)-MAX(ay0,by0),0.0_dbl)*MAX(MIN(az1,bz1)-MAX(az0,bz0),0.0_dbl)

        csum = c000 + c010 + c100 + c110 + c001 + c011 + c101 + c111

	Nbj = 0.0_dbl ! initialize Nbj - the number of moles of drug in the effective volume surrounding the particle
	Veff = 0.0_dbl ! initialize Veff - the eff. volume of each particle
	bulkconc = Cb_global
	!bulkconc = current%pardata%bulk_conc
	! Need to solve an equation for Rj/Reff in order to estimate Veff and Nbj (see notes form July 2015)
	CALL Find_Root(current%pardata%parid,bulkconc,current%pardata%par_conc &
		      ,current%pardata%gamma_cont,current%pardata%rp,Nbj,Veff)
	current%pardata%Veff = Veff ! store Veff in particle record
	current%pardata%Nbj = Nbj	! store Nbj in particle record
	Nbj = Nbj/zcf3 ! convert Nbj (number of moles) to a conc by division with cell volume (dimensional) 

        IF (csum.EQ.0.0) THEN ! csum = 0 then dump the drug molecules to the nearest fluid node
	        delphi_particle(ix0,iy0,iz0)   = delphi_particle(ix0,iy0,iz0)+current%pardata%delNBbyCV!(phi(ix0,iy0,iz0)/bulk_conc(i))
		tausgs_particle_x(ix0,iy0,iz0) = tausgs_particle_x(ix0,iy0,iz0) - current%pardata%up*Nbj*1.0_dbl
		tausgs_particle_y(ix0,iy0,iz0) = tausgs_particle_y(ix0,iy0,iz0) - current%pardata%vp*Nbj*1.0_dbl
		tausgs_particle_z(ix0,iy0,iz0) = tausgs_particle_z(ix0,iy0,iz0) - current%pardata%wp*Nbj*1.0_dbl

        ELSE ! if not, then distribute it according to the model. This helps us to conserve the total number of drug molecules

                c000 = c000/csum
        	c010 = c010/csum
        	c100 = c100/csum
        	c110 = c110/csum
        	c001 = c001/csum
        	c011 = c011/csum
        	c101 = c101/csum
        	c111 = c111/csum

        	!delphi_particle(ix0,iy0,iz0)=delphi_particle(ix0,iy0,iz0)+current%pardata%delNBbyCV!(phi(ix0,iy0,iz0)/bulk_conc(i))
                
        	delphi_particle(ix0,iy0,iz0)=delphi_particle(ix0,iy0,iz0)+current%pardata%delNBbyCV*c000
        	delphi_particle(ix0,iy1,iz0)=delphi_particle(ix0,iy1,iz0)+current%pardata%delNBbyCV*c010
        	delphi_particle(ix1,iy0,iz0)=delphi_particle(ix1,iy0,iz0)+current%pardata%delNBbyCV*c100
        	delphi_particle(ix1,iy1,iz0)=delphi_particle(ix1,iy1,iz0)+current%pardata%delNBbyCV*c110
        	delphi_particle(ix0,iy0,iz1)=delphi_particle(ix0,iy0,iz1)+current%pardata%delNBbyCV*c001
        	delphi_particle(ix0,iy1,iz1)=delphi_particle(ix0,iy1,iz1)+current%pardata%delNBbyCV*c011
        	delphi_particle(ix1,iy0,iz1)=delphi_particle(ix1,iy0,iz1)+current%pardata%delNBbyCV*c101
        	delphi_particle(ix1,iy1,iz1)=delphi_particle(ix1,iy1,iz1)+current%pardata%delNBbyCV*c111

        	tausgs_particle_x(ix0,iy0,iz0)=tausgs_particle_x(ix0,iy0,iz0)- current%pardata%up*Nbj*c000
		tausgs_particle_x(ix0,iy1,iz0)=tausgs_particle_x(ix0,iy1,iz0)- current%pardata%up*Nbj*c010
		tausgs_particle_x(ix1,iy0,iz0)=tausgs_particle_x(ix1,iy0,iz0)- current%pardata%up*Nbj*c100
        	tausgs_particle_x(ix1,iy1,iz0)=tausgs_particle_x(ix1,iy1,iz0)- current%pardata%up*Nbj*c110
        	tausgs_particle_x(ix0,iy0,iz1)=tausgs_particle_x(ix0,iy0,iz1)- current%pardata%up*Nbj*c001
        	tausgs_particle_x(ix0,iy1,iz1)=tausgs_particle_x(ix0,iy1,iz1)- current%pardata%up*Nbj*c011
        	tausgs_particle_x(ix1,iy0,iz1)=tausgs_particle_x(ix1,iy0,iz1)- current%pardata%up*Nbj*c101
        	tausgs_particle_x(ix1,iy1,iz1)=tausgs_particle_x(ix1,iy1,iz1)- current%pardata%up*Nbj*c111

        	tausgs_particle_y(ix0,iy0,iz0)=tausgs_particle_y(ix0,iy0,iz0)- current%pardata%vp*Nbj*c000
        	tausgs_particle_y(ix0,iy1,iz0)=tausgs_particle_y(ix0,iy1,iz0)- current%pardata%vp*Nbj*c010
        	tausgs_particle_y(ix1,iy0,iz0)=tausgs_particle_y(ix1,iy0,iz0)- current%pardata%vp*Nbj*c100
        	tausgs_particle_y(ix1,iy1,iz0)=tausgs_particle_y(ix1,iy1,iz0)- current%pardata%vp*Nbj*c110
        	tausgs_particle_y(ix0,iy0,iz1)=tausgs_particle_y(ix0,iy0,iz1)- current%pardata%vp*Nbj*c001
        	tausgs_particle_y(ix0,iy1,iz1)=tausgs_particle_y(ix0,iy1,iz1)- current%pardata%vp*Nbj*c011
        	tausgs_particle_y(ix1,iy0,iz1)=tausgs_particle_y(ix1,iy0,iz1)- current%pardata%vp*Nbj*c101
        	tausgs_particle_y(ix1,iy1,iz1)=tausgs_particle_y(ix1,iy1,iz1)- current%pardata%vp*Nbj*c111

        	tausgs_particle_z(ix0,iy0,iz0)=tausgs_particle_z(ix0,iy0,iz0)- current%pardata%wp*Nbj*c000
        	tausgs_particle_z(ix0,iy1,iz0)=tausgs_particle_z(ix0,iy1,iz0)- current%pardata%wp*Nbj*c010
        	tausgs_particle_z(ix1,iy0,iz0)=tausgs_particle_z(ix1,iy0,iz0)- current%pardata%wp*Nbj*c100
        	tausgs_particle_z(ix1,iy1,iz0)=tausgs_particle_z(ix1,iy1,iz0)- current%pardata%wp*Nbj*c110
        	tausgs_particle_z(ix0,iy0,iz1)=tausgs_particle_z(ix0,iy0,iz1)- current%pardata%wp*Nbj*c001
        	tausgs_particle_z(ix0,iy1,iz1)=tausgs_particle_z(ix0,iy1,iz1)- current%pardata%wp*Nbj*c011
        	tausgs_particle_z(ix1,iy0,iz1)=tausgs_particle_z(ix1,iy0,iz1)- current%pardata%wp*Nbj*c101
        	tausgs_particle_z(ix1,iy1,iz1)=tausgs_particle_z(ix1,iy1,iz1)- current%pardata%wp*Nbj*c111

        END IF


	! point to next node in the list
	current => next
ENDDO

!------------------------------------------------
END SUBROUTINE Interp_ParToNodes_Conc ! Interpolate Particle concentration release to node locations 
!------------------------------------------------





!------------------------------------------------
SUBROUTINE Find_Root(parid,conc,cs,gammaj,Rj,Nbj,Veff)
!------------------------------------------------
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
!write(*,*) 'test code', iter,parid,Nbj,conc,Veff,xnew,Rj,iter,error
!write(*,*) parid,iter,error,Rj,Reff,Nbj
!------------------------------------------------
END SUBROUTINE Find_Root
!------------------------------------------------





!------------------------------------------------
SUBROUTINE Particle_Track
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)   :: i,ipartition,ii,jj,kk
REAL(dbl)      :: xpold(1:np),ypold(1:np),zpold(1:np) ! old particle coordinates (working coordinates are stored in xp,yp,zp)
!REAL(dbl)      :: xp(1:np),yp(1:np),zp(1:np) 	      ! working particle coordinates (working coordinates are stored in xp,yp,zp)
!REAL(dbl)      :: xpnew(1:np),ypnew(1:np),zpnew(1:np) ! new particle coordinates (working coordinates are stored in xp,yp,zp)
REAL(dbl)      :: upold(1:np),vpold(1:np),wpold(1:np) ! old particle velocity components (new vales are stored in up, vp, wp)
!REAL(dbl)      :: up(1:np),vp(1:np),wp(1:np) 	      ! working particle velocity (working coordinates are stored in xp,yp,zp)
TYPE(ParRecord), POINTER :: current
TYPE(ParRecord), POINTER :: next

ParticleTransfer = .FALSE. ! AT this time we do not know if any particles need to be transferred.
delphi_particle = 0.0_dbl ! set delphi_particle to 0.0 before every time step, when the particle drug release happens. 
tausgs_particle_x = 0.0_dbl
tausgs_particle_y = 0.0_dbl
tausgs_particle_z = 0.0_dbl
	
IF (iter.GT.iter0+0_lng) THEN ! IF condition ensures that at the first step, the only part of this subroutine that operates is computing the partitions the particles belong to without releasing any drug.  
! Second order interpolation in time
!Backup particle data from previous time step
! Using a linked list of particle records
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node

	current%pardata%xpold = current%pardata%xp
	current%pardata%ypold = current%pardata%yp
	current%pardata%zpold = current%pardata%zp
	
	current%pardata%upold = current%pardata%up
	current%pardata%vpold = current%pardata%vp
	current%pardata%wpold = current%pardata%wp
	
	current%pardata%xp=current%pardata%xpold+current%pardata%up
	current%pardata%yp=current%pardata%ypold+current%pardata%vp
	current%pardata%zp=current%pardata%zpold+current%pardata%wp
	
	!IF(current%pardata%zp.GE.REAL(nz,dbl)) THEN
	!	current%pardata%zp = MOD(current%pardata%zp,REAL(nz,dbl))
	!ENDIF

	!yp(i)=Cj ! test

	! point to next node in the list
	current => next
	!write(*,*) i
ENDDO
CALL Interp_Parvel
! Using a linked list of particle records
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node

	current%pardata%xp=current%pardata%xpold+0.5*(current%pardata%up+current%pardata%upold)
	current%pardata%yp=current%pardata%ypold+0.5*(current%pardata%vp+current%pardata%vpold)
	current%pardata%zp=current%pardata%zpold+0.5*(current%pardata%wp+current%pardata%wpold)
	!xpnew(i)=current%pardata%xp
	!ypnew(i)=current%pardata%yp
	!zpnew(i)=current%pardata%zp

	!write(*,*) current%parid
	! point to next node in the list
	current => next
	!write(*,*) i
ENDDO

!IF (ASSOCIATED(ParListHead%next)) THEN
!	write(*,*) 'In Particle_Track','iter= ',iter,'SUBID =',mySub,'ParID =',ParListHead%next%pardata%parid,'wp =',ParListHead%next%pardata%wp*vcf,'Transfer FLAG =',ParticleTransfer,ParListHead%next%pardata%cur_part,ParListHead%next%pardata%new_part
!END IF

CALL Interp_Parvel ! interpolate final particle velocities after the final position is ascertained. 
CALL Interp_bulkconc ! interpolate final bulk_concentration after the final position is ascertained.
CALL Update_Sh ! Update the Sherwood number for each particle depending on the shear rate at the particle location. 
!delphi_particle = 0.0_dbl ! set delphi_particle to 0.0 before every time step, when the particle drug release happens. 
CALL Calc_Scalar_Release ! Updates particle radius, calculates new drug conc release rate delNBbyCV. 
CALL Interp_ParToNodes_Conc ! distributes released drug concentration to neightbouring nodes 
!drug molecules released by the particle at this new position
ENDIF


! Now update tausgs only for those cells that have non-zero values of tausgs
DO kk=0,nzSub+1
        DO jj=0,nySub+1
                DO ii=0,nxSub+1
			if (tausgs_particle_x(ii,jj,kk).ne.0.0_dbl) then
                        	tausgs_particle_x(ii,jj,kk) = u(ii,jj,kk)*phi(ii,jj,kk)
			endif
			if (tausgs_particle_y(ii,jj,kk).ne.0.0_dbl) then
                        	tausgs_particle_y(ii,jj,kk) = v(ii,jj,kk)*phi(ii,jj,kk)
			endif
			if (tausgs_particle_z(ii,jj,kk).ne.0.0_dbl) then
                        	tausgs_particle_z(ii,jj,kk) = w(ii,jj,kk)*phi(ii,jj,kk)
			endif
                ENDDO
        ENDDO
ENDDO



!IF (myid .EQ. master) THEN
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node
	
	! Wrappign around in z-direction for periodic BC in z
	IF(current%pardata%zp.GE.REAL(nz,dbl)) THEN
		current%pardata%zp = MOD(current%pardata%zp,REAL(nz,dbl))
	ENDIF
	IF(current%pardata%zp.LE.0.0_dbl) THEN
		current%pardata%zp = current%pardata%zp+REAL(nz,dbl)
	ENDIF

	! Wrappign around in y-direction for periodic BC in y
	IF(current%pardata%yp.GE.REAL(ny,dbl)) THEN
		current%pardata%yp = MOD(current%pardata%yp,REAL(ny,dbl))
	ENDIF
	IF(current%pardata%yp.LT.1.0_dbl) THEN
		current%pardata%yp = current%pardata%yp+REAL(ny,dbl)
	ENDIF


	! Estimate to which partition the updated position belongs to.
	!ParticleTransfer = .FALSE.
	DO ipartition = 1_lng,NumSubsTotal 


		IF ((current%pardata%xp.GE.REAL(iMinDomain(ipartition),dbl)-1.0_dbl).AND.&
		(current%pardata%xp.LT.(REAL(iMaxDomain(ipartition),dbl)+0.0_dbl)).AND. &
		(current%pardata%yp.GE.REAL(jMinDomain(ipartition),dbl)-1.0_dbl).AND. &
		(current%pardata%yp.LT.(REAL(jMaxDomain(ipartition),dbl)+0.0_dbl)).AND. &
		(current%pardata%zp.GE.REAL(kMinDomain(ipartition),dbl)-1.0_dbl).AND. &
		(current%pardata%zp.LT.(REAL(kMaxDomain(ipartition),dbl)+0.0_dbl))) THEN

			current%pardata%new_part = ipartition
		END IF
                !write(*,*) ipartition,kMinDomain(ipartition),kMaxDomain(ipartition)
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
	! point to next node in the list
	current => next
	!write(*,*) i
ENDDO
!ENDIF 

!------------------------------------------------
END SUBROUTINE Particle_Track
!------------------------------------------------






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

INTEGER(lng)	:: i,j,k,m					! index variables
REAL(dbl)		:: UU,ue,ve,we,Usum		! precalculated quantities for use in the feq equation
REAL(dbl)		:: feq						! equilibrium distribution function

! Balaji modified to change indices form 0 to nzSub+1
DO k=1,nzSub+0
  DO j=1,nySub+0
    DO i=1,nxSub+0

      IF(node(i,j,k) .EQ. FLUID) THEN

        UU = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)						! U . U
      
        DO m=0,NumDistDirs
        
          ue	= u(i,j,k)*ex(m)																			! u . e
          ve	= v(i,j,k)*ey(m)																			! v . e
          we	= w(i,j,k)*ez(m)																			! w . e

          Usum	= ue + ve + we																				! U . e
        
          feq	= (wt(m)*rho(i,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5*Usum*Usum - 1.5*uu)	! equilibrium distribution function
          f(m,i,j,k)		= f(m,i,j,k) - oneOVERtau*(f(m,i,j,k) - feq)							! collision
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
