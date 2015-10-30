!==================================================================================================
MODULE LBM				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================
USE SetPrecision
USE Setup
USE ICBC

IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE LBM_Setup	! setup the LBM simulation
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Initialize variables and arrays
f		= 0.0_dbl		! distribution functions
fplus	= 0.0_dbl		! post-collision distribution functions
u  	= 0.0_dbl		! x-velocity
v     = 0.0_dbl		! y-velocity
w     = 0.0_dbl		! z-velocity
rho   = 0.0_lng		! density
ex    = 0.0_dbl		! lattice discretized velocity vector (x-component)
ey    = 0.0_dbl		! lattice discretized velocity vector (x-component)
ez    = 0.0_dbl		! lattice discretized velocity vector (x-component)
bb    = 0_lng			! bounceback directions
wt    = 0.0_dbl		! weighting coefficients

! Fill out weighting coefficient array
wt(0)    = 2.0_dbl/9.0_dbl	
wt(1:6)  = 1.0_dbl/9.0_dbl
wt(7:14) = 1.0_dbl/72.0_dbl

! Fill out bounceback array
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

! Fill out symmetry array
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

! Fill velocity direction vector arrays
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

! Define other simulation parameters
nuL   		= (2.0_dbl*tau - 1.0_dbl)/6.0_dbl	! lattice kinematic viscosity
denL 			= 1.0_dbl									! arbitrary lattice density (1.0 for convenience)
oneOVERtau 	= 1.0_dbl/tau								! reciprical of tau
cs				= (1.0_dbl)/(SQRT(3.0_dbl))			! speed of sound on the lattice

! Initialize timestep
iter = 0_lng												! intialize the starting timestep to 0 - will get reset in 'ICs' in ICBCM.f90

! Calculate feq for initial condition
CALL Equilibrium

!------------------------------------------------
END SUBROUTINE LBM_Setup
!------------------------------------------------

!------------------------------------------------
SUBROUTINE Particle_Setup
!------------------------------------------------
IMPLICIT NONE
IF (restart) THEN
ELSE
	CALL Interp_Parvel
ENDIF
!------------------------------------------------
END SUBROUTINE Particle_Setup
!------------------------------------------------
!------------------------------------------------
SUBROUTINE Interp_Parvel_1 ! Using a crude interpolation approach
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)  :: i
REAL(dbl)     :: s1,s2,s3,s4,x1,y1,a,b,c,d

DO i=1,np
	up(i)=0.5*(u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+u(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	vp(i)=0.5*(v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+v(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	wp(i)=0.5*(w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+w(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	!up(i)=u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
	!vp(i)=v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
	!wp(i)=w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
ENDDO
!------------------------------------------------
END SUBROUTINE Interp_Parvel_1 ! Using a crude interpolation approach
!------------------------------------------------
!------------------------------------------------
SUBROUTINE Interp_Parvel ! Using Trilinear interpolation
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd

DO i=1,np
	!ix0=FLOOR(xp(i))
	!ix1=CEILING(xp(i))
	!iy0=FLOOR(yp(i))
	!iy1=CEILING(yp(i))
	!iz0=FLOOR(zp(i))
	!iz1=CEILING(zp(i))

	!yp(i)=Cj ! test

	ix0=FLOOR(xp(i))
	ix1=FLOOR(xp(i))+1_lng
	iy0=FLOOR(yp(i))
	iy1=FLOOR(yp(i))+1_lng
	iz0=FLOOR(zp(i))
	iz1=FLOOR(zp(i))+1_lng
	xd=(xp(i)-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
	yd=(yp(i)-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
	zd=(zp(i)-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
	!yd=0.0_dbl ! TEST: used to keep particle motion plainly 2-D ! Balaji added

!	! New interpolation to ensure symmetry
!	ix0=FLOOR(xp(i))
!	IF (xp(i).GT.Ci) THEN
!		ix1=FLOOR(xp(i))+1_lng
!	ELSE
!		ix1=FLOOR(xp(i))-1_lng
!	ENDIF
!	iy0=FLOOR(yp(i))
!	IF (yp(i).GT.Cj) THEN ! It is important to interpolate on the radial direction to ensure symmetry. 
!		iy1=FLOOR(yp(i))+1_lng
!	ELSE
!		iy1=FLOOR(yp(i))-1_lng
!	ENDIF
!	iz0=FLOOR(zp(i))
!	iz1=FLOOR(zp(i))+1_lng
!	xd=(xp(i)-REAL(ix0,dbl))/ABS(REAL(ix1,dbl)-REAL(ix0,dbl))	
!	yd=(yp(i)-REAL(iy0,dbl))/ABS(REAL(iy1,dbl)-REAL(iy0,dbl))	
!	zd=(zp(i)-REAL(iz0,dbl))/ABS(REAL(iz1,dbl)-REAL(iz0,dbl))


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
        up(i)=c

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
        vp(i)=c

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
        wp(i)=c

	!up(i)=0.5*(u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+u(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	!vp(i)=0.5*(v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+v(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	!wp(i)=0.5*(w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+w(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	!up(i)=u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
	!vp(i)=v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
	!wp(i)=w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
ENDDO
!------------------------------------------------
END SUBROUTINE Interp_Parvel ! Using Trilinear interpolation
!------------------------------------------------
!------------------------------------------------
SUBROUTINE Interp_Parvel_Crude ! Using a crde interpolation approach
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)  :: i
REAL(dbl)     :: s1,s2,s3,s4,x1,y1,a,b,c,d

DO i=1,np
	up(i)=0.5*(u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+u(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	vp(i)=0.5*(v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+v(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	wp(i)=0.5*(w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+w(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
	!up(i)=u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
	!vp(i)=v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
	!wp(i)=w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
ENDDO
!------------------------------------------------
END SUBROUTINE Interp_Parvel_Crude ! Using a crde interpolation approach
!------------------------------------------------
!------------------------------------------------
SUBROUTINE Interp_bulkconc ! Using Trilinear interpolation
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd

DO i=1,np
	!ix0=FLOOR(xp(i))
	!ix1=CEILING(xp(i))
	!iy0=FLOOR(yp(i))
	!iy1=CEILING(yp(i))
	!iz0=FLOOR(zp(i))
	!iz1=CEILING(zp(i))

	!yp(i)=Cj ! test

	ix0=FLOOR(xp(i))
	ix1=FLOOR(xp(i))+1_lng
	iy0=FLOOR(yp(i))
	iy1=FLOOR(yp(i))+1_lng
	iz0=FLOOR(zp(i))
	iz1=FLOOR(zp(i))+1_lng
	xd=(xp(i)-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
	yd=(yp(i)-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
	zd=(zp(i)-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
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
        bulk_conc(i)=c
        !bulk_conc(i)=phi(ix0,iy0,iz0)

ENDDO
!------------------------------------------------
END SUBROUTINE Interp_bulkconc ! Using Trilinear interpolation
!------------------------------------------------
!------------------------------------------------
SUBROUTINE Calc_Scalar_Release! Interpolate Particle concentration release to node locations 
! Called by Particle_Track (LBM.f90) to get delNBbyCV, update particle radius,
! Sh(t)- sherwood number
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd
REAL(dbl)     :: deltaR,bulkVolume


DO i=1,np
	rpold(i)=rp(i)
	!if (i.eq.1) then
	!write(*,*) rp(i),tcf,molarvol,diffm,sh(i),par_conc(i)
	!pause
	!endif
	sh(i)=1.0_dbl/(1.0_dbl-gamma_cont(i))
	rp(i)=0.5_dbl*(rpold(i)+sqrt(max(rpold(i)*rpold(i)-4.0_dbl*tcf*molarvol*diffm*sh(i)*(par_conc(i)-bulk_conc(i)),0.0_dbl)))
	deltaR=rpold(i)-rp(i)
	bulkVolume=xcf*ycf*zcf
	delNBbyCV(i)=4.0_dbl*PI*rp(i)*rp(i)*deltaR*par_conc(i)/(bulkVolume*(1.0_dbl-gamma_cont(i)))

ENDDO
!------------------------------------------------
END SUBROUTINE Calc_Scalar_Release
!------------------------------------------------

!------------------------------------------------
SUBROUTINE Interp_ParToNodes_Conc ! Interpolate Particle concentration release to node locations 
! Called by Particle_Track (LBM.f90) to get delphi_particle
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd

DO i=1,np
	ix0=FLOOR(xp(i))
	ix1=FLOOR(xp(i))+1_lng
	iy0=FLOOR(yp(i))
	iy1=FLOOR(yp(i))+1_lng
	iz0=FLOOR(zp(i))
	iz1=FLOOR(zp(i))+1_lng
	!delphi_particle(ix0,iy0,iz0)=delNBbyCV(i)*0.125!(phi(ix0,iy0,iz0)/bulk_conc(i))
	!delphi_particle(ix0,iy1,iz0)=delNBbyCV(i)*0.125!(phi(ix0,iy1,iz0)/bulk_conc(i))
	!delphi_particle(ix1,iy0,iz0)=delNBbyCV(i)*0.125!(phi(ix1,iy0,iz0)/bulk_conc(i))
	!delphi_particle(ix1,iy1,iz0)=delNBbyCV(i)*0.125!(phi(ix1,iy1,iz0)/bulk_conc(i))
	!delphi_particle(ix0,iy0,iz1)=delNBbyCV(i)*0.125!(phi(ix0,iy0,iz1)/bulk_conc(i))
	!delphi_particle(ix0,iy1,iz1)=delNBbyCV(i)*0.125!(phi(ix0,iy1,iz1)/bulk_conc(i))
	!delphi_particle(ix1,iy0,iz1)=delNBbyCV(i)*0.125!(phi(ix1,iy0,iz1)/bulk_conc(i))
	!delphi_particle(ix1,iy1,iz1)=delNBbyCV(i)*0.125!(phi(ix1,iy1,iz1)/bulk_conc(i))
	delphi_particle(ix0,iy0,iz0)=delNBbyCV(i)*(phi(ix0,iy0,iz0))/(bulk_conc(i))
	delphi_particle(ix0,iy1,iz0)=delNBbyCV(i)*(phi(ix0,iy1,iz0))/(bulk_conc(i))
	delphi_particle(ix1,iy0,iz0)=delNBbyCV(i)*(phi(ix1,iy0,iz0))/(bulk_conc(i))
	delphi_particle(ix1,iy1,iz0)=delNBbyCV(i)*(phi(ix1,iy1,iz0))/(bulk_conc(i))
	delphi_particle(ix0,iy0,iz1)=delNBbyCV(i)*(phi(ix0,iy0,iz1))/(bulk_conc(i))
	delphi_particle(ix0,iy1,iz1)=delNBbyCV(i)*(phi(ix0,iy1,iz1))/(bulk_conc(i))
	delphi_particle(ix1,iy0,iz1)=delNBbyCV(i)*(phi(ix1,iy0,iz1))/(bulk_conc(i))
	delphi_particle(ix1,iy1,iz1)=delNBbyCV(i)*(phi(ix1,iy1,iz1))/(bulk_conc(i))

ENDDO
!------------------------------------------------
END SUBROUTINE Interp_ParToNodes_Conc ! Interpolate Particle concentration release to node locations 
!------------------------------------------------
!------------------------------------------------
SUBROUTINE Particle_Track
!------------------------------------------------
IMPLICIT NONE
INTEGER(lng)   :: i
REAL(dbl)      :: xpold(1:np),ypold(1:np),zpold(1:np) ! old particle coordinates (working coordinates are stored in xp,yp,zp)
REAL(dbl)      :: xpnew(1:np),ypnew(1:np),zpnew(1:np) ! new particle coordinates (working coordinates are stored in xp,yp,zp)
REAL(dbl)      :: upold(1:np),vpold(1:np),wpold(1:np) ! old particle velocity components (new vales are stored in up, vp, wp)

! Second order interpolation in time

!Backup particle data from previous time step
xpold(1:np) = xp(1:np)
ypold(1:np) = yp(1:np)
zpold(1:np) = zp(1:np)
upold(1:np) = up(1:np)
vpold(1:np) = vp(1:np)
wpold(1:np) = wp(1:np)

DO i=1,np
	xp(i)=xpold(i)+up(i)
	yp(i)=ypold(i)+vp(i)
	zp(i)=zpold(i)+wp(i)
	IF(zp(i).GE.REAL(nz,dbl)) THEN
		zp(i) = MOD(zp(i),REAL(nz,dbl))
	ENDIF
	!yp(i)=Cj ! test
ENDDO
CALL Interp_Parvel
!CALL Interp_bulkconc
DO i=1,np
	xp(i)=xpold(i)+0.5*(up(i)+upold(i))
	yp(i)=ypold(i)+0.5*(vp(i)+vpold(i))
	zp(i)=zpold(i)+0.5*(wp(i)+wpold(i))
	xpnew(i)=xp(i)
	ypnew(i)=yp(i)
	zpnew(i)=zp(i)
	IF(zp(i).GE.REAL(nz,dbl)) THEN
		zp(i) = MOD(zp(i),REAL(nz,dbl))
	ENDIF
	!yp(i)=Cj ! test

	!xp(i)=xpold(i)
	!yp(i)=ypold(i)
	!zp(i)=zpold(i)
ENDDO


CALL Interp_Parvel ! interpolate final particle velocities after the final position is ascertained. 
CALL Interp_bulkconc ! interpolate final bulk_concentration after the final position is ascertained.
CALL Calc_Scalar_Release ! Updates particle radius, calculates new drug conc release rate delNBbyCV. 
CALL Interp_ParToNodes_Conc ! At this position, also calculate the number of
!drug molecules released by the particle at this new position 
DO i=1,np
	xp(i)=xpnew(i)
	yp(i)=ypnew(i)
	zp(i)=zpnew(i)

	!yp(i)=Cj ! test
ENDDO

! First Order Interpolation in time
!DO i=1,np
!	xp(i)=xp(i)+up(i)
!	yp(i)=yp(i)+vp(i)
!	zp(i)=zp(i)+wp(i)
!	!CALL INTERP_PARVEL
!	!WRITE(*,*) "Particle tracking ",i,xp(i),yp(i),zp(i)
!ENDDO
!CALL Interp_Parvel


IF (myid .EQ. master) THEN
      open(72,file='particle1-history.dat',position='append')
      write(72,*) iter,xp(1),yp(1),zp(1),up(1),vp(1),wp(1),sh(1),rp(1),bulk_conc(1),delNBbyCV(1)
      close(72)
      open(73,file='particle3-history.dat',position='append')
      write(73,*) iter, xp(3),yp(3),zp(3),up(3),vp(3),wp(3),sh(3),rp(3),bulk_conc(3),delNBbyCV(3)
      close(73) 
      open(74,file='particle5-history.dat',position='append')
      write(74,*) iter, xp(5),yp(5),zp(5),up(5),vp(5),wp(5),sh(5),rp(5),bulk_conc(5),delNBbyCV(5)
      close(74)
      open(75,file='particle7-history.dat',position='append')
      write(75,*) iter, xp(7),yp(7),zp(7),up(7),vp(7),wp(7),sh(7),rp(7),bulk_conc(7),delNBbyCV(7)
      close(75)
      open(76,file='particle9-history.dat',position='append')
      write(76,*) iter, xp(9),yp(9),zp(9),up(9),vp(9),wp(9),sh(9),rp(9),bulk_conc(9),delNBbyCV(9)
      close(76) 
      open(77,file='particle10-history.dat',position='append')
      write(77,*) iter, xp(10),yp(10),zp(10),up(10),vp(10),wp(10),sh(10),rp(10),bulk_conc(10),delNBbyCV(10)
      close(77)
      open(78,file='particle8-history.dat',position='append')
      write(78,*) iter, xp(8),yp(8),zp(8),up(8),vp(8),wp(8),sh(8),rp(8),bulk_conc(8),delNBbyCV(8)
      close(78)
      open(79,file='particle6-history.dat',position='append')
      write(79,*) iter, xp(6),yp(6),zp(6),up(6),vp(6),wp(6),sh(6),rp(6),bulk_conc(6),delNBbyCV(6)
      close(79)
      open(80,file='particle4-history.dat',position='append')
      write(80,*) iter, xp(4),yp(4),zp(4),up(4),vp(4),wp(4),sh(4),rp(4),bulk_conc(4),delNBbyCV(4)
      close(80)
      open(81,file='particle2-history.dat',position='append')
      write(81,*) iter, xp(2),yp(2),zp(2),up(2),vp(2),wp(2),sh(2),rp(2),bulk_conc(2),delNBbyCV(2)
      close(81)
ENDIF

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
!DO k=0,nzSub+1
!  DO j=0,nySub+1
!    DO i=0,nxSub+1

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
!DO k=0,nzSub+1
!  DO j=0,nySub+1
!    DO i=0,nxSub+1

      IF(node(i,j,k) .EQ. FLUID) THEN

        UU = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)						! U . U
      
        DO m=0,NumDistDirs
        
          ue	= u(i,j,k)*ex(m)																			! u . e
          ve	= v(i,j,k)*ey(m)																			! v . e
          we	= w(i,j,k)*ez(m)																			! w . e

          Usum	= ue + ve + we																				! U . e
        
          feq	= (wt(m)*rho(i,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5*Usum*Usum - 1.5*uu)	! equilibrium distribution function
          f(m,i,j,k)		= f(m,i,j,k) - oneOVERtau*(f(m,i,j,k) - feq)							! collision
        
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

fplus = f										! store the post-collision distribution function

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
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN															! macro- boundary
            !CALL BounceBack2(m,i,j,k,im1,jm1,km1,fbb)					  										! implement the bounceback BCs [MODULE: ICBC]
	    ! Balaji added after commenting out the earlier method
            CALL BounceBack2New(m,i,j,k,im1,jm1,km1,fbb)					  										! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi
            CALL BounceBackV2(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)							! implement the bounceback BCs [MODULE: ICBC]
            !CALL BounceBackVL(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)							! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM.f90 at Line 89: node(im1,jm1,km1) is out of range"
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
          !IF(km1.eq.nzSub+1) km1=1 ! Balaji added

          IF(node(im1,jm1,km1) .EQ. FLUID) THEN 
            f(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN															! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  												! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi
            CALL BounceBackVL(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)							! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node(im1,jm1,km1) is out of range"
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
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN									! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  						! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi
            CALL BounceBackVL(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)							! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node(im1,jm1,km1) is out of range"
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
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN									! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  						! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi
            CALL BounceBackVL(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),fbb)							! implement the bounceback BCs [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node(im1,jm1,km1) is out of range"
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

! Balaji modified to include 0 to nzSub+1
DO k=1,nzSub
  DO j=1,nySub
    DO i=1,nxSub
!DO k=0,nzSub+1
!  DO j=0,nySub+1
!    DO i=0,nxSub+1
      
      IF(node(i,j,k) .EQ. FLUID) THEN

        ! initialize arrays
        rho(i,j,k)		= 0.0_dbl								! density
        u(i,j,k)		= 0.0_dbl								! x-velocity
        v(i,j,k)		= 0.0_dbl								! y-velocity
        w(i,j,k)		= 0.0_dbl								! z-velocity     

        DO m=0,NumDistDirs  
          rho(i,j,k)	= rho(i,j,k) + f(m,i,j,k)			! density
          u(i,j,k)	= u(i,j,k)   + f(m,i,j,k)*ex(m)	! x-velocity
          v(i,j,k)	= v(i,j,k)   + f(m,i,j,k)*ey(m)	! y-velocity
          w(i,j,k)	= w(i,j,k)   + f(m,i,j,k)*ez(m)	! z-velocity
        END DO

        IF(rho(i,j,k) .NE. 0) THEN
          u(i,j,k) = u(i,j,k)/rho(i,j,k)					! x-velocity
          v(i,j,k) = v(i,j,k)/rho(i,j,k)					! y-velocity
          w(i,j,k) = w(i,j,k)/rho(i,j,k)					! z-velocity
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

          CALL PrintFieldsTEST										! output the velocity, density, and scalar fields [MODULE: Output]

          STOP

        END IF

      ELSE
      
        rho(i,j,k)	= denL									! density (zero gauge pressure)
        u(i,j,k)		= 0.0_dbl								! x-velocity
        v(i,j,k)		= 0.0_dbl								! y-velocity
        w(i,j,k)		= 0.0_dbl								! z-velocity
        phi(i,j,k)	= phiWall								! scalar
        
	

      END IF

    END DO
  END DO
END DO   

!! Balaji added for serial job ! need to do this for parallel using
!MPI_Transfer
!	!write(*,*) iter,MAXVAL(f(:,:,:,0)-f(:,:,:,nzSub))
!	!write(*,*) iter,MAXVAL(f(:,:,:,nzSub+1)-f(:,:,:,1))
!       u(0:nxSub+1,0:nySub+1,0)=u(0:nxSub+1,0:nySub+1,nzSub)
!       u(0:nxSub+1,0:nySub+1,nzSub+1)=u(0:nxSub+1,0:nySub+1,1)
!       u(0,0:nySub+1,0:nzSub+1)=0.0_dbl
!       u(nxSub+1,0:nySub+1,0:nzSub+1)=0.0_dbl
!       u(:,0,:)=0.0_dbl
!       u(:,nySub+1,:)=0.0_dbl
!	!write(*,*) iter,MAXVAL(f(:,:,:,0)-f(:,:,:,nzSub))
!	!write(*,*) iter,MAXVAL(f(:,:,:,nzSub+1)-f(:,:,:,1))
!
!       v(0:nxSub+1,0:nySub+1,0)=v(0:nxSub,0:nySub+1,nzSub)
!       v(0:nxSub+1,0:nySub+1,nzSub+1)=v(0:nxSub+1,0:nySub+1,1)
!       v(0,0:nySub+1,0:nzSub+1)=0.0_dbl
!       v(nxSub+1,0:nySub+1,0:nzSub+1)=0.0_dbl
!       v(:,0,:)=0.0_dbl
!       v(:,nySub+1,:)=0.0_dbl
!       
!       w(0:nxSub+1,0:nySub+1,0)=w(0:nxSub+1,0:nySub+1,nzSub)
!       w(0:nxSub+1,0:nySub+1,nzSub+1)=w(0:nxSub+1,0:nySub+1,1)
!       w(0,0:nySub+1,0:nzSub+1)=0.0_dbl
!       w(nxSub+1,0:nySub+1,0:nzSub+1)=0.0_dbl
!       w(:,0,:)=0.0_dbl
!       w(:,nySub+1,:)=0.0_dbl
!
!       rho(0:nxSub+1,0:nySub+1,0)=rho(0:nxSub+1,0:nySub+1,nzSub)
!       rho(0:nxSub+1,0:nySub+1,nzSub+1)=rho(0:nxSub+1,0:nySub+1,1)
!       !rho(0,0:nySub+1,0:nzSub+1)=0.0_dbl
!       !rho(nxSub+1,0:nySub+1,0:nzSub+1)=0.0_dbl
!       !rho(:,0,:)=0.0_dbl
!       !rho(:,nySub+1,:)=0.0_dbl


!------------------------------------------------
END SUBROUTINE Macro
!------------------------------------------------

!================================================
END MODULE LBM
!================================================
