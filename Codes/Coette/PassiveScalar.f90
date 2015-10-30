!==================================================================================================
MODULE PassiveScalar		! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs)
!==================================================================================================
USE SetPrecision
USE Setup
USE ICBC

IMPLICIT NONE 

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar_Setup	! sets up the passive scalar component
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! initialize arrays
phi    		= 0.0_dbl								! scalar
phiTemp		= 0.0_dbl								! temporary scalar

! scalar parameters
Dm   			= nuL/Sc									! binary molecular diffusivity (scalar in fluid)
Dmcf			= (zcf*zcf)/tcf						! conversion factor for diffusivity
!Delta			= 1.0_dbl - 6.0_dbl*Dm				! scalar diffusion parameter
Delta			= 1.0_dbl - 6.0_dbl*Dm				! scalar diffusion parameter

! set scalar values at sources/sinks
phiWall		= 0.0_dbl								! value of scalar at the boundary

! set the scalar standard devation for gaussian distributions
sigma = 0.1_dbl*D										! 1/10 of the Diameter

! determine scalar starting iteration
phiStart 	= NINT((phiPer*Tmix)/tcf)
IF (phiPer.EQ.0.0) THEN
	phiStart 	= NINT((phiPer*Tmix)/tcf)+1 ! Balaji changed this to add 1 as for phiPer=0, phiSTart=0. But iter never has a value 0.
ENDIF

phiInNodes = 0.0_dbl
phiOutNodes = 0.0_dbl

!------------------------------------------------
END SUBROUTINE Scalar_Setup
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar								! calculates the evolution of scalar in the domain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1		! index variables
REAL(dbl) :: phiBC				! scalar contribution from boundary
REAL(dbl) :: phiOutSurf,phiInSurf		! scalar contribution coming from and going into the boundary
REAL(dbl) :: tausgs				! contribution form tau_sgs term from particle closure

CALL ScalarDistribution						! sets/maintains initial distributions of scalar [MODULE: ICBC.f90]

! store the previous scalar values
phiTemp = phi

! Stream the scalar
DO k=1,nzSub
  DO j=1,nySub
    DO i=1,nxSub
      
      IF(node(i,j,k) .EQ. FLUID) THEN
      	!IF (iter.lt.1200) THEN
	!	phiTemp(i,j,k) = phiTemp(i,j,k)+ delphi_particle(i,j,k) ! Balaji added to introduce drug concentration release
	!END IF
	!phiTemp(i,j,k) = phiTemp(i,j,k)+ delphi_particle(i,j,k) ! Balaji added to introduce drug concentration release
        phi(i,j,k) = Delta*phiTemp(i,j,k) ! 
        !phi(i,j,k) = 0.0_dbl*Delta*phiTemp(i,j,k) ! 
	phi(i,j,k) = phi(i,j,k)+ delphi_particle(i,j,k) ! Balaji added to introduce drug concentration release
	tausgs = ((tausgs_particle_x(i+1,j,k)-tausgs_particle_x(i-1,j,k)) + &
		 (tausgs_particle_y(i,j+1,k)-tausgs_particle_y(i,j-1,k)) + &
		 (tausgs_particle_z(i,j,k+1)-tausgs_particle_z(i,j,k-1)))*0.5_dbl
	phi(i,j,k) = phi(i,j,k)+ tausgs ! Balaji added - to handle SGS particle effects.

        DO m=0,NumDistDirs
      
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)

          IF(node(im1,jm1,km1) .EQ. FLUID) THEN 
            phi(i,j,k) = phi(i,j,k) + (fplus(m,im1,jm1,km1)/rho(im1,jm1,km1) - wt(m)*Delta)*phiTemp(im1,jm1,km1)
            !phi(i,j,k) = phi(i,j,k) + (fplus(m,im1,jm1,km1)/rho(im1,jm1,km1) - wt(m)*Delta)*phiTemp(im1,jm1,km1) &
	    !								     + wt(m)*Delta*phiTemp(i,j,k)
          ELSE IF((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN		! macro- boundary
            !CALL ScalarBC(m,i,j,k,im1,jm1,km1,phiBC)! Gino's BC															! implement scalar boundary condition (using BB f's)	[MODULE: ICBC]
            !phi(i,j,k) = phi(i,j,k) + phiBC     
            !CALL AbsorbedScalarS(i,j,k,m,phiBC)																	! measure the absorption rate

            CALL ScalarBC2(m,i,j,k,im1,jm1,km1,phiBC,phiOutSurf,phiInSurf)	! Wang's BC														! implement scalar boundary condition (using BB f's)	[MODULE: ICBC]
	    !phiBC =(fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)
            phi(i,j,k) = phi(i,j,k) + phiBC     
            CALL AbsorbedScalarS2(i,j,k,m,phiOutSurf,phiInSurf)																	! measure the absorption rate
            !CALL AbsorbedScalarS(i,j,k,m,phiBC)																	! measure the absorption rate
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi
            CALL ScalarBCV(m,i,j,k,im1,jm1,km1,(-node(im1,jm1,km1)),phiBC)								! implement scalar boundary condition (using BB f's)	[MODULE: ICBC]
            phi(i,j,k) = phi(i,j,k) + phiBC     
            CALL AbsorbedScalarV(i,j,k,m,phiBC)																	! measure the absorption rate
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in PassiveScalar.f90 at Line 89: node(im1,jm1,km1) is out of range"
            WRITE(1000,*) "iter",iter
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

	!phi(i,j,k) = phi(i,j,k) + delphi_particle(i,j,k) ! Balaji added to introduce drug concentration release
       	!  fix spurious oscillations in moment propagation method for high Sc #s
        IF(phi(i,j,k) .LT. 0.0_dbl) THEN
          write(*,*) 'phi is negative'
          STOP
          phi(i,j,k) = 0.0_dbl
        END IF
      END IF

    END DO
  END DO
END DO

! Add the amount of scalar absorbed through the outer and villous surfaces
phiAbsorbed = 	phiAbsorbedS + phiAbsorbedV																		! total amount of scalar absorbed up to current time
      
!------------------------------------------------
END SUBROUTINE Scalar
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarS(i,j,k,m,phiBC)		! measures the total absorbed scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m				! index variables
REAL(dbl), INTENT(IN) :: phiBC     				! scalar contribution from the boundary condition
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

phiIN 	= phiBC																						! contribution from the wall to the crrent node (in)
phiOUT	= (fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)		! contribution to the wall from the current node (out)

phiAbsorbedS = phiAbsorbedS + (phiOUT - phiIN)	!- wt(m)*Delta*phiWall			! add the amount of scalar that has been absorbed at the current location in the current direction

!------------------------------------------------
END SUBROUTINE AbsorbedScalarS
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarS2(i,j,k,m,phiOutSurf,phiInSurf)		! measures the total absorbed scalar ! New method by Balaji, Dec 2014
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m				! index variables
REAL(dbl), INTENT(IN) :: phiOutSurf,phiInSurf     			! scalar contribution from the boundary condition
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

phiIN 	= phiInSurf		! contribution from the wall to the crrent node (in)
phiOUT	= phiOutSurf		! contribution to the wall from the current node (out)

phiAbsorbedS = phiAbsorbedS + (phiOUT - phiIN) !- wt(m)*Delta*phiWall											! add the amount of scalar that has been absorbed at the current location in the current direction

!------------------------------------------------
END SUBROUTINE AbsorbedScalarS2
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE AbsorbedScalarV(i,j,k,m,phiBC)		! measures the total absorbed scalar
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m				! index variables
REAL(dbl), INTENT(IN) :: phiBC     				! scalar contribution from the boundary condition
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

phiIN 	= phiBC																						! contribution from the wall to the crrent node (in)
phiOUT	= (fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)		! contribution to the wall from the current node (out)

phiAbsorbedV = phiAbsorbedV + (phiOUT - phiIN)												! add the amount of scalar that has been absorbed at the current location in the current direction

!------------------------------------------------
END SUBROUTINE AbsorbedScalarV
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarInOut								! measure scalar that has left or entered the domain through the inlet or outlet
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1,iComm		! index variables
REAL(dbl) :: phiOUT, phiIN							! scalar values exchanged with the wall

! XY Faces (z-direction)
DO iComm=5,6

  IF(SubID(iComm) .EQ. 0) THEN					

    k = XY_SendIndex(iComm)						! k index

    DO j=1,nySub
      DO i=1,nxSub

        IF((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN

          DO m=1,NumFs_face

            ! i,j,k location of neighboring node
            im1 = i - ex(bb(f_Comps(iComm,m)))
            jm1 = j - ey(bb(f_Comps(iComm,m)))
            km1 = k - ez(bb(f_Comps(iComm,m)))

            IF((node(im1,jm1,km1) .EQ. SOLID).OR.(node(im1,jm1,km1) .EQ. SOLID2)) THEN
              phiIN = (fplus(bb(f_Comps(iComm,m)),im1,jm1,km1)/rho(im1,jm1,km1) - wt(bb(f_Comps(iComm,m)))*Delta)	&								! scalar contribution from inlet/outlet to current node
                      *phiTemp(im1,jm1,km1)		
              phiOUT	= (fplus(f_Comps(iComm,m),i,j,k)/rho(i,j,k) - wt(f_Comps(iComm,m))*Delta)*phiTemp(i,j,k)										! scalar contribution from current node to inlet/outlet
              phiInOut = phiInOut + (phiOUT - phiIN)
            END IF

          END DO
 
        END IF

      END DO
    END DO
  
  END IF

END DO

!------------------------------------------------
END SUBROUTINE ScalarInOut
!------------------------------------------------
!!------------------------------------------------
!SUBROUTINE Interp_ParToNodes_Conc ! Interpolate Particle concentration release to node locations 
!! Called by Particle_Track (LBM.f90) to get delphi_particle
!!------------------------------------------------
!IMPLICIT NONE
!INTEGER(lng)  :: i,ix0,ix1,iy0,iy1,iz0,iz1
!REAL(dbl)     :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd
!
!DO i=1,np
!	!ix0=FLOOR(xp(i))
!	!ix1=CEILING(xp(i))
!	!iy0=FLOOR(yp(i))
!	!iy1=CEILING(yp(i))
!	!iz0=FLOOR(zp(i))
!	!iz1=CEILING(zp(i))
!
!	ix0=FLOOR(xp(i))
!	ix1=FLOOR(xp(i))+1_lng
!	iy0=FLOOR(yp(i))
!	iy1=FLOOR(yp(i))+1_lng
!	iz0=FLOOR(zp(i))
!	iz1=FLOOR(zp(i))+1_lng
!
!	xd=(xp(i)-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
!	yd=(yp(i)-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
!	zd=(zp(i)-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
!
!	delphi_particle(ix0,iy0,iz0)=delNBbyCV(i)
!	
!!	!yd=0.0_dbl ! TEST: used to keep particle motion plainly 2-D ! Balaji added
!!
!!	! u-interpolation
!!	! Do first level linear interpolation in x-direction
!!	c00 = u(ix0,iy0,iz0)*(1.0_dbl-xd)+u(ix1,iy0,iz0)*xd	
!!	c01 = u(ix0,iy0,iz1)*(1.0_dbl-xd)+u(ix1,iy0,iz1)*xd	
!!	c10 = u(ix0,iy1,iz0)*(1.0_dbl-xd)+u(ix1,iy1,iz0)*xd	
!!	c11 = u(ix0,iy1,iz1)*(1.0_dbl-xd)+u(ix1,iy1,iz1)*xd	
!!	! Do second level linear interpolation in y-direction
!!	c0  = c00*(1.0_dbl-yd)+c10*yd
!!	c1  = c01*(1.0_dbl-yd)+c11*yd
!!	! Do third level linear interpolation in z-direction
!!	c   = c0*(1.0_dbl-zd)+c1*zd
!!        up(i)=c
!!
!!	! v-interpolation
!!	! Do first level linear interpolation in x-direction
!!	c00 = v(ix0,iy0,iz0)*(1.0_dbl-xd)+v(ix1,iy0,iz0)*xd	
!!	c01 = v(ix0,iy0,iz1)*(1.0_dbl-xd)+v(ix1,iy0,iz1)*xd	
!!	c10 = v(ix0,iy1,iz0)*(1.0_dbl-xd)+v(ix1,iy1,iz0)*xd	
!!	c11 = v(ix0,iy1,iz1)*(1.0_dbl-xd)+v(ix1,iy1,iz1)*xd	
!!	! Do second level linear interpolation in y-direction
!!	c0  = c00*(1.0_dbl-yd)+c10*yd
!!	c1  = c01*(1.0_dbl-yd)+c11*yd
!!	! Do third level linear interpolation in z-direction
!!	c   = c0*(1.0_dbl-zd)+c1*zd
!!        vp(i)=c
!!
!!	! w-interpolation
!!	! Do first level linear interpolation in x-direction
!!	c00 = w(ix0,iy0,iz0)*(1.0_dbl-xd)+w(ix1,iy0,iz0)*xd	
!!	c01 = w(ix0,iy0,iz1)*(1.0_dbl-xd)+w(ix1,iy0,iz1)*xd	
!!	c10 = w(ix0,iy1,iz0)*(1.0_dbl-xd)+w(ix1,iy1,iz0)*xd	
!!	c11 = w(ix0,iy1,iz1)*(1.0_dbl-xd)+w(ix1,iy1,iz1)*xd	
!!	! Do second level linear interpolation in y-direction
!!	c0  = c00*(1.0_dbl-yd)+c10*yd
!!	c1  = c01*(1.0_dbl-yd)+c11*yd
!!	! Do third level linear interpolation in z-direction
!!	c   = c0*(1.0_dbl-zd)+c1*zd
!!        wp(i)=c
!!
!!	!up(i)=0.5*(u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+u(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
!!	!vp(i)=0.5*(v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+v(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
!!	!wp(i)=0.5*(w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))+w(CEILING(xp(i)),CEILING(yp(i)),CEILING(zp(i))))
!!	!up(i)=u(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
!!	!vp(i)=v(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
!!	!wp(i)=w(FLOOR(xp(i)),FLOOR(yp(i)),FLOOR(zp(i)))
!ENDDO
!!------------------------------------------------
!END SUBROUTINE Interp_ParToNodes_Conc ! Interpolate Particle concentration release to node locations 
!!------------------------------------------------

!================================================
END MODULE PassiveScalar
!================================================
