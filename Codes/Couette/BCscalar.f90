!==================================================================================================
MODULE BCscalar		! Sets Scalar Boundary Conditions
!==================================================================================================
USE SetPrecision
USE Setup  
USE MPI
USE IC
USE BClbm

CONTAINS


!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar_Fixed_Scalar(m,i,j,k,im1,jm1,km1,phiBC)				! implements the scalar BCs 
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1				! index variables
REAL(dbl), INTENT(OUT) :: phiBC     					! scalar contribution from the boundary condition
INTEGER(lng) :: ip1,jp1,kp1 						! neighboring nodes (2 away from the wall)
REAL(dbl) :: q, rhoB,phiB,feq_m						! distance ratio from the current node to the solid node
REAL(dbl) :: rhoAstar,phiAstar, PkAstar,feq_Astar,feq_Bstar 		! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: fPlusBstar, rhoBstar, phiBstar, PkBstar 			! Values interpolated to Bstar location
REAL(dbl) :: phiijk_m							! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)
REAL(dbl) :: cosTheta, sinTheta						! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb, ubb,vbb,wbb							! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk 							! radius of the solid node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt
INTEGER(lng) :: it
REAL(dbl) :: dPhiDn ! Gradient at the wall
REAL(dbl) :: Pw !Permeability

Pw = 2.0*Dm !Permeability is set equal to 2*Dm

CALL qCalcFarhad(i,q)		

IF (rijk .GE. rOut(k)) THEN
   ub = velOut(km1) 			
   vb = 0.0_dbl 				
   wb = 0.0_dbl 							
ELSE IF (rijk .LE. rIn(k)) THEN
   ub = velIn(km1) 				
   vb = 0.0_dbl					
   wb = 0.0_dbl 					
END IF		

!----- neighboring node (fluid side)	
ip1 = i + ex(m) 			
jp1 = j + ey(m)			
kp1 = k + ez(m)		

!------ if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF	

phiWall= ( (phiTemp(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) - (phiTemp(ip1,j,k)*q*q/(1.0+2.0*q)) ) / ( 1.0 + (Pw/Dm)*q*(1+q)/(1+2.0*q) )	! calculate phiWall for flux BC (eq. 28 in paper)

!----- Computing values at A* & scalar streamed from A* (Chpter 3 paper)
rhoAstar= (rho(i,j,k)- rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)	! extrapolate the density
CALL Equilibrium_LOCAL(m,rhoAstar,ub,vb,wb,feq_Astar)    		! calculate the equibrium distribution function in the mth direction
phiAstar= phiWall							! getting phi at the solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar			! contribution from the wall in mth direction (0 if phiWall=0)

!------ Computing values at B* & scalar streamed from B* (Chpter 3 paper)
rhoBstar=   (1-q)*rho(ip1,jp1,kp1)     + q* rho(i,j,k)
phiBstar=   (1-q)*phiTemp(ip1,jp1,kp1) + q* phiTemp(i,j,k)
fPlusBstar= (1-q)*fplus(m,ip1,jp1,kp1) + q* fplus(m,i,j,k)
PkBstar=    (fplusBstar/rhoBstar - wt(m)*Delta)*phiBstar


!------------------------------------------------
END SUBROUTINE Scalar_Fixed_Scalar
!------------------------------------------------






!--------------------------------------------------------------------------------------------------
SUBROUTINE Scalar_Fixed_Flux(m,i,j,k,im1,jm1,km1,phiBC)				! implements the scalar BCs 
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1				! index variables
REAL(dbl), INTENT(OUT) :: phiBC     					! scalar contribution from the boundary condition
INTEGER(lng) :: ip1,jp1,kp1 						! neighboring nodes (2 away from the wall)
REAL(dbl) :: q, rhoB,phiB,feq_m						! distance ratio from the current node to the solid node
REAL(dbl) :: rhoAstar,phiAstar, PkAstar,feq_Astar,feq_Bstar 		! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: fPlusBstar, rhoBstar, phiBstar, PkBstar 			! Values interpolated to Bstar location
REAL(dbl) :: phiijk_m							! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)
REAL(dbl) :: cosTheta, sinTheta						! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb, ubb,vbb,wbb							! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk 							! radius of the solid node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt
INTEGER(lng) :: it

CALL qCalcFarhad(i,q)		

!IF ((j.EQ.21).AND.(k.EQ.3)) THEN
!   write(*,*) iter,i,q
!END IF

IF (rijk .GE. rOut(k)) THEN
   ub = velOut(km1) 			
   vb = 0.0_dbl 				
   wb = 0.0_dbl 							
ELSE IF (rijk .LE. rIn(k)) THEN
   ub = velIn(km1) 				
   vb = 0.0_dbl					
   wb = 0.0_dbl 					
END IF		

!----- neighboring node (fluid side)	
ip1 = i + ex(m) 			
jp1 = j + ey(m)			
kp1 = k + ez(m)		

!------ if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF	

ubb= ub		 !0.0
vbb= vb 	!0.0
wbb= wb 	!0.0

!----- Computing values at A* & scalar streamed from A* (Chpter 3 paper)
rhoAstar= (rho(i,j,k)- rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)	! extrapolate the density
CALL Equilibrium_LOCAL(m,rhoAstar,ubb,vbb,wbb,feq_Astar)		! calculate the equibrium distribution function in the mth direction
phiAstar= phiWall							! getting phi at the solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar			! contribution from the wall in mth direction (0 if phiWall=0)

!------ Computing values at B* & scalar streamed from B* (Chpter 3 paper)
rhoBstar=   (1-q)*rho(ip1,jp1,kp1)     + q*rho(i,j,k)
phiBstar=   (1-q)*phi(ip1,jp1,kp1)     + q*phi(i,j,k)

fPlusBstar= (1-q)*fplus(m,ip1,jp1,kp1) + q*fplus(m,i,j,k)
PkBstar=    (fplusBstar/rhoBstar - wt(m)*Delta)*phiBstar

CALL Equilibrium_LOCAL(m,rhoBstar,ubb,vbb,wbb,feq_Bstar)
!PkBstar=    (feq_Bstar/rhoBstar - wt(m)*Delta)*phiBstar

phiBC= PkAstar+ (PkAstar- PkBstar)*(1-q)

!IF ((j.EQ.21).AND.(k.EQ.3)) THEN
!   write(*,*) iter, i, feq_Bstar, fPlusBstar, PkBstar 
!END IF


!------------------------------------------------
END SUBROUTINE Scalar_Fixed_Flux
!------------------------------------------------




!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarBC2(m,i,j,k,im1,jm1,km1,phiBC,phiOut,phiIn)								! implements the scalar BCs by Balaji Dec 2014.
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1								! index variables
REAL(dbl), INTENT(OUT) :: phiBC,phiOut,phiIn     											! scalar contribution from the boundary condition
INTEGER(lng) :: ip1,jp1,kp1 														! neighboring nodes (2 away from the wall)
INTEGER(lng) :: ip2,jp2,kp2 														! neighboring nodes (3 away from the wall)
REAL(dbl) :: q																! distance ratio from the current node to the solid node
REAL(dbl) :: rhoB,phiB															! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: rhoAst,PAstToBSt,ScAst													! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes to an interior node 1 lattice unit away. 
REAL(dbl) :: rhoBst,PBstToCSt,ScBst,fBst												! values of density and at the boundary, and contribution of scalar from the B* node to C* node (wee yanxing wang's 2010 paper)
REAL(dbl) :: PBstToASt,fBstopp														! values of density and at the boundary, and contribution of scalar from the B* node to A* node (see yanxing wang's 2010 paper)
REAL(dbl) :: PAtoB, PAtoO, PBtoA, PCtoB, PBtoC, PCtoD,POtoA, Pmovingwall,fbb,fmoving															! values of density and at the boundary, and contribution of scalar from the B* node to A* node (see yanxing wang's 2010 paper)
REAL(dbl) :: rhoX,ScX
REAL(dbl) :: feq_m,feq_bbm,feq1_m,feq2_m,fnoneq_m,fnoneq1_m,fnoneq2_m									! equilibrium distribution function in the mth direction
REAL(dbl) :: phiijk_m	! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)
REAL(dbl) :: phiAst	! Scalar value at the boundary - For both Dirichlet and Neumann BC. For Dirichlet BC it is phiWall, but for Neumann BC it is computed from the flux BC.  
REAL(dbl) :: dphidn	! Neumann flux BC.  
REAL(dbl) :: cosTheta, sinTheta													! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb																! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk													! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt				! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 
INTEGER(lng) :: it			! loop index variables

LOGICAL :: BC2FLAG

! neighboring node (fluid side)	
ip1 = i + ex(m) 	! i + 1
jp1 = j + ey(m)		! j + 1
kp1 = k + ez(m)		! k + 1
ip2 = ip1 + ex(m) 	! i + 2
jp2 = jp1 + ey(m)	! j + 2
kp2 = kp1 + ez(m)	! k + 2

BC2FLAG = .FALSE.
IF ((node(ip1,jp1,kp1) .EQ. FLUID).AND.(node(ip2,jp2,kp2) .EQ. FLUID)) THEN
	BC2FLAG = .TRUE.
END IF
BC2FLAG = .FALSE.

IF (BC2FLAG) THEN
! Calculate q using Yanxing's method
!*****************************************************************************
		!rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
		rijk = x(im1)								! height at current location
		 ! Initial fluid node guess
                 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                
		 ! Initial solid node guess
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
                 
	 IF (k.NE.km1) THEN
                 DO it=1,qitermax
		   ! guess of boundary location 
                   xt=(x1+x2)/2.0_dbl
                   yt=(y1+y2)/2.0_dbl
                   zt=(z1+z2)/2.0_dbl

      		   !rt = SQRT(xt*xt + yt*yt)
      		   rt = xt
		   !Write(*,*) 'test'
                   IF (node(im1,jm1,km1).EQ.SOLID2) THEN
		   	!ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
		   	ht = ((zt-z(k))*rOut(km1)+(z(km1)-zt)*rOut(k))/(z(km1)-z(k))
		   	!ht = (r(km1)+r(k))/2.0_dbl
		   ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
		   	!ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
		   	ht = ((zt-z(k))*rIn(km1)+(z(km1)-zt)*rIn(k))/(z(km1)-z(k))
		   	!ht = (r(km1)+r(k))/2.0_dbl
		   END IF

                   IF(rt.GT.ht) then
                     x2=xt
                     y2=yt
                     z2=zt
                   ELSE
                     x1=xt
                     y1=yt
                     z1=zt
                   END IF
				   
                 END DO
		 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                 
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
 
                 q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
		 !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
	 ELSE
		  DO it=1,qitermax
		   ! guess of boundary location 
                   xt=(x1+x2)/2.0_dbl
                   yt=(y1+y2)/2.0_dbl
                   zt=(z1+z2)/2.0_dbl

      		   !rt = SQRT(xt*xt + yt*yt)
      		   rt = xt
		   !Write(*,*) 'test'
		   IF (node(im1,jm1,km1).EQ.SOLID2) THEN
			   !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
			   !ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
			   ht = (rOut(km1)+rOut(k))/2.0_dbl
		   ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
			   !ht = (ABS(zt-z(k))*r(km1)+ABS(z(km1)-zt)*r(k))/ABS(z(km1)-z(k))
			   !ht = ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
			   ht = (rIn(km1)+rIn(k))/2.0_dbl
		   END IF

                   IF(rt.GT.ht) then
                     x2=xt
                     y2=yt
                     z2=zt
                   ELSE
                     x1=xt
                     y1=yt
                     z1=zt
                   END IF
				   
                 END DO
		 x1=x(i)
                 y1=y(j)
                 z1=z(k)
                 
                 x2=x(im1)
                 y2=y(jm1)
                 z2=z(km1)
 
                 q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
		 !write(*,*) 'q',q,zt,z1,z2,0.5*(z1+z2),rt,ht
	 ENDIF
		 cosTheta=xt/rt
		 sinTheta=yt/rt
	 IF (k.NE.km1) THEN
		   IF (node(im1,jm1,km1).EQ.SOLID2) THEN
			 !vt = (ABS(zt-z(k))*vel(km1)+ABS(z(km1)-zt)*vel(k))/ABS(z(km1)-z(k))
			 vt = ((zt-z(k))*velOut(km1)+(z(km1)-zt)*velOut(k))/(z(km1)-z(k))
		   ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
			 !vt = (ABS(zt-z(k))*vel(km1)+ABS(z(km1)-zt)*vel(k))/ABS(z(km1)-z(k))
			 vt = ((zt-z(k))*velIn(km1)+(z(km1)-zt)*velIn(k))/(z(km1)-z(k))
		   END IF
	 ELSE
		   IF (node(im1,jm1,km1).EQ.SOLID2) THEN
			 vt = (velOut(k)+velOut(km1))*0.5_dbl
		   ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
			 vt = (velIn(k)+velIn(km1))*0.5_dbl
		   END IF
	 ENDIF
		ub =  velIn(km1) !0.0_dbl!0.0!vel(km1)*cosTheta		! x-component of the velocity at i,j,k
		vb = 0.0_dbl	!0.0!vel(km1)*sinTheta			! y-component of the velocity at i,j,k
		wb = 0.0_dbl   	!vt!vel(km1)!0.0_dbl			! only z-component of velocity	
                write(*,*) 'ub', ub
	        ! make sure 0<q<1
        IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
          OPEN(1000,FILE="error.txt")
          WRITE(1000,*) "q=",q
          WRITE(1000,*) "m=",m
          WRITE(1000,*) "i=",i,"j=",j,"k=",k
          WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
          CLOSE(1000)
          STOP
        END IF	
ELSE
        q=0.5_dbl
	rijk = x(im1)								! height at current location
        
        !cosTheta = x(im1)/rijk!r(km1)	! COS(theta)
        !sinTheta = y(jm1)/rijk!r(km1)	! SIN(theta)
        
        
	IF (rijk .GE. rOut(k)) THEN
		ub =  velOut(km1) 	!0.0_dbl!0.0!vel(km1)*cosTheta			! x-component of the velocity at i,j,k
		vb = 0.0_dbl		!0.0!vel(km1)*sinTheta				! y-component of the velocity at i,j,k
		wb = 0.0_dbl 		!velOut(km1)!vel(km1)!0.0_dbl			! only z-component in this case			
	ELSE IF (rijk .LE. rIn(k)) THEN
		ub = velIn(km1) !0.0_dbl!0.0!vel(km1)*cosTheta					! x-component of the velocity at i,j,k
		vb = 0.0_dbl	!0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
		wb = 0.0_dbl 	!velIn(km1)!vel(km1)!0.0_dbl						! only z-component in this case	
	END IF			


END IF

! Using fbb to bounce back
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN

    IF (BC2FLAG) THEN

            rhoX = 2.0_dbl*q*rho(i,j,k) + (1.0_dbl-2.0_dbl*q)*rho(ip1,jp1,kp1)
            ScX  =  2.0_dbl*q*phiTemp(i,j,k) + (1.0_dbl-2.0_dbl*q)*phiTemp(ip1,jp1,kp1)
            rhoAst =  (rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)
            ScAst =  (phiTemp(i,j,k) - phiTemp(ip1,jp1,kp1))*(q) + phiTemp(i,j,k)
        
        ! 2nd order Lallemand and Luo
            !fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            !fmoving = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            !fmoving = 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))
            fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) &
                + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) & 
                - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) &
                + fmoving

    ELSE

        ! zeroth order half-way bounce back
            rhoX = rho(i,j,k)
            ScX  = phiTemp(i,j,k)
            rhoAst = rho(i,j,k)
            ScAst =  phiTemp(i,j,k)
            !fmoving = (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
            fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))!*(rhoAst/ScAst)
            !fmoving = 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m)))
            fbb = fplus(bb(m),i,j,k) &
          	  + fmoving

    END IF

  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN

    IF (BC2FLAG) THEN

            rhoX = 2.0_dbl*q*rho(i,j,k) + (1.0_dbl-2.0_dbl*q)*rho(ip1,jp1,kp1)
            ScX  =  2.0_dbl*q*phiTemp(i,j,k) + (1.0_dbl-2.0_dbl*q)*phiTemp(ip1,jp1,kp1)
            rhoAst =  (rho(i,j,k) - rho(ip1,jp1,kp1))*(q) + rho(i,j,k)
            ScAst =  (phiTemp(i,j,k) - phiTemp(ip1,jp1,kp1))*(q) + phiTemp(i,j,k)
        
        ! 2nd order Lallemand and Luo
            fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
            fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 	&
                + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q	&
                - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1) &
                + fmoving

    ELSE
        
        ! zeroth order half-way bounce back
            rhoX = rho(i,j,k)
            ScX  = phiTemp(i,j,k)
            rhoAst = rho(i,j,k)
            ScAst =  phiTemp(i,j,k)
            fmoving = (6.0_dbl*wt(m)*rhoAst*(ub*ex(m) + vb*ey(m) + wb*ez(m)))!*(rhoAst/ScAst)
            fbb = fplus(bb(m),i,j,k) &
        	+ fmoving

    END IF

  ELSE
    OPEN(1000,FILE='error-'//sub//'.txt')
    WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
    CLOSE(1000)
    STOP
  END IF



phiBC = 1.0_dbl*(1.0_dbl*(fbb-fmoving)/rhoX - wt(bb(m))*Delta)*ScX+(fmoving/rhoAst)*ScAst !+ wt(bb(m))*Delta*phiTemp(i,j,k)
phiOut = 1.0_dbl*(1.0_dbl*fplus(bb(m),i,j,k)/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)!POtoA*q + (1.0_dbl-q)*PAtoB ! phi going out of the surface
phiIn =  phiBC!PAtoO*q + (1.0_dbl-q)*PBtoA ! phi going into the surface

!------------------------------------------------
END SUBROUTINE ScalarBC2
!------------------------------------------------














































!--------------------------------------------------------------------------------------------------
SUBROUTINE GetPhiWallNew(ip0,jp0,kp0,m,q,dphidn,phiAst,rhoAst)	! Get phiAst for specified flux  New (Balaji, Jan 2015)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: ip0,jp0,kp0,m	! index variables
REAL(dbl), INTENT(IN) :: q     	! Boundary distance fraction
REAL(dbl) :: qlocal     	! Boundary distance fraction (local variable)
REAL(dbl), INTENT(IN) :: dphidn ! Specified boundary flux (gradient)
REAL(dbl), INTENT(OUT) :: phiAst,rhoAst ! scalar value at the wall
REAL(dbl):: phiA,phiB,phiC ! scalar values at the wall neighbours
REAL(dbl):: rhoA,rhoB,rhoC ! scalar values at the wall neighbours
REAL(dbl) :: delta ! lattice distance in direction m
REAL(dbl) :: nx1,ny1,nz1 ! normal direction 1
REAL(dbl) :: px1,py1,pz1 ! surface point 1
REAL(dbl) :: px2,py2,pz2 ! surface point 2
REAL(dbl) :: px3,py3,pz3 ! surface point 3
REAL(dbl) :: px4,py4,pz4 ! surface point 4
REAL(dbl) :: magtvec,magnvec ! magnitude of vectors
REAL(dbl) :: dx,dy,dz,dist  ! direction vector components
INTEGER(lng):: i,j,k,im1,jm1,km1,ip1,jp1,kp1,ip2,jp2,kp2	! index variables
REAL(dbl),DIMENSION(3,3) :: arr,arrinv
REAL(dbl),DIMENSION(3) :: brhs
REAL(dbl) :: epsx,epsy,epsz,epsmin,veclen ! temp variables


! Solid node
im1 = ip0 - ex(m)
jm1 = jp0 - ey(m)
km1 = kp0 - ez(m)
! neighboring node (fluid side)	
ip1 = ip0 + ex(m) 	! i + 1
jp1 = jp0 + ey(m)	! j + 1
kp1 = kp0 + ez(m)	! k + 1
ip2 = ip1 + ex(m) 	! i + 2
jp2 = jp1 + ey(m)	! j + 2
kp2 = kp1 + ez(m)	! k + 2


!open(9,file='testoutput.dat',position='append')

! get normal direction 
! to do this, we will need at the least 3 points on the surface plane. 

!! Approach 1: We get the 3-4 surface points and normal vector from analytical surface equations
!!!! Compute normal from analytical surface
epsmin = 1.0e-15_dbl


! Get first surface point from Ray tracing
epsx=0.0_dbl
epsy=0.0_dbl
epsz=0.0_dbl
CALL IntersectSurface(m,ip0,jp0,kp0,im1,jm1,km1,px1,py1,pz1,qlocal,epsz,epsy,epsx)

nx1 = 2.0_dbl*px1
ny1 = 2.0_dbl*py1
IF (k.GT.1) THEN
	IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	nz1 = -(2.0_dbl*((rOut(k)-rOut(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*rOut(k-1)/zcf+(pz1-(z(k-1)/zcf))*rOut(k)/zcf)/((z(k)-z(k-1))/zcf))
	ELSE
	nz1 = -(2.0_dbl*((rIn(k)-rIn(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*rIn(k-1)/zcf+(pz1-(z(k-1)/zcf))*rIn(k)/zcf)/((z(k)-z(k-1))/zcf))
	ENDIF
ELSE 
	IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	nz1 = -(2.0_dbl*((rOut(k+1)-rOut(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*rOut(k)/zcf+(pz1-(z(k)/zcf))*rOut(k+1)/zcf)/((z(k+1)-z(k))/zcf))
	ELSE
	nz1 = -(2.0_dbl*((rIn(k+1)-rIn(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*rIn(k)/zcf+(pz1-(z(k)/zcf))*rIn(k+1)/zcf)/((z(k+1)-z(k))/zcf))
	END IF

ENDIF
magnvec = SQRT(nx1**2+ny1**2+nz1**2)
nx1 = nx1/magnvec
ny1 = ny1/magnvec
nz1 = nz1/magnvec


! Get points along the normal vector in the fluid side of the surface. 
veclen = 0.1_dbl
px2 = px1 - veclen*nx1
py2 = py1 - veclen*ny1
pz2 = pz1 - veclen*nz1
CALL InterpolateProp(px2,py2,pz2,phiA,rhoA)

px3 = px2 - veclen*nx1
py3 = py2 - veclen*ny1
pz3 = pz2 - veclen*nz1
CALL InterpolateProp(px3,py3,pz3,phiB,rhoB)
!check if these are fluid points. 

delta = veclen
phiAst = (4.0_dbl*phiA - 1.0_dbl*phiB - 2.0*delta*dphidn)/3.0_dbl
rhoAst = (4.0_dbl*rhoA - 1.0_dbl*rhoB - 2.0*delta*0.0_dbl)/3.0_dbl

!------------------------------------------------
END SUBROUTINE GetPhiWallNew
!------------------------------------------------










!------------------------------------------------
SUBROUTINE InterpolateProp(xp,yp,zp,phip,rhop) ! Interpolate solution field w to location xp,yp,zp 
!------------------------------------------------
IMPLICIT NONE
REAL(dbl), INTENT(IN) :: xp,yp,zp	! coordinates
REAL(dbl) :: xps,yps,zps	! coordinates
REAL(dbl), INTENT(OUT) :: phip,rhop   	! Interpolated field value being returned 
INTEGER(lng) :: ix0,iy0,iz0,ix1,iy1,iz1,numFLUIDs
REAL(dbl) :: xd,yd,zd
REAL(dbl) :: c00,c01,c10,c11,c,c0,c1
REAL(dbl) :: rhoSum,phiSum

	xps = xp-REAL(iMin-2_lng)+REAL(Ci - 1_lng)
	ix0=FLOOR(xps) 
	ix1=ix0 + 1_lng
        yps = yp-REAL(jMin-2_lng)+REAL(Cj - 1_lng)
	iy0=FLOOR(yps)
	iy1=iy0 + 1_lng
        zps = zp + 0.5_dbl- REAL(kMin-1_lng)
	iz0=FLOOR(zps)
	iz1=iz0+1_lng
	xd=(xps-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))	
	yd=(yps-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))	
	zd=(zps-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))


	! w-interpolation
	! Do first level linear interpolation in x-direction
	IF ((node(ix0,iy0,iz0).EQ. FLUID).OR.(node(ix1,iy0,iz0).EQ. FLUID) &
	.OR.(node(ix0,iy0,iz1).EQ. FLUID).OR.(node(ix1,iy0,iz1).EQ. FLUID) &
	.OR.(node(ix0,iy1,iz0).EQ. FLUID).OR.(node(ix1,iy1,iz0).EQ. FLUID) &
	.OR.(node(ix0,iy1,iz1).EQ. FLUID).OR.(node(ix1,iy1,iz1).EQ. FLUID)) THEN
		c00 = phiTemp(ix0,iy0,iz0)*(1.0_dbl-xd)+phiTemp(ix1,iy0,iz0)*xd	
		c01 = phiTemp(ix0,iy0,iz1)*(1.0_dbl-xd)+phiTemp(ix1,iy0,iz1)*xd	
		c10 = phiTemp(ix0,iy1,iz0)*(1.0_dbl-xd)+phiTemp(ix1,iy1,iz0)*xd	
		c11 = phiTemp(ix0,iy1,iz1)*(1.0_dbl-xd)+phiTemp(ix1,iy1,iz1)*xd	
		! Do second level linear interpolation in y-direction
		c0  = c00*(1.0_dbl-yd)+c10*yd
		c1  = c01*(1.0_dbl-yd)+c11*yd
		! Do third level linear interpolation in z-direction
		c   = c0*(1.0_dbl-zd)+c1*zd
	        phip=c

		c00 = rho(ix0,iy0,iz0)*(1.0_dbl-xd)+rho(ix1,iy0,iz0)*xd	
		c01 = rho(ix0,iy0,iz1)*(1.0_dbl-xd)+rho(ix1,iy0,iz1)*xd	
		c10 = rho(ix0,iy1,iz0)*(1.0_dbl-xd)+rho(ix1,iy1,iz0)*xd	
		c11 = rho(ix0,iy1,iz1)*(1.0_dbl-xd)+rho(ix1,iy1,iz1)*xd	
		! Do second level linear interpolation in y-direction
		c0  = c00*(1.0_dbl-yd)+c10*yd
		c1  = c01*(1.0_dbl-yd)+c11*yd
		! Do third level linear interpolation in z-direction
		c   = c0*(1.0_dbl-zd)+c1*zd
	        rhop=c
	ELSE
		OPEN(1000,FILE='error-'//sub//'.txt')
		WRITE(1000,*) "Error in InterpolateProp in ICBC.f90 (line 1736): rhos are zero"
		WRITE(1000,*) "rho(i,j,k) is zero",rho(ix0,iy0,iz0),rho(ix1,iy0,iz0),rho(ix0,iy0,iz1) &
		,rho(ix1,iy0,iz1),rho(ix0,iy1,iz0),rho(ix1,iy1,iz0),rho(ix0,iy1,iz1),rho(ix1,iy1,iz1)
		WRITE(1000,*) "node(i,j,k) is ",node(ix0,iy0,iz0),node(ix1,iy0,iz0),node(ix0,iy0,iz1) &
		,node(ix1,iy0,iz1),node(ix0,iy1,iz0),node(ix1,iy1,iz0),node(ix0,iy1,iz1),node(ix1,iy1,iz1)
		WRITE(1000,*) "node index",ix0,iy0,iz0,ix1,iy1,iz1
		CLOSE(1000)
		STOP
	ENDIF
         

!------------------------------------------------
END SUBROUTINE InterpolateProp
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE GetPhiWall(ip0,jp0,kp0,m,q,dphidn,phiAst)	! Get phiAst for specified flux  New (Balaji, Dec 2014)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: ip0,jp0,kp0,m	! index variables
REAL(dbl), INTENT(IN) :: q     	! Boundary distance fraction
REAL(dbl) :: qlocal     	! Boundary distance fraction (local variable)
REAL(dbl), INTENT(IN) :: dphidn ! Specified boundary flux (gradient)
REAL(dbl), INTENT(OUT) :: phiAst ! scalar value at the wall
!REAL(dbl), INTENT(IN) :: phiA,phiB,phiC ! scalar values at the wall neighbours
REAL(dbl):: phiA,phiB,phiC ! scalar values at the wall neighbours
REAL(dbl) :: dphidtAst1,dphidtA1,dphidtB1,dphidtC1 ! tangential gradient of phi (1st tangent vector)
REAL(dbl) :: dphidtAst2,dphidtA2,dphidtB2,dphidtC2 ! tangential gradient of phi (2nd tangent vector)
REAL(dbl) :: delta ! lattice distance in direction m
REAL(dbl) :: tx1,ty1,tz1 ! tangent direction 1
REAL(dbl) :: tx2,ty2,tz2 ! tangent direction 2
REAL(dbl) :: nx1,ny1,nz1 ! normal direction 1
REAL(dbl) :: px1,py1,pz1 ! surface point 1
REAL(dbl) :: px2,py2,pz2 ! surface point 2
REAL(dbl) :: px3,py3,pz3 ! surface point 3
REAL(dbl) :: px4,py4,pz4 ! surface point 4
REAL(dbl) :: sx1,sy1,sz1 ! surface vector 1
REAL(dbl) :: sx2,sy2,sz2 ! surface vector 2
REAL(dbl) :: magtvec,magnvec ! magnitude of vectors
REAL(dbl) :: dx,dy,dz,dist  ! direction vector components
REAL(dbl) :: epsx,epsy,epsz,eps,epsmin  ! deviation vector components
REAL(dbl) :: dphidx,dphidy,dphidz,dphidk ! grad(phi) components
INTEGER(lng):: i,j,k,im1,jm1,km1,ip1,jp1,kp1,ip2,jp2,kp2	! index variables
REAL(dbl),DIMENSION(3,3) :: arr,arrinv
REAL(dbl),DIMENSION(3) :: brhs
REAL(dbl) :: temp,aa1,aa2,aa3,bb1,bb2,bb3,cc1,cc2,cc3,rb1,rb2,rb3,t1,t2,sz,sx,sy


! Solid node
im1 = ip0 - ex(m)
jm1 = jp0 - ey(m)
km1 = kp0 - ez(m)
! neighboring node (fluid side)	
ip1 = ip0 + ex(m) 	! i + 1
jp1 = jp0 + ey(m)	! j + 1
kp1 = kp0 + ez(m)	! k + 1
ip2 = ip1 + ex(m) 	! i + 2
jp2 = jp1 + ey(m)	! j + 2
kp2 = kp1 + ez(m)	! k + 2


!open(9,file='testoutput.dat',position='append')

! get normal direction 
! to do this, we will need at the least 3 points on the surface plane. 

!! Approach 1: We get the 3-4 surface points and normal vector from analytical surface equations
!!!! Compute normal from analytical surface
epsmin = 1.0e-15_dbl

! Get first surface point from Ray tracing
epsx=0.0_dbl
epsy=0.0_dbl
epsz=0.0_dbl
CALL IntersectSurface(m,ip0,jp0,kp0,im1,jm1,km1,px1,py1,pz1,qlocal,epsz,epsy,epsx)

nx1 = 2.0_dbl*px1
ny1 = 2.0_dbl*py1
IF (k.GT.1) THEN
	IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	nz1 = -(2.0_dbl*((rOut(k)-rOut(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*rOut(k-1)/zcf+(pz1-(z(k-1)/zcf))*rOut(k)/zcf)/((z(k)-z(k-1))/zcf))
	ELSE
	nz1 = -(2.0_dbl*((rIn(k)-rIn(k-1))/(z(k)-z(k-1)))*(((z(k)/zcf)-pz1)*rIn(k-1)/zcf+(pz1-(z(k-1)/zcf))*rIn(k)/zcf)/((z(k)-z(k-1))/zcf))
	ENDIF
ELSE 
	IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	nz1 = -(2.0_dbl*((rOut(k+1)-rOut(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*rOut(k)/zcf+(pz1-(z(k)/zcf))*rOut(k+1)/zcf)/((z(k+1)-z(k))/zcf))
	ELSE
	nz1 = -(2.0_dbl*((rIn(k+1)-rIn(k))/(z(k+1)-z(k)))*(((z(k+1)/zcf)-pz1)*rIn(k)/zcf+(pz1-(z(k)/zcf))*rIn(k+1)/zcf)/((z(k+1)-z(k))/zcf))
	END IF

ENDIF
magnvec = SQRT(nx1**2+ny1**2+nz1**2)
nx1 = nx1/magnvec
ny1 = ny1/magnvec
nz1 = nz1/magnvec

! get 2nd point analytically
epsx = 0.0001_dbl
IF (abs(ny1+nz1).GT.epsmin) THEN
	epsy = -(nx1*epsx)/(ny1+nz1)
ELSE
	epsy = -(nx1*epsx)/sign(epsmin,ny1+nz1)
ENDIF
px2 = px1 + epsx
py2 = py1 + epsy
pz2 = pz1 + epsz


! Calculate first surface vector after obtaining the points
sx1 = px1-px2
sy1 = py1-py2
sz1 = pz1-pz2
magtvec = SQRT(sx1**2+sy1**2+sz1**2)
sx1 = sx1/magtvec
sy1 = sy1/magtvec
sz1 = sz1/magtvec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!! Set 1st tangent vector as the first surface vector
tx1 = sx1
ty1 = sy1
tz1 = sz1
!! Calculate 2nd tangent vector on tangent plane using surface normal direction and one of the surface vector
tx2 = (ny1*tz1-nz1*ty1)
ty2 = -(nx1*tz1-nz1*tx1)
tz2 = (nx1*ty1-ny1*tx1)
magtvec = SQRT(tx2**2+ty2**2+tz2**2)
tx2 = tx2/magtvec
ty2 = ty2/magtvec
tz2 = tz2/magtvec

!! Re-calculate 1st tangent vector on tangent plane using surface normal direction and the other surface vector
tx1 = -(ny1*tz2-nz1*ty2)
ty1 = -(-(nx1*tz2-nz1*tx2))
tz1 = -(nx1*ty2-ny1*tx2)
magtvec = SQRT(tx1**2+ty1**2+tz1**2)
tx1 = tx1/magtvec
ty1 = ty1/magtvec
tz1 = tz1/magtvec


!! We have dphi/dn. Need to get dphi/dt from the fluid
i = ip0
j = jp0
k = kp0
IF ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i-1,j,k))/(x(i+1)-x(i-1))
ELSEIF ((node(i+1,j,k).NE.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i,j,k)-phiTemp(i-1,j,k))/(x(i)-x(i-1))
ELSEIF  ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).NE.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i,j,k))/(x(i+1)-x(i))
ENDIF
IF ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j-1,k))/(y(j+1)-y(j-1))
ELSEIF ((node(i,j+1,k).NE.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j,k)-phiTemp(i,j-1,k))/(y(j)-y(j-1))
ELSEIF  ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).NE.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j,k))/(y(j+1)-y(j))
ENDIF
IF ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k-1))/(z(k+1)-z(k-1))
ELSEIF ((node(i,j,k+1).NE.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k)-phiTemp(i,j,k-1))/(z(k)-z(k-1))
ELSEIF  ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).NE.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k))/(z(k+1)-z(k))
ENDIF

dphidtA1 = dphidx*tx1 + dphidy*ty1 + dphidz*tz1
dphidtA2 = dphidx*tx2 + dphidy*ty2 + dphidz*tz2
!! get dphidtB 
!
i = ip1
j = jp1
k = kp1
IF ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i-1,j,k))/(x(i+1)-x(i-1))
ELSEIF ((node(i+1,j,k).NE.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i,j,k)-phiTemp(i-1,j,k))/(x(i)-x(i-1))
ELSEIF  ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).NE.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i,j,k))/(x(i+1)-x(i))
ENDIF
IF ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j-1,k))/(y(j+1)-y(j-1))
ELSEIF ((node(i,j+1,k).NE.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j,k)-phiTemp(i,j-1,k))/(y(j)-y(j-1))
ELSEIF  ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).NE.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j,k))/(y(j+1)-y(j))
ENDIF
IF ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k-1))/(z(k+1)-z(k-1))
ELSEIF ((node(i,j,k+1).NE.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k)-phiTemp(i,j,k-1))/(z(k)-z(k-1))
ELSEIF  ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).NE.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k))/(z(k+1)-z(k))
ENDIF

dphidtB1 = dphidx*tx1 + dphidy*ty1 + dphidz*tz1
dphidtB2 = dphidx*tx2 + dphidy*ty2 + dphidz*tz2
!! get dphidtC 
!
i = ip2
j = jp2
k = kp2
IF ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i-1,j,k))/(x(i+1)-x(i-1))
ELSEIF ((node(i+1,j,k).NE.FLUID).AND.(node(i-1,j,k).EQ.FLUID)) THEN
	dphidx = (phiTemp(i,j,k)-phiTemp(i-1,j,k))/(x(i)-x(i-1))
ELSEIF  ((node(i+1,j,k).EQ.FLUID).AND.(node(i-1,j,k).NE.FLUID)) THEN
	dphidx = (phiTemp(i+1,j,k)-phiTemp(i,j,k))/(x(i+1)-x(i))
ENDIF
IF ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j-1,k))/(y(j+1)-y(j-1))
ELSEIF ((node(i,j+1,k).NE.FLUID).AND.(node(i,j-1,k).EQ.FLUID)) THEN
	dphidy = (phiTemp(i,j,k)-phiTemp(i,j-1,k))/(y(j)-y(j-1))
ELSEIF  ((node(i,j+1,k).EQ.FLUID).AND.(node(i,j-1,k).NE.FLUID)) THEN
	dphidy = (phiTemp(i,j+1,k)-phiTemp(i,j,k))/(y(j+1)-y(j))
ENDIF
IF ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k-1))/(z(k+1)-z(k-1))
ELSEIF ((node(i,j,k+1).NE.FLUID).AND.(node(i,j,k-1).EQ.FLUID)) THEN
	dphidz = (phiTemp(i,j,k)-phiTemp(i,j,k-1))/(z(k)-z(k-1))
ELSEIF  ((node(i,j,k+1).EQ.FLUID).AND.(node(i,j,k-1).NE.FLUID)) THEN
	dphidz = (phiTemp(i,j,k+1)-phiTemp(i,j,k))/(z(k+1)-z(k))
ENDIF

dphidtC1 = dphidx*tx1 + dphidy*ty1 + dphidz*tz1
dphidtC2 = dphidx*tx2 + dphidy*ty2 + dphidz*tz2

!! get dphidt
dphidtAst1 = 0.5_dbl*(q+1.0_dbl)*(q+2.0_dbl)*dphidtA1+0.5_dbl*(q+1.0_dbl)*(q)*dphidtC1-(q+2.0_dbl)*(q)*dphidtB1
dphidtAst2 = 0.5_dbl*(q+1.0_dbl)*(q+2.0_dbl)*dphidtA2+0.5_dbl*(q+1.0_dbl)*(q)*dphidtC2-(q+2.0_dbl)*(q)*dphidtB2

! We already know dphidn. Now get dphidk.
! get grad(phi) vector
! Set up linear systems to get components of grad(phi) at Astar


!write(*,*) nx1,ny1,nz1,tx1,ty1,tz1,tx2,ty2,tz2 
aa1 = nx1; 
aa2 = ny1; 
aa3 = nz1;
bb1 = tx1; 
bb2 = ty1; 
bb3 = tz1;
cc1 = tx2; 
cc2 = ty2;
cc3 = tz2;
rb1 = dphidn;
rb2 = dphidtAst1;
rb3 = dphidtAst2;

t1 = (cc1*aa2-cc2*aa1)*(cc1*bb3-cc3*bb1)-(cc1*aa3-cc3*aa1)*(cc1*bb2-cc2*bb1)
t2 = (cc1*rb2-cc2*rb1)*(cc1*bb3-cc3*bb1)-(cc1*rb3-cc3*rb1)*(cc1*bb2-cc2*bb1)
sx = t2/(MAX(abs(t1),epsmin)*sign(1.0_dbl,t1))
t1 = cc1*bb2-cc2*bb1
sy = ((cc1*rb2-cc2*rb1)-sx*(cc1*aa2-cc2*aa1))/(MAX(abs(t1),epsmin)*sign(1.0_dbl,t1))
t1 = cc1
sz = (rb1 - aa1*sx - bb1*sy)/(MAX(abs(t1),epsmin)*sign(1.0_dbl,t1))

dphidx = sx
dphidy = sy
dphidz = sz


!! get dphidk = gradphi.dot.k_vec
dist = SQRT(ex(bb(m))**2+ey(bb(m))**2+ez(bb(m))**2)
dx=ex(bb(m))/dist	! i direction
dy=ey(bb(m))/dist	! j direction
dz=ez(bb(m))/dist	! k direction
dphidk = (dphidx*dx+dphidy*dy+dphidz*dz)
!! now get phiAst
phiA = phiTemp(ip0,jp0,kp0)
phiB = phiTemp(ip1,jp1,kp1)
phiC = phiTemp(ip2,jp2,kp2)
delta = SQRT(ex(bb(m))**2 + ex(bb(m))**2 + ex(bb(m))**2)*xcf
phiAst = ((q+1.0_dbl)*(q+1.0_dbl)*phiA-q*q*phiB-dphidk*delta*q*(q+1.0_dbl))/(2.0_dbl*q+1.0_dbl)
!------------------------------------------------
END SUBROUTINE GetPhiWall
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE GetPerturbedVector(dx,dy,dz,epsx,epsy,epsz)	! Perturb a vector with components dx,dy and dz within a cone (using method from stackoverflow, Ray Tracing News)
!--------------------------------------------------------------------------------------------------
REAL(dbl), INTENT(IN) :: dx,dy,dz  ! direction vector components
REAL(dbl), INTENT(OUT) :: epsx,epsy,epsz ! deviation vector components
REAL(dbl) :: eps,epsmin  ! local variables
REAL(dbl) :: halfangle  ! cone halfangle
REAL(dbl) :: term1,term2,term3  ! temp variables
REAL(dbl) :: randx1,randy1,randz1 ! 1st random unit vector on tangent plane to (dx,dy,dz) vector
REAL(dbl) :: randx2,randy2,randz2 ! 2nd random unit vector on tangent plane to (dx,dy,dz) vector
REAL(dbl) :: s,rv,h,dist,theta,phiang,sinT,xc,yc,zc,newdx,newdy,newdz,randnum ! temp variables

halfangle = 3.0_dbl ! cone half angle in degrees
epsmin = 1.0e-8_dbl

! Algorithm to calculate perturbed vector within a cone

!dist = SQRT(dx**2+dy**2+dz**2)
!dx=dx/dist	! i direction
!dy=dy/dist	! j direction
!dz=dz/dist	! k direction


CALL RANDOM_NUMBER(randnum)
randx1 = 0.1_dbl*(0.0_dbl+randnum)+0.0_dbl
CALL RANDOM_NUMBER(randnum)
randy1 = 0.1_dbl*(0.0_dbl+randnum)+0.0_dbl
CALL RANDOM_NUMBER(randnum)
randz1 = 0.1_dbl*(0.0_dbl+randnum)+0.0_dbl
!write(*,*) randx1,randy1,randz1
if (abs(dx).GT.1.0e-8_dbl) then
randx1 = -(dz*randz1+dy*randy1)/(dx+epsmin)
elseif (abs(dy).GT.1.0e-8_dbl) then
randy1 = -(dz*randz1+dx*randx1)/(dy+epsmin)
else
randz1 = -(dx*randx1+dy*randy1)/(dz+epsmin)
endif
dist = SQRT(randx1**2+randy1**2+randz1**2)
randx1=randx1/dist	! i direction
randy1=randy1/dist	! j direction
randz1=randz1/dist	! k direction

randx2 = (dy*randz1-dz*randy1)
randy2 = -(dx*randz1-dz*randx1)
randz2 = (dx*randy1-dy*randx1)

dist = SQRT(randx2**2+randy2**2+randz2**2)
randx2=randx2/dist	! i direction
randy2=randy2/dist	! j direction
randz2=randz2/dist	! k direction

CALL RANDOM_NUMBER(s)
CALL RANDOM_NUMBER(rv)
s = s*1.0_dbl + 0.0_dbl
rv = rv + 0.0_dbl
theta = halfangle*PI/(180.0_dbl) ! create a cone with coning angle 10 degrees
h = cos(theta)
phiang = 2.0_dbl*PI*s
zc = h + ( 1.0_dbl - h ) * rv
sinT = sqrt( 1.0_dbl - zc * zc )
xc = cos( phiang ) * sinT
yc = sin( phiang ) * sinT
! Estimate the perturbed vector
!perturbed = rand_vector * xc + cross_vector * yc + orig_vector * zc
newdx = randx1*xc + randx2*yc + dx*zc
newdy = randy1*xc + randy2*yc + dy*zc
newdz = randz1*xc + randz2*yc + dz*zc

dist = SQRT(newdx**2+newdy**2+newdz**2)
newdx = newdx/dist	! i direction
newdy = newdy/dist	! j direction
newdz = newdz/dist	! k direction


! We return the perturbations instead of the perturbed vector
epsx = newdx - dx
epsy = newdy - dy
epsz = newdz - dz

!------------------------------------------------
END SUBROUTINE GetPerturbedVector
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE RotateVector(dx,dy,dz,ang1,ang2,ang3)	! Rotate a vector with components dx,dy and dz by angles ang1,ang2,ang3 (in radians)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

REAL(dbl), INTENT(INOUT) :: dx,dy,dz	! original vector
REAL(dbl) :: tempdx,tempdy,tempdz	! teporary vector components
REAL(dbl), INTENT(IN) :: ang1,ang2,ang3	! rotation angles to rotate (perturb) the vector joining fluid and solid nodes for identifying nearby nodes on the solid surface

! ROtate a vector using a 3D ortation matrix. I directly compute the matrix multiplied vector

tempdx = dx*(cos(ang2)*cos(ang1)) + dy*(cos(ang3)*sin(ang1)+sin(ang3)*sin(ang2)*cos(ang1))  &
		+ dz*(sin(ang3)*sin(ang1)-cos(ang3)*sin(ang2)*cos(ang1))
tempdy = dx*(-cos(ang2)*sin(ang1)) + dy*(cos(ang3)*cos(ang1)-sin(ang3)*sin(ang2)*sin(ang1)) &
		+ dz*(sin(ang3)*cos(ang1)+cos(ang3)*sin(ang2)*sin(ang1))
tempdz = dx*(sin(ang2)) + dy*(-sin(ang3)*cos(ang2)) + dz*(cos(ang3)*cos(ang2))

dx=tempdx
dy=tempdy
dz=tempdz

!------------------------------------------------
END SUBROUTINE RotateVector
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE RotateVectorAboutAxis(dx,dy,dz,ang1,ang2,ang3)	! Rotate a vector with components dx,dy and dz by angles ang1,ang2,ang3 (in radians)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

REAL(dbl), INTENT(INOUT) :: dx,dy,dz	! original vector
REAL(dbl) :: tempdx,tempdy,tempdz	! teporary vector components
REAL(dbl), INTENT(IN) :: ang1,ang2,ang3	! rotation angles to rotate (perturb) the vector joining fluid and solid nodes for identifying nearby nodes on the solid surface

! ROtate a vector using a 3D ortation matrix. I directly compute the matrix multiplied vector

tempdx = dx*(cos(ang1)+dx*dx*(1.0_dbl-cos(ang1))) + dy*(dx*dy*(1.0_dbl-cos(ang1))-dz*sin(ang1))  &
		+ dz*(dx*dz*(1.0_dbl-cos(ang1))+dy*sin(ang1))
tempdy = dx*(dy*dx*(1.0_dbl-cos(ang1))+dz*sin(ang1)) + dy*(cos(ang1)+dy*dy*(1.0_dbl-cos(ang1))) &
		+ dz*(dy*dz*(1.0_dbl-cos(ang1))-dx*sin(ang1))
tempdz = dx*(dz*dx*(1.0_dbl-cos(ang1))-dy*sin(ang1)) + dy*(dz*dy*(1.0_dbl-cos(ang1))+dx*sin(ang1)) &
		+ dz*(cos(ang1)+dz*dz*(1.0_dbl-cos(ang1)))

dx=tempdx
dy=tempdy
dz=tempdz

!------------------------------------------------
END SUBROUTINE RotateVectorAboutAxis
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE IntersectSurface(m,i,j,k,im1,jm1,km1,px,py,pz,q,ang1,ang2,ang3)			! Intersects boudnary surface using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
REAL(dbl), INTENT(INOUT) :: ang1,ang2,ang3	! rotation angles to rotate (perturb) the vector joining fluid and solid nodes for identifying nearby nodes on the solid surface
REAL(dbl), INTENT(OUT) :: px,py,pz			! boundary points
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Bx,By,Bz									! solid node
REAL(dbl) :: AB,AP1,AP2,AP										! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: dxorig,dyorig,dzorig				! unit vector pointing from A to B
REAL(dbl) :: r1,r2,z1,z2,slope,intercept			! radius and z location at k and km1, slope of line connecting those two points, z-intercept of the r-equation
REAL(dbl) :: slope2,term1,term2						! terms used in calculation
REAL(dbl), INTENT(OUT) :: q						! distance from fluid point to surface using Ray Tracing


!open(9,file='testoutput.dat',position='append')

! RAY
! point A (current node)
Ax = x(i)/xcf
Ay = y(j)/ycf
Az = z(k)/zcf

! point B (solid node)
Bx = x(im1)/xcf
By = y(jm1)/ycf
Bz = z(km1)/zcf

! distance from A to B
AB = SQRT((Bx - Ax)*(Bx - Ax) + (By - Ay)*(By - Ay) + (Bz - Az)*(Bz - Az))

! unit vector (d) from point A to point B
dxorig = (x(im1)-x(i))/AB									! i direction
dyorig = (y(jm1)-y(j))/AB									! j direction
dzorig = (z(km1)-z(k))/AB									! k direction

10 dx=dxorig+ang3
dy=dyorig+ang2
dz=dzorig+ang1

AB = SQRT(dx**2+dy**2+dz**2)
dx=dx/AB	! i direction
dy=dy/AB	! j direction
dz=dz/AB	! k direction

AB = SQRT((Bx - Ax)*(Bx - Ax) + (By - Ay)*(By - Ay) + (Bz - Az)*(Bz - Az))

! SURFACE
IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	r1 = rOut(k)													! radius at k (distance from CL)
	r2 = rOut(km1)													! radius at km1 (distance from CL)
ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
	r1 = rIn(k)													! radius at k (distance from CL)
	r2 = rIn(km1)													! radius at km1 (distance from CL)
END IF
z1 = z(k)/zcf													! z-coordinate at k
z2 = z(km1)/zcf													! z-coordinate at km1

IF(k .NE. km1) THEN
  slope = (r2-r1)/(z2-z1)								! approximate the surface as a conincal shell (linear between k values)
ELSE
  slope = 0.0_dbl
END IF

intercept = r1 - slope*z1								! z-intercept of the linearly approximated r-equation

! terms used in calculation
slope2 = slope*slope															! slope^2
term1 = Ax*dx + Ay*dy - Az*dz*slope2 - intercept*dz*slope	! reoccuring term
term2 = dx*dx + dy*dy - dz*dz*slope*slope			! reoccuring term

! calculate the distance from the current node (point A) to the wall (point P)
AP1 = (1.0_dbl/(2.0_dbl*term2)) * &
     (-2.0_dbl*term1					&
   + SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))
AP2 = (1.0_dbl/(2.0_dbl*term2)) * &
     (-2.0_dbl*term1					&
   - SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))
IF ((AP1.GE.0.0_dbl).AND.(AP2.GE.0.0_dbl)) THEN
 AP = MIN(AP1,AP2)
ELSE
 AP = MAX(AP1,AP2)
ENDIF
q = AP/AB	! distance ratio

! calculate coordinates of the boudnayr interseciton point
px = Ax + AP*dx
py = Ay + AP*dy
pz = Az + AP*dz


! balaji added
!q=0.5

! make sure 0<q<1
IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 2.00000001_dbl)) THEN 
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) " ERROR in IntersectSurface"
  WRITE(1000,*) "q=",q
  WRITE(1000,*) "m=",m
  WRITE(1000,*) "i=",i,"j=",j,"k=",k
  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
  WRITE(1000,*) "Ax=",Ax,"Ay=",Ay,"Az=",Az
  WRITE(1000,*) "Bx=",Bx,"By=",By,"Bz=",Bz
  WRITE(1000,*) "Px=",px,"Py=",py,"Pz=",pz
  WRITE(1000,*) "dxorig=",dxorig,"dyorig=",dyorig,"dzorig=",dzorig
  WRITE(1000,*) "dx=",dx,"dy=",dy,"dz=",dz
  WRITE(1000,*) "epsx=",ang3,"epsy=",ang2,"epsz=",ang1
  WRITE(1000,*) "r1=",r1,"r2=",r2
  WRITE(1000,*) "z1=",z1,"z2=",z2
  WRITE(1000,*) "slope=",slope
  WRITE(1000,*) "term1=",term1,"term2=",term2
  WRITE(1000,*) "intercept=",intercept
  WRITE(1000,*) "AB=",AB,"AP=",AP
  CLOSE(1000)
  STOP
END IF																																
!CLOSE(9)

!------------------------------------------------
END SUBROUTINE IntersectSurface
!------------------------------------------------










! -------------------------------------------------------------------- 
SUBROUTINE Gauss (a,ainv,n)       ! Invert matrix by Gauss method 
! -------------------------------------------------------------------- 
IMPLICIT NONE 
INTEGER :: n 
REAL(dbl), INTENT(IN) :: a(n,n) 
REAL(dbl), INTENT(OUT) :: ainv(n,n) 
! - - - Local Variables - - - 
REAL(dbl) :: b(n,n), c, d, temp(n) 
INTEGER(lng) :: i, j, k, m, imx(1), ipvt(n) 
! - - - - - - - - - - - - - - 
b = a 
ipvt = (/ (i, i = 1, n) /) 
DO k = 1,n 
   imx = MAXLOC(ABS(b(k:n,k))) 
   m = k-1+imx(1) 
   IF (m /= k) THEN 
      ipvt( (/m,k/) ) = ipvt( (/k,m/) ) 
      b((/m,k/),:) = b((/k,m/),:) 
   END IF 
   d = 1/b(k,k) 
   temp = b(:,k) 
   DO j = 1, n 
      c = b(k,j)*d 
      b(:,j) = b(:,j)-temp*c 
      b(k,j) = c 
   END DO 
   b(:,k) = temp*(-d) 
   b(k,k) = d 
END DO 
ainv(:,ipvt) = b 
!--------------------------------------------------------------------------------------------------
END SUBROUTINE Gauss
!--------------------------------------------------------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE ScalarBCV(m,i,j,k,im1,jm1,km1,vNum,phiBC)						! implements the scalar BCs for the villi
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1								! index variables
INTEGER(lng), INTENT(IN) :: vNum													! number of the current villus
REAL(dbl), INTENT(OUT) :: phiBC     											! scalar contribution from the boundary condition
INTEGER(dbl) :: ip1,jp1,kp1 														! neighboring nodes (2 away from the wall)
REAL(dbl) :: Cx,Cy,Cz,Vx,Vy,Vz													! vector between villous base and current node, vector between villous base and villous tip
REAL(dbl) :: uV,vV,wV																! villus wall velocity at actual coordinate system
REAL(dbl) :: q																			! distance ratio from the current node to the solid node
REAL(dbl) :: rhoB,phiB																! values of density and at the boundary, and contribution of scalar from the boundary and solid nodes
REAL(dbl) :: feq_m																	! equilibrium distribution function in the mth direction
REAL(dbl) :: phiijk_m																! contribution of scalar streamed in the mth direction to (ip1,jp1,kp1)

! find C and V
CALL qCalcV(m,i,j,k,im1,jm1,km1,vNum,q,Cx,Cy,Cz,Vx,Vy,Vz,2)

! find the influence of villous velocity on the current point
CALL VilliVelocity(vNum,Cx,Cy,Cz,uV,vV,wV)

! neighboring node (fluid side)	
ip1 = i + ex(m) 																		! i + 1
jp1 = j + ey(m)																		! j + 1
kp1 = k + ez(m)																		! k + 1		

! if (ip1,jp1,kp1) is not in the fluid domain, use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF	

! assign values to boundary (density, scalar, f)
rhoB = (rho(i,j,k) - rho(ip1,jp1,kp1))*(1+q) + rho(ip1,jp1,kp1)		! extrapolate the density
CALL Equilibrium_LOCAL(m,rhoB,uV,vV,wV,feq_m)								! calculate the equibrium distribution function in the mth direction

! find the contribution of scalar streamed from the wall to the current node (i,j,k), and from the current node to the next neighboring node (ip1,jp1,kp1)
phiB		= (feq_m/rhoB - wt(m)*Delta)*phiWall								! contribution from the wall in the mth direction (zero if phiWall=0)
phiijk_m	= (fplus(m,i,j,k)/rho(i,j,k) - wt(m)*Delta)*phiTemp(i,j,k)	! contribution from the current node to the next node in the mth direction

! if q is too small, the extrapolation to phiBC can create a large error...
IF(q .LT. 0.25) THEN
  q = 0.25_dbl  																		! approximate the distance ratio as 0.25
END IF

! extrapolate using phiB and phijk_m to obtain contribution from the solid node to the current node
phiBC		= ((phiB - phiijk_m)/q) + phiB										! extrapolated scalar value at the solid node, using q

!------------------------------------------------
END SUBROUTINE ScalarBCV
!------------------------------------------------










!================================================
END MODULE BCscalar
!================================================
