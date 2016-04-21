!==================================================================================================
MODULE BClbm		! Sets LBM Boundary Conditions
!==================================================================================================
USE SetPrecision
USE Setup  
USE MPI
USE IC

CONTAINS





!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			! implements the (moving) bounceback boundary conditions (1st order accurate - Ladd)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb											! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk													! radius of current node

rijk = x(im1)								! height at current location

IF (rijk .GE. rOut(k)) THEN
	ub = velOut(km1)	 ! 0.0_dbl							! x-component of the velocity at i,j,k
	vb = 0.0_dbl		 !0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
	wb = 0.0_dbl 		 !velOut(km1)	!vel(km1)!0.0_dbl						! only z-component in this case			
ELSE IF (rijk .LE. rIn(k)) THEN
	ub = velIn(km1)		 !0.0_dbl!0.0!vel(km1)*cosTheta						! x-component of the velocity at i,j,k
	vb = 0.0_dbl		 !0.0!vel(km1)*sinTheta						! y-component of the velocity at i,j,k
	wb = 0.0_dbl 		 !velIn(km1)!vel(km1)!0.0_dbl						! only z-component in this case	
END IF				

fbb = fplus(bb(m),i,j,k) + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m))	! bounced back distribution function with added momentum
fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
!------------------------------------------------
END SUBROUTINE BounceBackL
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBack2(m,i,j,k,im1,jm1,km1,fbb)	! implements the (moving) bounceback boundary conditions (2nd order accurate - Lallemand)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2					! index variables
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb											! wall velocity (x-, y-, z- components)
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]
REAL(dbl) :: rijk													! radius of current node

ip1 = i + ex(m)													! i location of 1st neighbor in the m direction
jp1 = j + ey(m)													! j location of 1st neighbor in the m direction
kp1 = k + ez(m)													! k location of 1st neighbor in the m direction

ip2 = i + 2_lng*ex(m)											! i location of 2nd neighbor in the m direction
jp2 = j + 2_lng*ey(m)											! j location of 2nd neighbor in the m direction
kp2 = k + 2_lng*ez(m)											! k location of 2nd neighbor in the m direction


IF((node(ip1,jp1,kp1) .EQ. FLUID) .AND. (node(ip2,jp2,kp2) .EQ. FLUID)) THEN		! continue with 2nd order BB if the two positive neighbors are in the fluid (most cases)

  rijk = x(im1)										! height at current location
	IF (rijk .GE. rOut(k)) THEN
		ub = velOut(km1) 	!0.0_dbl 						! x-component of the velocity at i,j,k
		vb = 0.0_dbl									! y-component of the velocity at i,j,k
		wb = 0.0_dbl		 !  velOut(km1)!vel(km1)!0.0_dbl						! only z-component in this case			
	ELSE IF (rijk .LE. rIn(k)) THEN
		ub = velIn(km1) !0.0_dbl 				! x-component of the velocity at i,j,k
		vb = 0.0_dbl 						! y-component of the velocity at i,j,k
		wb = 0.0_dbl 	!velIn(km1) 				! only z-component in this case	
	END IF	

  CALL qCalc(m,i,j,k,im1,jm1,km1,q)							! calculate q					

  ! bounced back distribution function with added momentum
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
    fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 															&
        + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) 													& 
        - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 													&
        !+ 6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m))
        + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)) ! Set rho = 1.0
        !+ 6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m))) ! use actual rho that fluctates

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
    fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q																&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)							&
        !+ (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
        + (6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))	! Set rho = 1.0
        !+ (6.0_dbl*wt(bb(m))*rho(i,j,k)*(ub*ex(bb(m)) + vb*ey(bb(m)) + wb*ez(bb(m))))/(q*(2.0_dbl*q + 1.0_dbl)) ! Use actual rho that fluctuateѕ

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
  ELSE
    OPEN(1000,FILE='error-'//sub//'.txt')
    WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
    CLOSE(1000)
    STOP
  END IF

ELSE

  CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)

END IF

!------------------------------------------------
END SUBROUTINE BounceBack2
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBack2New(m,i,j,k,im1,jm1,km1,fbb)	! implements the (moving) bounceback boundary conditions (2nd order accurate - Lallemand)
! Implemented by Balaji 10/28/2014 using a method similar to Yanxing
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2					! index variables
REAL(dbl) :: cosTheta, sinTheta								! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb											! wall velocity (x-, y-, z- components)
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]
REAL(dbl) :: rijk													! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt				! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 
INTEGER(lng) :: it			! loop index variables

ip1 = i + ex(m)													! i location of 1st neighbor in the m direction
jp1 = j + ey(m)													! j location of 1st neighbor in the m direction
kp1 = k + ez(m)													! k location of 1st neighbor in the m direction

ip2 = i + 2_lng*ex(m)											! i location of 2nd neighbor in the m direction
jp2 = j + 2_lng*ey(m)											! j location of 2nd neighbor in the m direction
kp2 = k + 2_lng*ez(m)											! k location of 2nd neighbor in the m direction


IF((node(ip1,jp1,kp1) .EQ. FLUID) .AND. (node(ip2,jp2,kp2) .EQ. FLUID)) THEN		! continue with 2nd order BB if the two positive neighbors are in the fluid (most cases)

!*****************************************************************************
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

      		   rt = xt
                   IF (rijk .GE. rOut(k)) THEN
		   	ht = ((zt-z(k))*rOut(km1)+(z(km1)-zt)*rOut(k))/(z(km1)-z(k))
		   ELSE
		   	ht = ((zt-z(k))*rIn(km1)+(z(km1)-zt)*rIn(k))/(z(km1)-z(k))
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

      		   rt = xt
		   IF (rijk .GE. rOut(k)) THEN
			   ht = (rOut(km1)+rOut(k))/2.0_dbl
		   ELSE
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
		   IF (rijk .GE. rOut(k)) THEN
			 vt = ((zt-z(k))*velOut(km1)+(z(km1)-zt)*velOut(k))/(z(km1)-z(k))
		   ELSE
			 vt = ((zt-z(k))*velIn(km1)+(z(km1)-zt)*velIn(k))/(z(km1)-z(k))
		   END IF
	 ELSE
		   IF (rijk .GE. rOut(k)) THEN
			 vt = (velOut(k)+velOut(km1))*0.5_dbl
		   ELSE
			 vt = (velIn(k)+velIn(km1))*0.5_dbl
		   END IF
	 ENDIF
		ub = velIn(km1)	!0.0_dbl!0.0!vel(km1)*cosTheta								! x-component of the velocity at i,j,k
		vb = 0.0_dbl	!0.0!vel(km1)*sinTheta								! y-component of the velocity at i,j,k
		wb = 0.0_dbl	 !vt!vel(km1)!0.0_dbl									! only z-component in this case)
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
	!q = max(q, 0.25_dbl)
!*****************************************************************************
!*****************************************************************************
!! Original method used by Gino
!  rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))				! radius at current location
!
!  cosTheta = x(im1)/rijk										! COS(theta)
!  sinTheta = y(jm1)/rijk										! SIN(theta)
!  !cosTheta = x(im1)/r(km1)									! COS(theta)
!  !sinTheta = y(jm1)/r(km1)									! SIN(theta)
!
!  ub = vel(km1)*cosTheta										! x-component of the velocity at i,j,k
!  vb = vel(km1)*sinTheta										! y-component of the velocity at i,j,k
!  wb = 0.0_dbl														! no z-component in this case)
!
!  CALL qCalc(m,i,j,k,im1,jm1,km1,q)							! calculate q					
!*****************************************************************************
		 !write(*,*) 'q-Gino',q,i,j,k,im1,jm1,km1


  ! bounced back distribution function with added momentum
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
    fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 															&
        + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) 													& 
        - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 													&
        !+ 6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)) ! use actual rho
        + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)) ! use rho = 1.0

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
    fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q																&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)							&
        !+ (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))	! Use actual rho
        + (6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))		! Use rho = 1.0

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))

  ELSE
    OPEN(1000,FILE='error-'//sub//'.txt')
    WRITE(1000,*) "Error in BounceBack2() in ICBC.f90 (line 137): q is not (0<=q<=1)...? Aborting."
    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
    CLOSE(1000)
    STOP
  END IF

ELSE

  CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)

END IF

!------------------------------------------------
END SUBROUTINE BounceBack2New
!------------------------------------------------







!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalc(m,i,j,k,im1,jm1,km1,q)			! calculates q (boundary distance ratio) using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
REAL(dbl), INTENT(OUT) :: q							! distance ratio
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Bx,By,Bz									! solid node
REAL(dbl) :: AB,AP										! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: r1,r2,z1,z2,slope,intercept			! radius and z location at k and km1, slope of line connecting those two points, z-intercept of the r-equation
REAL(dbl) :: slope2,term1,term2						! terms used in calculation
REAL(dbl) :: h1,h2, time ! Height of the lower and upper pistons and the physical time


!Get location of upper and lower pistons

time = iter*tcf

!h1 = lowerPistonVel * time		 !Ganesh
!h2 = upperPistonVel * time + initHeight !Ganesh


h2 = -0.0 *D + s1*time 
h1 = -0.4 *D + s1*time
 

! RAY
! point A (current node)
Ax = x(i)
Ay = y(j)
Az = z(k)

! point B (solid node)
Bx = x(im1)
By = y(jm1)
Bz = z(km1)

if (Bx .ge. h2) then

   q = (h2-Ax)/(Bx-Ax)
   
else if(Bx .le. h1) then

   q = (Ax-h1)/(Ax-Bx)

end if

! make sure 0<q<1
IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) "q    =",q
 WRITE(1000,*) "h1,h2 =",h1,h2
  WRITE(1000,*) "m=",m
  WRITE(1000,*) "i=",i,"j=",j,"k=",k
  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
  WRITE(1000,*) "Ax=",Ax,"Ay=",Ay,"Az=",Az
  WRITE(1000,*) "Bx=",Bx,"By=",By,"Bz=",Bz
  WRITE(1000,*) "dx=",dx,"dy=",dy,"dz=",dz
  WRITE(1000,*) "r1=",r1,"r2=",r2
  WRITE(1000,*) "z1=",z1,"z2=",z2
  WRITE(1000,*) "slope=",slope
  WRITE(1000,*) "term1=",term1,"term2=",term2
  WRITE(1000,*) "intercept=",intercept
  WRITE(1000,*) "AB=",AB,"AP=",AP
  CLOSE(1000)
  STOP
END IF																																

!------------------------------------------------
END SUBROUTINE qCalc
!------------------------------------------------




!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalcBalaji(m,i,j,k,im1,jm1,km1,q)			! calculates q (boundary distance ratio) using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
REAL(dbl), INTENT(OUT) :: q							! distance ratio
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Bx,By,Bz									! solid node
REAL(dbl) :: AB,AP										! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: r1,r2,z1,z2,slope,intercept			! radius and z location at k and km1, slope of line connecting those two points, z-intercept of the r-equation
REAL(dbl) :: slope2,term1,term2						! terms used in calculation

! RAY
! point A (current node)
Ax = x(i)
Ay = y(j)
Az = z(k)

! point B (solid node)
Bx = x(im1)
By = y(jm1)
Bz = z(km1)

! distance from A to B
AB = SQRT((Bx - Ax)*(Bx - Ax) + (By - Ay)*(By - Ay) + (Bz - Az)*(Bz - Az))

! unit vector (d) from point A to point B
dx = (x(im1)-x(i))/AB									! i direction
dy = (y(jm1)-y(j))/AB									! j direction
dz = (z(km1)-z(k))/AB									! k direction

! SURFACE
IF (node(im1,jm1,km1).EQ.SOLID2) THEN
	r1 = rOut(k)													! radius at k (distance from CL)
	r2 = rOut(km1)													! radius at km1 (distance from CL)
ELSE IF (node(im1,jm1,km1).EQ.SOLID) THEN
	r1 = rIn(k)													! radius at k (distance from CL)
	r2 = rIn(km1)													! radius at km1 (distance from CL)
END IF
z1 = z(k)													! z-coordinate at k
z2 = z(km1)													! z-coordinate at km1

IF(k .NE. km1) THEN
  slope = (r2-r1)/(z2-z1)								! approximate the surface as a conincal shell (linear between k values)
ELSE
  slope = 0.0_dbl
END IF

intercept = r1 - slope*z1								! z-intercept of the linearly approximated r-equation

! terms used in calculation
slope2 = slope*slope															! slope^2
term1 = Ax*dx + Ay*dy - Az*dz*slope2 - intercept*dz*slope		! reoccuring term
term2 = dx*dx + dy*dy - dz*dz*slope*slope								! reoccuring term

! calculate the distance from the current node (point A) to the wall (point P)
AP = (1.0_dbl/(2.0_dbl*term2)) * &
     (-2.0_dbl*term1					&
   + SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))

q = AP/AB													! distance ratio
q = max(q,0.001_dbl)

! balaji added
!q=0.5

! make sure 0<q<1
IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) "Farhad"
  WRITE(1000,*) "q=",q
  WRITE(1000,*) "m=",m
  WRITE(1000,*) "i=",i,"j=",j,"k=",k
  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
  WRITE(1000,*) "Ax=",Ax,"Ay=",Ay,"Az=",Az
  WRITE(1000,*) "Bx=",Bx,"By=",By,"Bz=",Bz
  WRITE(1000,*) "dx=",dx,"dy=",dy,"dz=",dz
  WRITE(1000,*) "r1=",r1,"r2=",r2
  WRITE(1000,*) "z1=",z1,"z2=",z2
  WRITE(1000,*) "slope=",slope
  WRITE(1000,*) "term1=",term1,"term2=",term2
  WRITE(1000,*) "intercept=",intercept
  WRITE(1000,*) "AB=",AB,"AP=",AP
  CLOSE(1000)
  STOP
END IF																																

!------------------------------------------------
END SUBROUTINE qCalcBalaji
!------------------------------------------------




!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalcFarhad(i,q)  
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i
REAL(dbl)    :: h1,h2,time,D_X,D_Y,q

time = iter*tcf
D_X= 20*D 
D_Y= D
h2=(-0.38* D_x)+ (s1*time)    
h1=(-0.48* D_x)+ (s1*time) 

IF (x(i) .LT. 0.5*(h1+h2) ) then 		  			!Left    
   q= (x(i)-h1)/xcf
ELSE						  	 	!Right 
   q= (h2-x(i))/xcf
END IF

IF ((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) "Farhad q", q
  CLOSE(1000)
  STOP
END IF
!------------------------------------------------
END SUBROUTINE qCalcFarhad
!------------------------------------------------






!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBackVL(m,i,j,k,im1,jm1,km1,vNum,fbb)	! implements the (moving) bounceback boundary conditions (1st order accurate - Ladd)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
INTEGER(lng), INTENT(IN) :: vNum								! number of the current villus
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
REAL(dbl) :: Cx,Cy,Cz,Vx,Vy,Vz								! vector between villous base and current node, vector between villous base and villous tip
REAL(dbl) :: uV,vV,wV											! villus wall velocity at actual coordinate system
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]

! find C and V
CALL qCalcV(m,i,j,k,im1,jm1,km1,vNum,q,Cx,Cy,Cz,Vx,Vy,Vz,1)

! find the influence of villous velocity on the current point
CALL VilliVelocity(vNum,Cx,Cy,Cz,uV,vV,wV)

! calculate bounced back distribution function with added momentum
fbb = fplus(bb(m),i,j,k) + 6.0_dbl*wt(m)*rho(i,j,k)*(uV*ex(m) + vV*ey(m) + wV*ez(m))	

!------------------------------------------------
END SUBROUTINE BounceBackVL
!------------------------------------------------











!--------------------------------------------------------------------------------------------------
SUBROUTINE BounceBackV2(m,i,j,k,im1,jm1,km1,vNum,fbb)	! implements the (moving) bounceback boundary conditions (2nd order accurate - Lallemand)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
INTEGER(lng), INTENT(IN) :: vNum								! number of the current villus
REAL(dbl), INTENT(OUT) :: fbb									! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2					! index variables
REAL(dbl) :: Cx,Cy,Cz,Vx,Vy,Vz								! vector between villous base and current node, vector between villous base and villous tip
REAL(dbl) :: uV,vV,wV											! villus wall velocity at actual coordinate system
REAL(dbl) :: q														! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]

! neighboring nodes
ip1 = i + ex(m)													! i location of 1st neighbor in the m direction
jp1 = j + ey(m)													! j location of 1st neighbor in the m direction
kp1 = k + ez(m)													! k location of 1st neighbor in the m direction

ip2 = i + 2_lng*ex(m)											! i location of 2nd neighbor in the m direction
jp2 = j + 2_lng*ey(m)											! j location of 2nd neighbor in the m direction
kp2 = k + 2_lng*ez(m)											! k location of 2nd neighbor in the m direction

IF((node(ip1,jp1,kp1) .EQ. FLUID) .AND. (node(ip2,jp2,kp2) .EQ. FLUID)) THEN		! continue with 2nd order BB if the two positive neighbors are in the fluid (most cases)

  ! find q (and C and V)
  CALL qCalcV(m,i,j,k,im1,jm1,km1,vNum,q,Cx,Cy,Cz,Vx,Vy,Vz,1)

  ! find the influence of villous velocity on the current point
  CALL VilliVelocity(vNum,Cx,Cy,Cz,uV,vV,wV)

  ! bounced back distribution function with added momentum
  IF((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
    fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 															&
        + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) 													& 
        - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 													&
        + 6.0_dbl*wt(m)*rho(i,j,k)*(uV*ex(m) + vV*ey(m) + wV*ez(m))
  ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
    fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 														&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q																&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)							&
        + (6.0_dbl*wt(m)*rho(i,j,k)*(uV*ex(m) + vV*ey(m) + wV*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
  ELSE
    OPEN(1000,FILE='error-'//sub//'.txt')
    WRITE(1000,*) "Error in BounceBackV2() in ICBC.f90 (line 516): q is not (0<=q<=1)...? Aborting."
    WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
    CLOSE(1000)
    STOP
  END IF

ELSE			! must be a case close to two boundaries (villi and surface) - use 1st order BB

  CALL BounceBackVL(m,i,j,k,im1,jm1,km1,vNum,fbb)

END IF

!------------------------------------------------
END SUBROUTINE BounceBackV2
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE VilliVelocity(vNum,Cx,Cy,Cz,uV,vV,wV)					! calculates the influence of the velocity of the villi on point (i,j,k)
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: vNum											! number of the current villus
REAL(dbl), INTENT(IN) :: Cx,Cy,Cz										! vector between villous base and current node
REAL(dbl), INTENT(OUT) :: uV,vV,wV										! villus velocity
REAL(dbl) :: ubx,uby,rijk													! wall velocity, radius at current villous base location
REAL(dbl) :: omegaX, omegaY, omegaZ, omegaT							! angular velocity components - x,y,z, and azimuthal directions
REAL(dbl) :: alpha															! angle made between the x-axis and the base of the villus													

! find the angular velocity components
omegaT	= (villiLoc(vNum,5) - villiLoc(vNum,9))/tcf				! angular velocity in the aziumthal direction
alpha		= ATAN(villiLoc(vNum,2)/villiLoc(vNum,1))					! angle made between the x-axis and the base of the villus
omegaX	= -omegaT*SIN(alpha)												! angular velocity in the x-direction
omegaY	= omegaT*COS(alpha)												! angular velocity in the y-direction
omegaZ	= (villiLoc(vNum,4) - villiLoc(vNum,10))/tcf				! angular velocity in the z-direction

! find the villous velocity components from the cross product of omega and C
uV = (omegaY*Cz - omegaZ*Cy)/vcf											! x-component
vV = (omegaZ*Cx - omegaX*Cz)/vcf											! y-component
wV = (omegaX*Cy - omegaY*Cx)/vcf											! z-component

! add the components from the wall to the u and v velocity
!rijk = SQRT(villiLoc(vNum,1)**2 + villiLoc(vNum,2)**2)						! radius at the villous base
!ubx = (velDom(NINT(villiLoc(vNum,3)/zcf))/vcf)*(villiLoc(vNum,1)/rijk)	! x-component of the wall velocity
!uby = (velDom(NINT(villiLoc(vNum,3)/zcf))/vcf)*(villiLoc(vNum,2)/rijk)	! y-component of the wall velocity

uV = uV + (velDomOut(NINT(villiLoc(vNum,3)/zcf))/vcf)*COS(alpha)	! active villous velocity component + wall velocity component (x-direction)
vV = vV + (velDomOut(NINT(villiLoc(vNum,3)/zcf))/vcf)*SIN(alpha)	! active villous velocity component + wall velocity component (y-direction)

!uV = uV + ubx																	! active villous velocity component + wall velocity component (x-direction)
!vV = vV + uby																	! active villous velocity component + wall velocity component (y-direction)

!------------------------------------------------
END SUBROUTINE VilliVelocity
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalcV(m,i,j,k,im1,jm1,km1,vNum,q,Cx,Cy,Cz,Vx,Vy,Vz,scalarORfluid)	! calculates q (boundary distance ratio) using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
INTEGER(lng), INTENT(IN) :: vNum						! number of the current villus
REAL(dbl), INTENT(OUT) :: q,Cx,Cy,Cz,Vx,Vy,Vz	! distance ratio
REAL(dbl) :: Ax,Ay,Az									! current node
REAL(dbl) :: Bx,By,Bz									! solid node
REAL(dbl) :: ABx,ABy,ABz								! vector between nodes
REAL(dbl) :: dx,dy,dz									! unit vector pointing from A to B
REAL(dbl) :: AB,AP										! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: term1,term2,term3						! terms used in the solution of AP
REAL(dbl) :: t1,t2										! two solutions for AP (AP is the smaller one)
REAL(dbl) :: dotAV,AMag,cosThetaAV,sinThetaAV	! dot product of A and V, magnitudes of A and V, cosine and sine of the angle between A and V
REAL(dbl) :: dotBV,BMag,cosThetaBV					! dot product of B and V, magnitudes of B, angle between B and V
REAL(dbl) :: VMag											! magnitudes of V
INTEGER(lng) :: cylORsphere							! flag to tell whether the solid node is in the cylinder or the sphere (1 for cyl, 2 for sphere)
INTEGER(lng) :: scalarORfluid							! flag to tell whether qCalcV was called for scalar or fluid BC
INTEGER(lng) :: mpierr									! MPI standard error variable

! find the vector between the villous base and the current node (shift coordinate system)
Ax = (x(i) - villiLoc(vNum,1))						! x-coordinate
Ay = (y(j) - villiLoc(vNum,2))						! y-coordinate
Az = (z(k) - villiLoc(vNum,3))						! z-coordinate

! find the vector between the villous base and the solid node (shift coordinate system)
Bx = (x(im1) - villiLoc(vNum,1))						! x-coordinate
By = (y(jm1) - villiLoc(vNum,2))						! y-coordinate
Bz = (z(km1) - villiLoc(vNum,3))						! z-coordinate

! find the vector between the current node and the solid node
ABx = x(im1) - x(i)										! x-coordinate 
ABy = y(jm1) - y(j)										! y-coordinate 				
ABz = z(km1) - z(k)										! z-coordinate 
AB = SQRT(ABx*ABx + ABy*ABy + ABz*ABz)				! magnitude of AB

! find the vector between the villous base and the villous tip
Vx = (villiLoc(vNum,6)-villiLoc(vNum,1))			! x-coordinate
Vy = (villiLoc(vNum,7)-villiLoc(vNum,2))			! y-coordinate
Vz = (villiLoc(vNum,8)-villiLoc(vNum,3))			! z-coordinate

! find the (unit) vector between the current node and solid node
dx = ABx/AB													! x-coordinate
dy = ABy/AB													! y-coordinate
dz = ABz/AB													! z-coordinate

! determine if the point is above or below the top of the villous cylinder (touching a solid node within the hemi-sphere or the within the cylinder)
!		and calculate the distance from the current node to the wall accordingly
dotAV	= Ax*Vx + Ay*Vy + Az*Vz							! dot product of A and V
dotBV	= Bx*Vx + By*Vy + Bz*Vz							! dot product of A and V
AMag	= SQRT(Ax*Ax + Ay*Ay + Az*Az)					! magnitude of A
BMag	= SQRT(Bx*Bx + By*By + Bz*Bz)					! magnitude of B
VMag	= SQRT(Vx*Vx + Vy*Vy + Vz*Vz)					! magnitude of V
cosThetaAV = dotAV/(AMag*VMag)						! cosine of angle between A and V
sinThetaAV = SQRT(1-cosThetaAV*cosThetaAV)		! sine of the angle between A and V
cosThetaBV = dotBV/(BMag*VMag)						! cosine of angle between B and V
IF((AMag*cosThetaAV .GE. VMag) .AND. ((BMag*cosThetaBV .GE. VMag) .OR. (AMag*sinThetaAV .LE. Rv))) THEN	! the solid node is in the sphere
  ! Calculate the distance from the current node (point A) to the wall (point P) (mathematica solution)  						
  term1 = (Ax - Vx)*dx + (Ay - Vy)*dy + (Az - Vz)*dz
  term2 = term1*term1
  term3 = (Ax - Vx)*(Ax - Vx) + (Ay - Vy)*(Ay - Vy) + (Az - Vz)*(Az - Vz) - Rv*Rv
  t1 = -term1 + SQRT(term2 - term3)
  t2 = -term1 - SQRT(term2 - term3)
  AP = MIN(t1,t2)											! distance from the current node (point A) to the villus (point P)							
  q = AP/AB													! distance ratio  
  cylORsphere = 2_lng
ELSE															! solid node is in the cylinder
  ! Calculate the distance from the current node (point A) to the wall (point P) (mathematica solution)
  term1 = Ay*dy*Vx*Vx + Az*dz*Vx*Vx - Ay*dx*Vx*Vy - Ax*dy*Vx*Vy + Ax*dx*Vy*Vy + Az*dz*Vy*Vy - Az*dx*Vx*Vz - Ax*dz*Vx*Vz	&
        - Az*dy*Vy*Vz - Ay*dz*Vy*Vz + Ax*dx*Vz*Vz + Ay*dy*Vz*Vz
  term2 = 0.5_dbl*SQRT((2.0_dbl*Ay*dy*Vx*Vx + 2.0_dbl*Az*dz*Vx*Vx - 2.0_dbl*Ay*dx*Vx*Vy - 2.0_dbl*Ax*dy*Vx*Vy + 2.0_dbl*Ax*dx*Vy*Vy &
        + 2.0_dbl*Az*dz*Vy*Vy - 2.0_dbl*Az*dx*Vx*Vz - 2.0_dbl*Ax*dz*Vx*Vz - 2.0_dbl*Az*dy*Vy*Vz - 2.0_dbl*Ay*dz*Vy*Vz &
        + 2.0_dbl*Ax*dx*Vz*Vz + 2.0_dbl*Ay*dy*Vz*Vz)**2 - 4.0_dbl*(dy*dy*Vx*Vx + dz*dz*Vx*Vx - 2.0_dbl*dx*dy*Vx*Vy + dx*dx*Vy*Vy &
        + dz*dz*Vy*Vy - 2.0_dbl*dx*dz*Vx*Vz - 2.0_dbl*dy*dz*Vy*Vz + dx*dx*Vz*Vz + dy*dy*Vz*Vz)*(Ay*Ay*Vx*Vx + Az*Az*Vx*Vx &
        - 2.0_dbl*Ax*Ay*Vx*Vy + Ax*Ax*Vy*Vy + Az*Az*Vy*Vy - 2.0_dbl*Ax*Az*Vx*Vz - 2.0_dbl*Ay*Az*Vy*Vz + Ax*Ax*Vz*Vz + Ay*Ay*Vz*Vz &
        - Rv*Rv*Vx*Vx - Rv*Rv*Vy*Vy - Rv*Rv*Vz*Vz))
  term3 = dy*dy*Vx*Vx + dz*dz*Vx*Vx - 2.0_dbl*dx*dy*Vx*Vy + dx*dx*Vy*Vy + dz*dz*Vy*Vy - 2.0_dbl*dx*dz*Vx*Vz - 2.0_dbl*dy*dz*Vy*Vz 	&
        + dx*dx*Vz*Vz + dy*dy*Vz*Vz
  t1 = -(term1 + term2)/term3
  t2 = -(term1 - term2)/term3
  AP = MIN(t1,t2)											! distance from the current node (point A) to the villus (point P)							
  q = AP/AB													! distance ratio
  cylORsphere = 1_lng
END IF

! fix slight discrepencies in calculation of q and node tracking
IF((q .LT. 0_lng) .AND. (q .GT. -0.3_dbl)) THEN
  q = 0.00000001_dbl
ELSE IF ((q .GT. 1_lng) .AND. (q .LT. 1.3_dbl)) THEN
  q = 0.99999999_dbl
!ELSE IF(((q .LT. -0.3_dbl) .AND. (q .GT. -2.0_dbl)) .OR. ((q .GT. 1.3_dbl) .AND. (q .LT. 2.0_dbl))) THEN
ELSE
  q = 0.5_dbl
END IF 

!! make sure 0<q<1
!IF((q .LT. 0_lng) .OR. (q .GT. 1_lng)) THEN 
!
!  OPEN(1000,FILE='error.'//sub//'.txt')
!  WRITE(1000,'(A50)') "error in ICBC.f90 at Line 753: q is out of range"
!  WRITE(1000,*) "iter",iter
!  WRITE(1000,*) "vNum",vNum
!  WRITE(1000,*) "m=",m
!  WRITE(1000,*) "i=",i,"j=",j,"k=",k
!  WRITE(1000,*) "x(i)=",x(i),"y(j)=",y(j),"z(k)=",z(k)
!  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
!  WRITE(1000,*) "x(im1)=",x(im1),"y(jm1)=",y(jm1),"z(km1)=",z(km1)
!  WRITE(1000,*) "node(i,j,k)=",node(i,j,k)
!  WRITE(1000,*) "node(im1,jm1,km1)=",node(im1,jm1,km1)
!  WRITE(1000,*) "Ax=",Ax,"Ay=",Ay,"Az=",Az
!  WRITE(1000,*) "Bx=",Bx,"By=",By,"Bz=",Bz
!  WRITE(1000,*) "ABx=",ABx,"ABy=",ABy,"ABz=",ABz
!  WRITE(1000,*) "Vx=",Vx,"Vy=",Vy,"Vz=",Vz
!  WRITE(1000,*) "villiLoc(vNum,1)=",villiLoc(vNum,1),"villiLoc(vNum,2)=",villiLoc(vNum,2),"villiLoc(vNum,3)=",villiLoc(vNum,3)  
!  WRITE(1000,*) "villiLoc(vNum,6)=",villiLoc(vNum,6),"villiLoc(vNum,7)=",villiLoc(vNum,7),"villiLoc(vNum,8)=",villiLoc(vNum,8)
!  WRITE(1000,*) "dx=",dx,"dy=",dy,"dz=",dz
!  WRITE(1000,*) "term1=",term1,"term2=",term2,"term3=",term3
!  WRITE(1000,*) "t1=",t1,"t2=",t2
!  WRITE(1000,*) "AB=",AB,"AP=",AP
!  WRITE(1000,*) "q=",q
!  WRITE(1000,*) "cylORsphere=",cylORsphere
!  WRITE(1000,*) "scalarORfluid=",scalarORfluid
!  CLOSE(1000)
!
!  OPEN(1001,FILE='ABV.'//sub//'.dat')
!  WRITE(1001,'(A60)') 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
!  WRITE(1001,'(3E15.5,6I6)') x(i), y(j), z(k), 0, 0, 0, 0, 0, -vNum
!  WRITE(1001,'(3E15.5,6I6)') x(im1), y(jm1), z(km1), 0, 0, 0, 0, 0, -vNum
!  WRITE(1001,'(3E15.5,6I6)') villiLoc(vNum,1), villiLoc(vNum,2), villiLoc(vNum,3), 0, 0, 0, 0, 0, -vNum
!  WRITE(1001,'(3E15.5,6I6)') villiLoc(vNum,6), villiLoc(vNum,7), villiLoc(vNum,8), 0, 0, 0, 0, 0, -vNum
!  CLOSE(1001)
!
!  CALL PrintFieldsTEST										! output the velocity, density, and scalar fields [MODULE: Output]
!
!  STOP
!
!END IF						

! find the vector from the base of the villus to the point P		
Cx = Ax + AP*dx											! x-coordinate
Cy = Ay + AP*dy											! y-coordinate
Cz = Az + AP*dz											! z-coordinate								

!------------------------------------------------
END SUBROUTINE qCalcV
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE PrintFieldsTEST	! print velocity, density, and scalar to output files
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng)	:: i,j,k,ii,jj,kk,n		! index variables (local and global)
CHARACTER(7)	:: iter_char				! iteration stored as a character

  ! scale the iteration by 1/10 such that the numbers used in the output file aren't too large
  WRITE(iter_char(1:7),'(I7.7)') iter

  ! store the current iteration in "filenum"
  filenum(fileCount) = iter
  fileCount = fileCount + 1_lng

  ! open the proper output file
  OPEN(60,FILE='out-TEST-'//iter_char//'-'//sub//'.dat')
  WRITE(60,*) 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
  WRITE(60,'(A10,E15.5,A5,I4,A5,I4,A5,I4,A8)') 'ZONE T="',iter/(nt/nPers),'" I=',nxSub,' J=',nySub,' K=',nzSub,'F=POINT'

  DO k=1,nzSub
    DO j=1,nySub
      DO i=1,nxSub

         ! convert local i,j,k, to global ii,jj,kk
         ii = ((iMin - 1_lng) + i)
         jj = ((jMin - 1_lng) + j)
         kk = ((kMin - 1_lng) + k)

         WRITE(60,'(8E15.5,I6)') xx(ii), yy(jj), zz(kk), u(i,j,k)*vcf, v(i,j,k)*vcf, w(i,j,k)*vcf, (rho(i,j,k)-denL)*dcf*pcf,	&
                                 phi(i,j,k), node(i,j,k)

      END DO
    END DO
  END DO

  CLOSE(60)

!  ! print villi locations
!  IF(myid .EQ. master) THEN
!
!    OPEN(607,FILE='villi-'//iter_char//'.dat')
!    OPEN(608,FILE='villi2-'//iter_char//'.dat')
!    WRITE(607,'(A60)') 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
!    WRITE(608,'(A60)') 'VARIABLES = "x" "y" "z" "u" "v" "w" "P" "phi" "node"'
!
!    DO n=1,numVilli
!
!      WRITE(607,'(3E15.5,6I4)') villiLoc(n,1), villiLoc(n,2), villiLoc(n,3), 0, 0, 0, 0, 0, 0
!      WRITE(608,'(3E15.5,6I4)') villiLoc(n,6), villiLoc(n,7), villiLoc(n,8), 0, 0, 0, 0, 0, 0
!
!    END DO
!
!    CLOSE(607)
!    CLOSE(608)
!
!  END IF

!------------------------------------------------
END SUBROUTINE PrintFieldsTEST
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE SymmetryBC							! implements symmetry boundary conditions
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) 	:: i,j,k,m,ii,jj,iComm		! index variables
INTEGER(lng)	:: mpierr						! MPI standard error variable

! Loop through the subdomain faces
! -YZ Face
iComm = 2
IF(SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   
  i = YZ_RecvIndex(iComm)	   				! i location of phantom nodes
  ii = YZ_SendIndex(iComm) + 1_lng			! i location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng
    DO j=0,nySub+1_lng
      DO m=1,NumFs_Face

        f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),ii,j,k)			! symmetry BC for 'f'

      END DO

      rho(i,j,k) = rho(ii,j,k)    			! symmetry BC for density 
      phi(i,j,k) = phi(ii,j,k)    			! symmetry BC for scalar

    END DO
  END DO

END IF

! -ZX Face
iComm = 4
IF(SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   
  j = ZX_RecvIndex(iComm)						! j location of phantom nodes
  jj = ZX_SendIndex(iComm) + 1_lng			! j location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng
    DO i=0,nxSub+1_lng
      DO m=1,NumFs_Face

        f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),i,jj,k)        ! symmetry BC

      END DO
 
      rho(i,j,k) = rho(i,jj,k)    			! symmetry BC for density 
      phi(i,j,k) = phi(i,jj,k)    			! symmetry BC for scalar

    END DO
  END DO

END IF

! Z Axis
iComm = 8
IF((SubID(2) .EQ. 0) .AND. (SubID(4) .EQ. 0) .AND. (SubID(iComm) .EQ. 0)) THEN									! if no neighbor in the iComm communication direction exists, iplement symmetry BC

  i = Z_RecvIndex(iComm,1)						! i location of phantom nodes	
  ii = Z_SendIndex(iComm,1) + 1_lng			! i location of 1 row in from the boundary nodes

  j = Z_RecvIndex(iComm,2)						! j location of phantom nodes	
  jj = Z_SendIndex(iComm,2) + 1_lng			! j location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng

    DO m=1,NumFs_Side

        f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),ii,jj,k)        ! symmetry BC

    END DO

    rho(i,j,k) = rho(ii,jj,k)  		  		! symmetry BC for density 
    phi(i,j,k) = phi(ii,jj,k)    			! symmetry BC for scalar

  END DO

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)	! synchronize all processing units

!------------------------------------------------
END SUBROUTINE SymmetryBC
!------------------------------------------------










!--------------------------------------------------------------------------------------------------
SUBROUTINE SymmetryBC_NODE						! implements symmetry boundary conditions
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) 	:: i,j,k,m,ii,jj,iComm		! index variables
INTEGER(lng)	:: mpierr						! MPI standard error variable

! Loop through the subdomain faces
! -YZ Face
iComm = 2
IF(SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   
  i = YZ_RecvIndex(iComm)	   				! i location of phantom nodes
  ii = YZ_SendIndex(iComm) + 1_lng			! i location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng
    DO j=0,nySub+1_lng

      node(i,j,k) = node(ii,j,k)    		! symmetry BC for node flag 

    END DO
  END DO

END IF

! -ZX Face
iComm = 4
IF(SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   
  j = ZX_RecvIndex(iComm)						! j location of phantom nodes
  jj = ZX_SendIndex(iComm) + 1_lng			! j location of 1 row in from the boundary nodes

  DO k=0,nzSub+1_lng
    DO i=0,nxSub+1_lng
 
      node(i,j,k) = node(i,jj,k)    		! symmetry BC for node flag 

    END DO
  END DO

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)	! synchronize all processing units

!------------------------------------------------
END SUBROUTINE SymmetryBC_NODE
!------------------------------------------------







!================================================
END MODULE BClbm
!================================================
