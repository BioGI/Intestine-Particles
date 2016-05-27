!==================================================================================================
MODULE BClbm							     ! Sets LBM Boundary Conditions
!==================================================================================================
USE SetPrecision
USE Setup  
USE MPI
USE IC

CONTAINS




!==================================================================================================
SUBROUTINE BounceBackL(m,i,j,k,im1,jm1,km1,fbb)		
!==================================================================================================
! implements the (moving) bounceback boundary conditions (1st order accurate - Ladd)

IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT)   :: fbb					! bounced back distribution function
REAL(dbl) :: cosTheta, sinTheta					! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb						! wall velocity (x-, y-, z- components)
REAL(dbl) :: rijk						! radius of current node

rijk = SQRT(x(im1)*x(im1) + y(jm1)*y(jm1))			! radius at current location

cosTheta = x(im1)/rijk						! COS(theta)
sinTheta = y(jm1)/rijk						! SIN(theta)

ub = vel(km1)*cosTheta						! x-component of the velocity at i,j,k
vb = vel(km1)*sinTheta						! y-component of the velocity at i,j,k
wb = 0.0_dbl							! no z-component in this case			

fbb = fplus(bb(m),i,j,k) + 6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m))

!==================================================================================================
END SUBROUTINE BounceBackL
!==================================================================================================










!==================================================================================================
SUBROUTINE BounceBack2(m,i,j,k,im1,jm1,km1,fbb)	      !(moving) bounceback BC (2nd order-Lallemand)
!==================================================================================================
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1			! index variables
REAL(dbl), INTENT(OUT) :: fbb					! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2				! index variables
REAL(dbl) :: cosTheta, sinTheta					! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb						! wall velocity (x-, y-, z- components)
REAL(dbl) :: q							! local wall distance ratio [(distance from current node to wall)/(distance to next node in that direction)]
REAL(dbl) :: rijk						! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt		! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 

ip1= i + ex(m)							! i location of 1st neighbor in the m direction
jp1= j + ey(m)							! j location of 1st neighbor in the m direction
kp1= k + ez(m)							! k location of 1st neighbor in the m direction
ip2= i + 2_lng*ex(m)						! i location of 2nd neighbor in the m direction
jp2= j + 2_lng*ey(m)						! j location of 2nd neighbor in the m direction
kp2= k + 2_lng*ez(m)						! k location of 2nd neighbor in the m direction

!----- 2nd order BB if two positive neighbors are in fluid (most cases)
IF ((node(ip1,jp1,kp1) .EQ. FLUID) .AND. (node(ip2,jp2,kp2) .EQ. FLUID)) THEN		
   CALL qCalc_iter(m,i,j,k,im1,jm1,km1,xt,yt,zt,rt,q)

   cosTheta= xt/rt
   sinTheta= yt/rt
   IF (k.NE.km1) THEN
      vt = ((zt-z(k))*vel(km1)+(z(km1)-zt)*vel(k))/(z(km1)-z(k))
   ELSE
      vt = (vel(k)+vel(km1))*0.5_dbl
   ENDIF
   ub = vt* cosTheta						! x-component of the velocity at i,j,k
   vb = vt* sinTheta						! y-component of the velocity at i,j,k
   wb = 0.0_dbl							! no z-component in this case)

   !------ bounced back distribution function with added momentum
   IF ((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN	! use rho = 1.0
      fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 							&
          +   (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1) 						& 
          - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 						&
          + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)) 
      fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
      fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
   ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN ! Use rho = 1.0 
      fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 						&
          + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q							&
          - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)				&
          + (6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))		
      fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
      fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))/(q*(2.0_dbl*q + 1.0_dbl))
   ELSE
      OPEN(1000,FILE='error-'//sub//'.txt')
      WRITE(1000,*) "Error in BounceBack2: q is not between 0 and 1. Aborting."
      WRITE(1000,*) "q=",q,"(i,j,k):",i,j,k
      CLOSE(1000)
      STOP
   END IF
ELSE
   CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)
END IF

!==================================================================================================
END SUBROUTINE BounceBack2
!==================================================================================================








!==================================================================================================
SUBROUTINE qCalc_iter(m,i,j,k,im1,jm1,km1,xt,yt,zt,rt,q)	! calculates q itteratively
!==================================================================================================
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1                 ! index variables
REAL(dbl), INTENT(OUT)   :: q 	    				! ilocal wall distance ratio
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2                         ! index variables
INTEGER(lng) :: it                                              ! loop index variables
REAL(dbl)    :: rijk                                            ! radius of current node
REAL(dbl)    :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt             ! temporary coordinates to search for exact boundary coordinate 

!----- Initial fluid node guess
x1= x(i)
y1= y(j)
z1= z(k)
                
!----- Initial solid node guess
x2= x(im1)
y2= y(jm1)
z2= z(km1)
                 
IF (k.NE.km1) THEN
    DO it=1,15
       !----- guess of boundary location 
       xt= (x1+x2)/2.0_dbl
       yt= (y1+y2)/2.0_dbl
       zt= (z1+z2)/2.0_dbl
       rt= SQRT(xt*xt + yt*yt)
       ht= ((zt-z(k))*r(km1)+(z(km1)-zt)*r(k))/(z(km1)-z(k))
       IF (rt .GT. ht) then
          x2= xt
          y2= yt
          z2= zt
       ELSE
          x1= xt
          y1= yt
          z1= zt
       END IF
    END DO
    x1= x(i)
    y1= y(j)
    z1= z(k)
    x2= x(im1)
    y2= y(jm1)
    z2= z(km1)
    q= sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
 ELSE
    DO it=1,15
       !----- guess of boundary location 
       xt= (x1+x2)/2.0_dbl
       yt= (y1+y2)/2.0_dbl
       zt= (z1+z2)/2.0_dbl
       rt = SQRT(xt*xt + yt*yt)
       ht = (r(km1)+r(k))/2.0_dbl
       IF (rt.GT.ht) then
          x2= xt
          y2= yt
          z2= zt
       ELSE
          x1= xt
          y1= yt
          z1= zt
       END IF
    END DO
    x1= x(i)
    y1= y(j)
    z1= z(k)
    x2= x(im1)
    y2= y(jm1)
    z2= z(km1)
    q= sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
END IF

!==================================================================================================
END SUBROUTINE qCalc_iter 
!==================================================================================================








!==================================================================================================
SUBROUTINE qCalc(m,i,j,k,im1,jm1,km1,q)			! calculates q (boundary distance ratio) using "ray tracing" - see wikipedia article
!==================================================================================================
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1		! current node, and neighboring node
REAL(dbl), INTENT(OUT) :: q				! distance ratio
REAL(dbl) :: Ax,Ay,Az					! current node
REAL(dbl) :: Bx,By,Bz					! solid node
REAL(dbl) :: AB,AP					! distances between current and solid nodes, and between current node and the wall
REAL(dbl) :: dx,dy,dz					! unit vector pointing from A to B
REAL(dbl) :: r1,r2,z1,z2,slope,intercept		! radius and z location at k and km1, slope of line connecting those two points, z-intercept of the r-equation
REAL(dbl) :: slope2,term1,term2				! terms used in calculation

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
dx = (x(im1)-x(i))/AB						! i direction
dy = (y(jm1)-y(j))/AB						! j direction
dz = (z(km1)-z(k))/AB						! k direction

! SURFACE
r1 = r(k)							! radius at k (distance from CL)
r2 = r(km1)							! radius at km1 (distance from CL)
z1 = z(k)							! z-coordinate at k
z2 = z(km1)							! z-coordinate at km1

IF(k .NE. km1) THEN
  slope = (r2-r1)/(z2-z1)					! approximate the surface as a conincal shell (linear between k values)
ELSE
  slope = 0.0_dbl
END IF

intercept = r1 - slope*z1					! z-intercept of the linearly approximated r-equation

!----- terms used in calculation
slope2 = slope*slope						! slope^2
term1 = Ax*dx + Ay*dy - Az*dz*slope2 - intercept*dz*slope	! reoccuring term
term2 = dx*dx + dy*dy - dz*dz*slope*slope			! reoccuring term

!------ calculate the distance from the current node (point A) to the wall (point P)
AP = (1.0_dbl/(2.0_dbl*term2)) * &
     (-2.0_dbl*term1					&
   + SQRT(4.0_dbl*(term1*term1 - (Ax*Ax + Ay*Ay - intercept*intercept - 2.0_dbl*Az*intercept*slope - Az*Az*slope2)*term2)))

q = AP/AB							! distance ratio
q = max(q,0.001_dbl)

!----- balaji added
!q=0.5

! make sure 0<q<1
IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
  OPEN(1000,FILE="error.txt")
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

!==================================================================================================
END SUBROUTINE qCalc
!==================================================================================================









!==================================================================================================
SUBROUTINE PrintFieldsTEST		       ! printis velocity, density & scalar to output files
!==================================================================================================
IMPLICIT NONE

INTEGER(lng)	:: i,j,k,ii,jj,kk,n		! index variables (local and global)
CHARACTER(7)	:: iter_char			! iteration stored as a character

  !----- scale the iteration by 1/10 such that the numbers used in the output file aren't too large
  WRITE(iter_char(1:7),'(I7.7)') iter

  !----- store the current iteration in "filenum"
  filenum(fileCount) = iter
  fileCount = fileCount + 1_lng

  !----- open the proper output file
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

!==================================================================================================
END SUBROUTINE PrintFieldsTEST
!==================================================================================================








!==================================================================================================
SUBROUTINE SymmetryBC			   				   ! implements symmetry BC
!==================================================================================================
IMPLICIT NONE

INTEGER(lng) 	:: i,j,k,m,ii,jj,iComm				! index variables
INTEGER(lng)	:: mpierr					! MPI standard error variable

!----- Loop through the subdomain faces
!----- YZ Face
iComm = 2
IF (SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   i = YZ_RecvIndex(iComm)					! i location of phantom nodes
   ii = YZ_SendIndex(iComm) + 1_lng				! i location of 1 row in from the boundary nodes
   DO k=0,nzSub+1_lng
      DO j=0,nySub+1_lng
         DO m=1,NumFs_Face
            f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),ii,j,k)	! symmetry BC for 'f'
         END DO
         rho(i,j,k) = rho(ii,j,k)    				! symmetry BC for density 
         phi(i,j,k) = phi(ii,j,k)    				! symmetry BC for scalar
      END DO
   END DO
END IF

!----- ZX Face
iComm = 4
IF (SubID(iComm) .EQ. 0) THEN					! if no neighbor in the iComm communication direction exists, implement symmetry BC
   j = ZX_RecvIndex(iComm)					! j location of phantom nodes
   jj = ZX_SendIndex(iComm) + 1_lng				! j location of 1 row in from the boundary nodes
   DO k=0,nzSub+1_lng
      DO i=0,nxSub+1_lng
         DO m=1,NumFs_Face
            f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),i,jj,k)        ! symmetry BC
         END DO
         rho(i,j,k) = rho(i,jj,k)    				! symmetry BC for density 
         phi(i,j,k) = phi(i,jj,k)    				! symmetry BC for scalar
      END DO
   END DO
END IF

!----- Z Axis
iComm = 8
IF ((SubID(2) .EQ. 0) .AND. (SubID(4) .EQ. 0) .AND. (SubID(iComm) .EQ. 0)) THEN		! if no neighbor in the iComm communication direction exists, iplement symmetry BC
   i = Z_RecvIndex(iComm,1)					! i location of phantom nodes	
   ii = Z_SendIndex(iComm,1) + 1_lng				! i location of 1 row in from the boundary nodes
   j = Z_RecvIndex(iComm,2)					! j location of phantom nodes	
   jj = Z_SendIndex(iComm,2) + 1_lng				! j location of 1 row in from the boundary nodes

   DO k=0,nzSub+1_lng
      DO m=1,NumFs_Side
         f(f_Comps(OppCommDir(iComm),m),i,j,k) = f(sym(f_Comps(OppCommDir(iComm),m),iComm),ii,jj,k)        ! symmetry BC
      END DO
      rho(i,j,k) = rho(ii,jj,k)  	  			! symmetry BC for density 
      phi(i,j,k) = phi(ii,jj,k)    				! symmetry BC for scalar
   END DO
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)				! synchronize all processing units

!==================================================================================================
END SUBROUTINE SymmetryBC
!==================================================================================================








!==================================================================================================
SUBROUTINE SymmetryBC_NODE				! implements symmetry boundary conditions
!==================================================================================================
IMPLICIT NONE

INTEGER(lng) 	:: i,j,k,m,ii,jj,iComm			! index variables
INTEGER(lng)	:: mpierr				! MPI standard error variable

!----- Loop through the subdomain faces
!------YZ Face
iComm = 2
IF (SubID(iComm) .EQ. 0) THEN				! if no neighbor in the iComm communication direction exists, implement symmetry BC
   i = YZ_RecvIndex(iComm)	   			! i location of phantom nodes
   ii = YZ_SendIndex(iComm) + 1_lng			! i location of 1 row in from the boundary nodes
   DO k=0,nzSub+1_lng
      DO j=0,nySub+1_lng
         node(i,j,k) = node(ii,j,k)    			! symmetry BC for node flag 
      END DO
   END DO
END IF

!----- ZX Face
iComm = 4
IF (SubID(iComm) .EQ. 0) THEN				! if no neighbor in the iComm communication direction exists, implement symmetry BC
   j = ZX_RecvIndex(iComm)				! j location of phantom nodes
   jj = ZX_SendIndex(iComm) + 1_lng			! j location of 1 row in from the boundary nodes
    DO k=0,nzSub+1_lng
       DO i=0,nxSub+1_lng
           node(i,j,k) = node(i,jj,k)    		! symmetry BC for node flag 
       END DO
    END DO
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)	! synchronize all processing units
!==================================================================================================
END SUBROUTINE SymmetryBC_NODE
!==================================================================================================






!==================================================================================================
END MODULE BClbm
!==================================================================================================
