!==================================================================================================
MODULE LBM				! LBM Subroutines (Equilibrium, Collision, Stream, Macro, Scalar, ScalarBCs,ParticleTracking)
!==================================================================================================
USE SetPrecision
USE Setup
USE IC
USE BClbm
USE MPI

IMPLICIT NONE

CONTAINS

!===================================================================================================
SUBROUTINE LBM_Setup	! setup the LBM simulation
!===================================================================================================
IMPLICIT NONE

!--Initialize variables and arrays------------------------------------------------------------------
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
!----iComm=2, -ZY FACE
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

!----iComm=3, -ZX FACE
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

!----iComm=8, Z AXIS
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

!----Fill velocity direction vector arrays --------------------------------------------------------
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
nuL   	   = (2.0_dbl*tau - 1.0_dbl)/6.0_dbl	! lattice kinematic viscosity
denL 	   = 1.0_dbl				! arbitrary lattice density (1.0 for convenience)
oneOVERtau = 1.0_dbl/tau			! reciprical of tau
cs	   = (1.0_dbl)/(SQRT(3.0_dbl))		! speed of sound on the lattice

!---- Initialize timestep --------------------------------------------------------------------------
iter = 0_lng					! intialize the starting timestep to 0 - will get reset in 'ICs' in ICBCM.f90

!---- Calculate feq for initial condition ----------------------------------------------------------
CALL Equilibrium

!---- Set f-from wall motion sums to zero at initial timestep --------------------------------------
fmovingsum = 0.0_dbl
fmovingrhosum = 0.0_dbl
!===================================================================================================
END SUBROUTINE LBM_Setup
!===================================================================================================








!===================================================================================================
SUBROUTINE Equilibrium	     !calculates the equilibrium distribution function and set f to feq (IC)
!===================================================================================================
IMPLICIT NONE

INTEGER(lng):: i,j,k,m				! index variables
REAL(dbl)   :: uu,ue,ve,we,Usum			! precalculated quantities for use in feq equation
REAL(dbl)   :: feq				! equilibrium distribution function

DO k=1,nzSub+0
  DO j=1,nySub+0
    DO i=1,nxSub+0
      IF(node(i,j,k) .EQ. FLUID) THEN
        uu = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)		! u . u
        DO m=0,NumDistDirs
          ue  = u(i,j,k)*ex(m)							! u . e
          ve  = v(i,j,k)*ey(m)							! v . e
          we  = w(i,j,k)*ez(m)							! w . e
          Usum= ue + ve + we							! U . e
          feq = (wt(m)*rho(i,j,k))*(1.0_dbl + 3.0_dbl*Usum + 4.5_dbl*Usum*Usum- 1.5_dbl*uu)	
          f(m,i,j,k) = feq    
        END DO
      END IF
    END DO
  END DO
END DO
!===================================================================================================
END SUBROUTINE Equilibrium
!===================================================================================================








!===================================================================================================
SUBROUTINE Collision		  ! equilibrium distribution function & collision step for each node
!===================================================================================================
IMPLICIT NONE

INTEGER(lng):: i,j,k,m		  ! index variables
REAL(dbl)   :: UU,ue,ve,we,Usum	  ! precalculated quantities for use in the feq equation
REAL(dbl)   :: feq		  ! equilibrium distribution function

DO k=1,nzSub+0
   DO j=1,nySub+0
      DO i=1,nxSub+0
         IF (node(i,j,k) .EQ. FLUID) THEN
            UU = u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)		! U . U
            DO m=0,NumDistDirs
               ue  = u(i,j,k)*ex(m)							! u . e
               ve  = v(i,j,k)*ey(m)							! v . e
               we  = w(i,j,k)*ez(m)							! w . e
               Usum= ue + ve + we							! U . e
               feq = (wt(m)*rho(i,j,k))*(1.0_dbl + 3.0_dbl*Usum+ 4.5*Usum*Usum- 1.5*uu)	! equilibrium distribution function
               f(m,i,j,k)= f(m,i,j,k) - oneOVERtau*(f(m,i,j,k) - feq)		    	! collision
            END DO 
         END IF
      END DO
   END DO
END DO
!===================================================================================================
END SUBROUTINE Collision
!===================================================================================================








!===================================================================================================
SUBROUTINE Stream	
!===================================================================================================
! stream the distribution functions between neighboring nodes (using Lallemand 2nd order moving BB)
IMPLICIT NONE

INTEGER(lng) :: i,j,k,m,im1,jm1,km1		! index variables
REAL(dbl) :: fbb				! bounced back distribution function

fplus = f					! store the post-collision distribution function
!----- interior nodes (away from other subdomains)
DO k=2,nzSub-1
  DO j=2,nySub-1
    DO i=2,nxSub-1
      IF(node(i,j,k) .EQ. FLUID) THEN
        DO m=1,NumDistDirs
          !----- i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)
          IF(node(im1,jm1,km1) .EQ. FLUID) THEN 
            f(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN						! macro- boundary
!           CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)						! implement the bounceback BCs 
            CALL BounceBack2(m,i,j,k,im1,jm1,km1,fbb)                                           ! implement the bounceback BCs 
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN	! villi
          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM.f90 at Line 274: node(im1,jm1,km1) is out of range"
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

!----- XY faces
DO k=1,nzSub,(nzSub-1)
  DO j=1,nySub
    DO i=1,nxSub
      IF(node(i,j,k) .EQ. FLUID) THEN
        DO m=1,NumDistDirs
          ! i,j,k location of neighboring node
          im1 = i - ex(m)
          jm1 = j - ey(m)
          km1 = k - ez(m)
          IF(node(im1,jm1,km1) .EQ. FLUID) THEN 
            f(m,i,j,k) = fplus(m,im1,jm1,km1)
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN							! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  				! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi

          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM.f90 at Line 312: node(im1,jm1,km1) is out of range"
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

!----- XZ faces
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
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN							! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  				! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi

          ELSE
            OPEN(1000,FILE="error.txt")
            WRITE(1000,'(A75)') "error in LBM.f90 at Line :350 node(im1,jm1,km1) is out of range"
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

!----- YZ faces
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
          ELSE IF(node(im1,jm1,km1) .EQ. SOLID) THEN		! macro- boundary
            CALL BounceBackL(m,i,j,k,im1,jm1,km1,fbb)			  				! implement the bounceback BCs (Ladd BB) [MODULE: ICBC]
            f(m,i,j,k) = fbb
          ELSE	IF((node(im1,jm1,km1) .LE. -1) .AND. (node(im1,jm1,km1) .GE. -numVilli)) THEN		! villi

          ELSE
            OPEN(1000,FILE="error.txt")
	    WRITE(1000,'(A75)') "error in LBM.f90 at Line 388: node(im1,jm1,km1) is out of range"
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
!===================================================================================================
END SUBROUTINE Stream
!===================================================================================================








!===================================================================================================
SUBROUTINE Macro	! calculate the macroscopic quantities
!===================================================================================================
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
          CALL PrintFieldsTEST						! output the velocity, density, and scalar fields [MODULE: Output]
          STOP
        END IF
      ELSE
        rho(i,j,k)= denL						! density (zero gauge pressure)
        u(i,j,k)  = 0.0_dbl						! x-velocity
        v(i,j,k)  = 0.0_dbl						! y-velocity
        w(i,j,k)  = 0.0_dbl						! z-velocity
        phi(i,j,k)= phiWall						! scalar
      END IF
    END DO
  END DO
END DO   
!===================================================================================================
END SUBROUTINE Macro
!===================================================================================================





!===================================================================================================
END MODULE LBM
!===================================================================================================
