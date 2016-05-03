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

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! index variables
REAL(dbl), INTENT(OUT) :: fbb			! bounced back distribution function
INTEGER(lng) :: ip1,jp1,kp1,ip2,jp2,kp2		! index variables
REAL(dbl) :: cosTheta, sinTheta			! COS(theta), SIN(theta)
REAL(dbl) :: ub, vb, wb				! wall velocity (x-, y-, z- components)
REAL(dbl) :: q					! local wall distance ratio [(distance to wall)/(distance to next node in that direction)]
REAL(dbl) :: rijk				! radius of current node
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt,vt! temporary coordinates to search for exact boundary coordinate (instead of ray tracing) 
INTEGER(lng) :: it				! loop index variables

ip1 = i + ex(m)					! i location of 1st neighbor in the m direction
jp1 = j + ey(m)					! j location of 1st neighbor in the m direction
kp1 = k + ez(m)					! k location of 1st neighbor in the m direction

ip2 = i + 2_lng*ex(m)				! i location of 2nd neighbor in the m direction
jp2 = j + 2_lng*ey(m)				! j location of 2nd neighbor in the m direction
kp2 = k + 2_lng*ez(m)				! k location of 2nd neighbor in the m direction


IF ((node(ip1,jp1,kp1) .EQ. FLUID) .AND. (node(ip2,jp2,kp2) .EQ. FLUID)) THEN 	!2nd order BB if two positive neighbors are in fluid (most cases)
   CALL qCalc(m,i,j,k,im1,jm1,km1,q) 
!  CALL qCalcFarhad(i,q)

   ub = velIn(km1)								! x-component of the velocity at i,j,k
   vb = 0.0_dbl									! y-component of the velocity at i,j,k
   wb = 0.0_dbl									! only z-component in this case)

   !------ make sure 0<q<1
   IF ((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
      OPEN(1000,FILE="error.txt")
      WRITE(1000,*) "q=",q
      WRITE(1000,*) "m=",m
      WRITE(1000,*) "i=",i,"j=",j,"k=",k
      WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
      CLOSE(1000)
      STOP
   END IF	


   !----- bounced back distribution function with added momentum
   IF ((q .LT. 0.5_dbl) .AND. (q .GT. -0.00000001_dbl)) THEN
      fbb = q*(1.0_dbl + 2.0_dbl*q)*fplus(bb(m),i,j,k) 				&
          + (1.0_dbl - 4.0_dbl*q*q)*fplus(bb(m),ip1,jp1,kp1)    	  	& 
          - q*(1.0_dbl - 2.0_dbl*q)*fplus(bb(m),ip2,jp2,kp2) 			&
          + 6.0_dbl*wt(m)*1.0_dbl*(ub*ex(m) + vb*ey(m) + wb*ez(m)) 		! use rho = 1.0

	fmovingsum = fmovingsum + (6.0_dbl*wt(m)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
	fmovingrhosum = fmovingrhosum + (6.0_dbl*wt(m)*rho(i,j,k)*(ub*ex(m) + vb*ey(m) + wb*ez(m)))
 
   ELSE IF((q .GE. 0.5_dbl) .AND. (q .LT. 1.00000001_dbl)) THEN
      fbb = fplus(bb(m),i,j,k)/(q*(2.0_dbl*q + 1.0_dbl)) 			&
        + ((2.0_dbl*q - 1.0_dbl)*fplus(m,i,j,k))/q	    		   	&
        - ((2.0_dbl*q - 1.0_dbl)/(2.0_dbl*q + 1.0_dbl))*fplus(m,ip1,jp1,kp1)	&
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
END SUBROUTINE BounceBack2
!------------------------------------------------





!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalc(m,i,j,k,im1,jm1,km1,q)			! calculates q (boundary distance ratio) using "ray tracing" - see wikipedia article
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1	! current node, and neighboring node
REAL(dbl), INTENT(OUT) :: q							! distance ratio
REAL(dbl) :: rijk ,htt,xo
REAL(dbl) :: x1,y1,z1,x2,y2,z2,xt,yt,zt,ht,rt
INTEGER(lng) :: it				! loop index variables


   rijk = x(im1)                    		! height at current location
   !----- Initial fluid node guess
   x1=x(i)
   y1=y(j)
   z1=z(k)
                
   !----- Initial solid node guess
   x2=x(im1)
   y2=y(jm1)
   z2=z(km1)
   
   xo= (rOut(k) +rIn(k)) / 2.0_dbl               
   
   IF (k.NE.km1) THEN
      DO it=1,qitermax
         !----- guess of boundary location 
         xt=(x1+x2)/2.0_dbl
         yt=(y1+y2)/2.0_dbl
         zt=(z1+z2)/2.0_dbl


         IF (rijk .GE. rOut(k)) THEN
            ht = ((zt-z(k))*rOut(km1)+(z(km1)-zt)*rOut(k))/(z(km1)-z(k))
         ELSE
            ht = ((zt-z(k))*rIn(km1)+(z(km1)-zt)*rIn(k))/(z(km1)-z(k))
         END IF

         htt= abs(ht-xo)
         rt = abs(xt-xo)

         IF (rt.GT.htt) then
            x2=xt
            y2=yt
            z2=zt
         ELSE
            x1=xt
            y1=yt
            z1=zt
         END IF
         IF ((I.EQ.97) .AND. (j.EQ.1) .AND. (k.EQ.1))THEN
            !write(*,*) 'A:',x1,x2,xt,rt,ht
         END IF

      END DO

      x1=x(i)
      y1=y(j)
      z1=z(k)
                 
      x2=x(im1)
      y2=y(jm1)
      z2=z(km1)
 
      q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

      IF ((I.EQ.97) .AND. (j.EQ.1) .AND. (k.EQ.1))THEN
         !write(*,*) 'B:', x1,y1,z1,x2,y2,z2,xt,yt,zt,rt,ht,q
      END IF


   ELSE
      DO it=1,qitermax
         !----- guess of boundary location 
         xt=(x1+x2)/2.0_dbl
         yt=(y1+y2)/2.0_dbl
         zt=(z1+z2)/2.0_dbl


         IF (rijk .GE. rOut(k)) THEN
            ht = (rOut(km1)+rOut(k))/2.0_dbl
         ELSE
            ht = (rIn(km1)+rIn(k))/2.0_dbl
         END IF

         htt= abs(ht-xo)
         rt = abs(xt-xo) 

         IF (rt.GT.htt) then
            x2=xt
            y2=yt
            z2=zt
         ELSE
            x1=xt
            y1=yt
            z1=zt
         END IF

         IF ((I.EQ.18) .AND. (j.EQ.1) .AND. (k.EQ.1))THEN
            !write(*,*) 'A:',x1,x2,xt,rt,ht
         END IF
      END DO

      x1=x(i)
      y1=y(j)
      z1=z(k)
                 
      x2=x(im1)
      y2=y(jm1)
      z2=z(km1)
 
      q=sqrt((xt-x1)**2+(yt-y1)**2+(zt-z1)**2)/sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
      IF ((I.EQ.18) .AND. (j.EQ.1) .AND. (k.EQ.1))THEN
         !write(*,*) 'B:', x1,y1,z1,x2,y2,z2,xt,yt,zt,rt,ht,q
      END IF

   ENDIF


! make sure 0<q<1
IF((q .LT. -0.00000001_dbl) .OR. (q .GT. 1.00000001_dbl)) THEN 
  OPEN(1000,FILE="error.txt")
  WRITE(1000,*) "q    =",q
  WRITE(1000,*) "m=",m
  WRITE(1000,*) "i=",i,"j=",j,"k=",k
  WRITE(1000,*) "im1=",im1,"jm1=",jm1,"km1=",km1
  CLOSE(1000)
  STOP
END IF																																

!------------------------------------------------
END SUBROUTINE qCalc
!------------------------------------------------






!--------------------------------------------------------------------------------------------------
SUBROUTINE qCalcFarhad(i,q)  
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: i
REAL(dbl)    :: h1,h2,time,D_X,D_Y,q

time = iter*tcf
D_X= 20.0_dbl *D 
D_Y= 0.50_dbl *D

h2= 0.5_dbl*D_x - 0.38_dbl*D_x + 5.0000e-5 + (s1*time)    
h1= 0.5_dbl*D_x - 0.48_dbl*D_x + 5.0000e-5 + (s1*time) 
IF (x(i) .LT. 0.5*(h1+h2) ) then 		  			!Left    
   q= (x(i)-h1)/xcf
ELSE							  	 	!Right 
   q= (h2-x(i))/xcf
   !write(*,*) i,x(i),h2,q 
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
