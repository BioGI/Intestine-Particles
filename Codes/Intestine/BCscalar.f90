!==================================================================================================
MODULE BCscalar		! Sets Scalar Boundary Conditions
!==================================================================================================
USE SetPrecision
USE Setup  
USE IC
USE BClbm
USE MPI

CONTAINS






!===================================================================================================
SUBROUTINE BC_Scalar(m,i,j,k,im1,jm1,km1,phiBC)				! implements the scalar BCs 
!===================================================================================================
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: m,i,j,k,im1,jm1,km1				! index variables
REAL(dbl),    INTENT(OUT):: phiBC     					! scalar contribution from the boundary condition
INTEGER(lng) :: mm,ip1,jp1,kp1,iB,jB,kB,xaxis   				! First neighboring node location
REAL(dbl)    :: rhoAstar, phiAstar, feq_Astar,  PkAstar,PkA 		! density at boundary and contribution of scalar from boundary
REAL(dbl)    :: rhoBstar, phiBstar, fPlusBstar, PkBstar 		! Values interpolated to Bstar location
REAL(dbl)    :: cosTheta, sinTheta					! COS(theta), SIN(theta)
REAL(dbl)    :: ub, vb, wb						! wall velocity (x-, y-, z- components)
REAL(dbl)    :: rijk 							! radius of the solid node
REAL(dbl)    :: Geom_norm_x,Geom_norm_y,Geom_norm_z
REAL(dbl)    :: q, q1, n_prod, n_prod_max
REAL(dbl)    :: xt,yt,zt,rt,vt						! Location of the boundary between i,j,k node and im1,jm1,km1 node
!===========================================================================
! HELP: How to set different boundary conditions
! BC Scalar-zero:       coeffPhi=1      coeffGrad= 0    coeffConst= 0
! BC Scalar-Non-zero:   coeffPhi=1      coeffGrad= 0    coeffConst= phi_BC
! BC-Flux-zero:         coeffPhi=0      coeffGrad= 1    coeffConst= 0
! BC-Flux-Non-zero:     coeffPhi=0      coeffGrad= 1    coeffConst= dphi/dn
! BC-Permeability:      coeffPhi=Pw/Dm  coeffGrad=-1    coeffConst= 0
!===========================================================================

!CALL qCalc(m,i,j,k,im1,jm1,km1,q)
!cosTheta = x(im1)/r(km1) 
!sinTheta = y(jm1)/r(km1)  
!ub = vel(km1)*cosTheta  		      			! x-component of the velocity at i,j,k
!vb = vel(km1)*sinTheta        					! y-component of the velocity at i,j,k
!wb = 0.0_dbl

CALL qCalc_iter(m,i,j,k,im1,jm1,km1,xt,yt,zt,rt,q)

IF (Flag_Couette) THEN
   ub= 0.0_dbl
   vb= 0.0_dbl
   xaxis= ANINT(0.5_dbl*(nx+1))
   IF ((iMin-1+im1) .GT. xaxis) THEN 
      wb = vel(k)
   ELSE 
      wb= -vel(k)
   END IF
ELSE
   cosTheta= xt/rt
   sinTheta= yt/rt
   IF (k.NE.km1) THEN
      vt = ((zt-z(k))*vel(km1)+(z(km1)-zt)*vel(k))/(z(km1)-z(k))
   ELSE
      vt = (vel(k)+vel(km1))*0.5_dbl
   ENDIF
   ub = vt* cosTheta						! x-component of the velocity at i,j,k
   vb = vt* sinTheta						! y-component of the velocity at i,j,k
   wb = -s_movingF/vcf 							! no z-component in this case)
END IF

!---------------------------------------------------------------------------------------------------
!----- Computing phi at the wall in case of Dirichlet BC -------------------------------------------
!---------------------------------------------------------------------------------------------------
IF (coeffGrad .EQ. 0.0) then   
   phiWall= coeffConst/coeffPhi
END IF      

!---------------------------------------------------------------------------------------------------
!----- Estimating  phi at the wall in case of Neumann or Mixed BC ----------------------------------
!----- Two Fluid nodes are needed for extrapolation ------------------------------------------------(Counter1.LT.4)) THEN

!----- 1st node; A, is adjacent to the boundary ----------------------------------------------------
!----- 2nd node; B, is a fluid neighbor of A, in the direction closest to geometry-normal ----------
!---------------------------------------------------------------------------------------------------
IF (coeffGrad .NE. 0.0) then
   Geom_norm_x= 1.0 
   Geom_norm_y= 0.0
   Geom_norm_z= 0.0
   n_prod_max= 0.0_dbl
   !----- Finding the mth direction closest to normal vector
   !----- Only one fluid neighboring node is needed for interpolation
   DO mm=0,NumDistDirs
      ip1= i+ ex(mm)
      jp1= j+ ey(mm)
      kp1= k+ ez(mm)
      IF (node(ip1,jp1,kp1) .EQ. FLUID) THEN
         n_prod= abs( Geom_norm_x * (ip1-i) + Geom_norm_y * (jp1-j) + Geom_norm_z * (kp1-k))
         IF (n_prod_max .LT. n_prod) THEN
            n_prod_max= n_prod
            iB= ip1
            jB= jp1
            kB= kp1
         END IF
      END IF
   END DO
   phiWall= ((phiTemp(i,j,k)*(1.0+q)*(1.0+q)/(1.0+2.0*q)) -    	&
             (phiTemp(iB,jB,kB)*q*q/(1.0+2.0*q)) -          	&
             (q*(1+q)/(1+2.0*q))* (coeffConst/coeffGrad)    )  	&
           /(1.0- (q*(1+q)/(1+2.0*q))*(coeffPhi/coeffGrad))
END IF

!----- neighboring node (fluid side) ---------------------------------------------------------------
ip1 = i + ex(m)
jp1 = j + ey(m)
kp1 = k + ez(m)
!------ This rarely happens (both neighboring nodes over a line are solid)
!------  use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF

!----- Computing values at A* & the scalar streamed from A* (Chpter 3 paper) -----------------------
rhoAstar= (rho(i,j,k)-rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)		! Extrapolate density
CALL Equilibrium_LOCAL(m,rhoAstar,ub,vb,wb,feq_Astar)    		! f_eq in mth direction
phiAstar= phiWall							! phi at solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar			! Contribution from A* to B*  

!----- Computing values at B* & the scalar streamed from B* (Chpter 3 paper) -----------------------
!rhoBstar=   (1-q)*rho(ip1,jp1,kp1)     + q* rho(i,j,k)
!phiBstar=   (1-q)*phiTemp(ip1,jp1,kp1) + q* phiTemp(i,j,k)
!fPlusBstar= (1-q)*fplus(m,ip1,jp1,kp1) + q* fplus(m,i,j,k)
!PkBstar=    (fplusBstar/rhoBstar - wt(m)*Delta)*phiBstar

!----- Scalar contribution from wall to the node----------------------------------------------------
!phiBC=      PkAstar+ (PkAstar- PkBstar)*(1-q)

!----- Using only A and A* for interpolation (instead of A* and B*) 
PkA= (fplus(m,i,j,k)/rho(i,j,k) - wt(m)*Delta)*phiTemp(i,j,k)		! contribution from current node to next in the mth direction
IF(q .LT. 0.25) THEN
  q = 0.25_dbl
END IF
phiBC	= ((PkAstar - PkA)/q) + PkAstar	

!===================================================================================================
END SUBROUTINE BC_Scalar
!===================================================================================================










!===================================================================================================
SUBROUTINE BC_Scalar_Permeability_1storder(m,i,j,k,im1,jm1,km1,phiBC)		 ! implements the scalar BCs 
!===================================================================================================
IMPLICIT NONE

REAL(dbl),   INTENT(OUT):: phiBC           				 ! scalar contribution from the boundary condition
INTEGER(lng),INTENT(IN) :: m,i,j,k,im1,jm1,km1		 ! index variables
INTEGER(lng) :: ix0,ix1,iy0,iy1,iz0,iz00,iz1,iz11	 ! Trilinear interpolation parameters
INTEGER(lng) :: P1_N_Solid_nodes, P2_N_Solid_nodes
INTEGER(lng) :: ip1,jp1,kp1,ii,jj,kk,mm,N_Case            ! First neighboring node location
INTEGER(lng) :: xaxis,yaxis
INTEGER(lng) :: phi_N,key1,key2,key3
INTEGER(lng) :: Counter1,Counter2,Interpolation_Dim,Communication_Dim
REAL(dbl)    :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd,dd! Trilinear interpolation parameters
REAL(dbl)    :: xt,yt,zt,rt,vt,q                   ! Location of the boundary between i,j,k node and im1,jm1,km1 node
REAL(dbl)    :: Geom_nx,Geom_ny,Geom_nz,Geom_n_mag
REAL(dbl)    :: Ax,Ay,Az,A_mag
REAL(dbl)    :: Prod
REAL(dbl)    :: rhoAstar, phiAstar, feq_Astar,  PkAstar,PkA	! density at boundary and contribution of scalar from boundary
REAL(dbl)    :: rhoBstar, phiBstar, fPlusBstar, PkBstar 		! Values interpolated to Bstar location
REAL(dbl)    :: rhoWall,rho1,rho2,rho3,rho_sum
REAL(dbl)    :: cosTheta, sinTheta					 
REAL(dbl)    :: ub, vb, wb
REAL(dbl)    :: P1_x,P1_y,P1_z,P2_x,P2_y,P2_z
REAL(dbl)    :: P1_phi,P2_phi,phiWall_new, Del_phiWall,DphiDn
REAL(dbl)    :: alpha,Diffusivity,Pww
REAL(dbl)    :: phi_sum, phi_ave,phi1,phi2,phi3

CALL qCalc_iter(m,i,j,k,im1,jm1,km1,xt,yt,zt,rt,q)

IF (Flag_Couette) THEN
   ub= 0.0_dbl
   vb= 0.0_dbl
   xaxis= ANINT(0.5_dbl*(nx+1))
   IF ((iMin-1+im1) .GT. xaxis) THEN 
      wb = vel(k)
   ELSE 
      wb= -vel(k)
   END IF
ELSE
   cosTheta= xt/rt
   sinTheta= yt/rt
   IF (k.NE.km1) THEN
      vt = ((zt-z(k))*vel(km1)+(z(km1)-zt)*vel(k))/(z(km1)-z(k))
   ELSE
      vt = (vel(k)+vel(km1))*0.5_dbl
   ENDIF
   ub = vt* cosTheta						
   vb = vt* sinTheta					
   wb = -s_movingF/vcf 
END IF

!--- Computing the normal to the geometry based on the equation for boundary 
Geom_nx= -xt/rt
Geom_ny= -yt/rt
Geom_nz= -amp1 *(2.0_dbl*PI/lambda1) *SIN(PI+(2.0_dbl *PI*zt/lambda1)) 
!---normalizing the geometry normal vector 
Geom_n_mag=sqrt(Geom_nx**2.0_dbl + Geom_ny**2.0_dbl + Geom_nz**2.0_dbl)
Geom_nx=Geom_nx/Geom_n_mag
Geom_ny=Geom_ny/Geom_n_mag
Geom_nz=Geom_nz/Geom_n_mag
!--------------------------------------------------------------------------------------------------
!--- Finding location of the point, P1, which is one mesh size away from (xt,yt,zt) at the boundary
alpha=1.0_dbl               ! The coefficient which defines the distance to walk away from the boundary
xaxis=ANINT(0.5_dbl*(nx+1))
yaxis=ANINT(0.5_dbl*(ny+1))
P1_x= ((xt + alpha*Geom_nx*xcf)/xcf) - iMin + xaxis + 1
P1_y= ((yt + alpha*Geom_ny*xcf)/xcf) - jMin + yaxis + 1
P1_z= ((zt + alpha*Geom_nz*xcf)/xcf) - kMin + 1

ix0= FLOOR(P1_x)
ix1= CEILING(P1_x)
iy0= FLOOR(P1_y)
iy1= CEILING(P1_y)
iz0= FLOOR(P1_z)
iz1= CEILING(P1_z)

P1_N_Solid_nodes =   node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1) 

IF(P1_N_Solid_nodes.EQ.0)THEN
  N0=N0+1
ELSEIF(P1_N_Solid_nodes.EQ.1)THEN
  N1=N1+1
ELSEIF(P1_N_Solid_nodes.EQ.2)THEN
  N2=N2+1
ELSEIF(P1_N_Solid_nodes.EQ.3)THEN
  N3=N3+1
ELSEIF(P1_N_Solid_nodes.EQ.4)THEN
  N4=N4+1
ELSEIF(P1_N_Solid_nodes.EQ.5)THEN
  N5=N5+1
ELSEIF(P1_N_Solid_nodes.EQ.6)THEN
  N6=N6+1
ELSEIF(P1_N_Solid_nodes.EQ.7)THEN
  N7=N7+1
ELSEIF(P1_N_Solid_nodes.EQ.8)THEN
  N8=N8+1
ENDIF


Communication_Dim= abs(i-im1) + abs(j-jm1) + abs(k-km1) !--- Interpolation based on a line (1D), face (2D) or cube (3D) 

IF (Communication_Dim.EQ.1) THEN ! P1 is located on one of the edges of the cube
   phiWall= phiTemp(i,j,k)
   rhoWall= rho(i,j,k)
   N_Case=1
   GOTO 300
ELSEIF (Communication_Dim.EQ.2) THEN ! P1 is located on one of the faces of the cube
   Interpolation_Dim = Node(im1,j,k) + Node(i,jm1,k)+ Node(i,j,km1) - 1
   IF (Interpolation_Dim.EQ.2) THEN
      phiTemp(im1,jm1,km1)= phiTemp(im1,j,k) *(2*abs(im1-i)-1) &
                          + phiTemp(i,jm1,k) *(2*abs(jm1-j)-1) &
                          + phiTemp(i,j,km1) *(2*abs(km1-k)-1) 
      rho(im1,jm1,km1)    = rho(im1,j,k) *(2*abs(im1-i)-1) &
                          + rho(i,jm1,k) *(2*abs(jm1-j)-1) &
                          + rho(i,j,km1) *(2*abs(km1-k)-1) 
      dd=sqrt((P1_x-i)**2 +(P1_y-j)**2 +(P1_z-k)**2) 
      phiwall = ( (sqrt(2.0_dbl)-dd)*phiTemp(i,j,k) + dd*phiTemp(im1,jm1,km1)) /sqrt(2.0_dbl) 
      rhoWall = ( (sqrt(2.0_dbl)-dd)*    rho(i,j,k) + dd*    rho(im1,jm1,km1)) /sqrt(2.0_dbl) 
      N_Case=2
   ELSEIF (Interpolation_Dim.EQ.1) THEN
      phiTemp(im1,jm1,km1)= abs(im1-i)* node(im1,j,k) * phiTemp(im1,j,k) &
                          + abs(jm1-j)* node(i,jm1,k) * phiTemp(i,jm1,k) &
                          + abs(km1-k)* node(i,j,km1) * phiTemp(i,j,km1) 
      rho(im1,jm1,km1)    = abs(im1-i)* node(im1,j,k) * rho(im1,j,k) &
                          + abs(jm1-j)* node(i,jm1,k) * rho(i,jm1,k) &
                          + abs(km1-k)* node(i,j,km1) * rho(i,j,km1) 
      dd=sqrt((P1_x-i)**2 +(P1_y-j)**2 +(P1_z-k)**2) 
      phiWall= ( (sqrt(2.0_dbl)-dd)*phiTemp(i,j,k) + dd*phiTemp(im1,jm1,km1)) /sqrt(2.0_dbl) 
      rhoWall= ( (sqrt(2.0_dbl)-dd)*rho(i,j,k)     + dd*rho(im1,jm1,km1)    ) /sqrt(2.0_dbl) 
      N_Case=3
   ELSEIF (Interpolation_Dim.EQ.0) THEN
      phiWall= phiTemp(i,j,k) 
      rhoWall= rho(i,j,k) 
      N_Case=4
   ENDIF
   GOTO 300
ELSEIF (Communication_Dim.EQ.3) THEN ! P1 is located on one of the faces of the cube

   !===================================================================================================================================================
   !=== Treating the Solid nodes in a cell to assing drug concentration values to them, since all 8 nodes are going to be used 
   !=== for trilinear interpolation of the drug concentration (to the location of the A* at the boundary)
   !===================================================================================================================================================
   opX(ix0)=ix1
   opX(ix1)=ix0
   opY(iy0)=iy1
   opY(iy1)=iy0
   opZ(iz0)=iz1
   opZ(iz1)=iz0

   P1_N_Solid_nodes = node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1) 

   !--- 1 Solid node: Case A --------------------------------------------------------------------------------------------------------------------------
   IF (P1_N_Solid_nodes .EQ. 1) THEN
      DO ii=ix0,ix1
         DO jj=iy0,iy1
            DO kk=iz0,iz1
               IF (node(ii,jj,kk).EQ.1) THEN 
                   phiTemp(ii,jj,kk) =(2.0_dbl*(phiTemp(opX(ii),jj,kk) + phiTemp(ii,opY(jj),kk) + phiTemp(ii,jj,opZ(kk))) & 
                                     -(phiTemp(opX(ii),opY(jj),kk) + phiTemp(ii,opY(jj),opZ(kk)) + phiTemp(opX(ii),jj,opZ(kk))))/3.0_dbl
                   rho(ii,jj,kk)     =(2.0_dbl*(rho(opX(ii),jj,kk) + rho(ii,opY(jj),kk) + rho(ii,jj,opZ(kk))) & 
                                     -(rho(opX(ii),opY(jj),kk) + rho(ii,opY(jj),opZ(kk)) + rho(opX(ii),jj,opZ(kk))))/3.0_dbl
                   N_Case=5
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   !--- 2 Solid nodes:  Cases B, C, D -------------------------------------------------------------------------------------------------------------------
   IF (P1_N_Solid_nodes .EQ. 2) THEN
      DO ii=ix0,ix1
         DO jj=iy0,iy1
            DO kk=iz0,iz1
               IF (node(ii,jj,kk).EQ.1) THEN 
                  !--- face1
                  key1= (1.0_lng-node(opX(ii),jj,kk)) *(1.0_lng-node(ii,opY(jj),kk)) *(1.0_lng-node(opX(ii),opY(jj),kk)) 
                  phi1= phiTemp(opX(ii),jj,kk) + phiTemp(ii,opY(jj),kk) - phiTemp(opX(ii),opY(jj),kk)
                  rho1=     rho(opX(ii),jj,kk) +     rho(ii,opY(jj),kk) -     rho(opX(ii),opY(jj),kk)
                  !--- face2
                  key2= (1.0_lng-node(opX(ii),jj,kk)) *(1.0_lng-node(ii,jj,opZ(kk))) *(1.0_lng-node(opX(ii),jj,opZ(kk))) 
                  phi2= phiTemp(opX(ii),jj,kk) + phiTemp(ii,jj,opZ(kk)) - phiTemp(opX(ii),jj,opZ(kk))
                  rho2=     rho(opX(ii),jj,kk) +     rho(ii,jj,opZ(kk)) -     rho(opX(ii),jj,opZ(kk))
                  !--- face3
                  key3= (1.0_lng-node(ii,opY(jj),kk)) *(1.0_lng-node(ii,jj,opZ(kk))) *(1.0_lng-node(ii,opY(jj),opZ(kk))) 
                  phi3= phiTemp(ii,opY(jj),kk) + phiTemp(ii,jj,opZ(kk)) - phiTemp(ii,opY(jj),opZ(kk))
                  rho3=     rho(ii,opY(jj),kk) +     rho(ii,jj,opZ(kk)) -     rho(ii,opY(jj),opZ(kk))

                  phiTemp(ii,jj,kk) = (key1*phi1+key2*phi2+key3*phi3)/(key1+key2+key3)
                  rho(ii,jj,kk)     = (key1*rho1+key2*rho2+key3*rho3)/(key1+key2+key3)
                  N_Case=6
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   !--- 3,4,5,6 Solid nodes: Cases E-U ------------------------------------------------------------------------------------------------------------------
   IF ((P1_N_Solid_nodes .EQ. 3).OR.(P1_N_Solid_nodes .EQ. 4).OR.(P1_N_Solid_nodes .EQ. 5).OR.(P1_N_Solid_nodes .EQ. 6))  THEN
      nnode(ix0,iy0,iz0)= node(ix0,iy0,iz0)
      nnode(ix1,iy0,iz0)= node(ix1,iy0,iz0)
      nnode(ix0,iy1,iz0)= node(ix0,iy1,iz0)
      nnode(ix0,iy0,iz1)= node(ix0,iy0,iz1)
      nnode(ix1,iy1,iz0)= node(ix1,iy1,iz0)
      nnode(ix1,iy0,iz1)= node(ix1,iy0,iz1)
      nnode(ix0,iy1,iz1)= node(ix0,iy1,iz1)
      nnode(ix1,iy1,iz1)= node(ix1,iy1,iz1)
100   Counter1= 0
      DO ii=ix0,ix1
         DO jj=iy0,iy1
            DO kk=iz0,iz1
               IF (nnode(ii,jj,kk).EQ.1) THEN
                  !--- face1 ------------------------------------------------------------------------
                  key1= (1-nnode(opX(ii),jj,kk))*(1-nnode(ii,opY(jj),kk))*(1-nnode(opX(ii),opY(jj),kk)) 
                  phi1=  phiTemp(opX(ii),jj,kk) + phiTemp(ii,opY(jj),kk) - phiTemp(opX(ii),opY(jj),kk)
                  rho1=      rho(opX(ii),jj,kk) +     rho(ii,opY(jj),kk) -     rho(opX(ii),opY(jj),kk)
                  !--- face2 ------------------------------------------------------------------------
                  key2= (1-nnode(opX(ii),jj,kk))*(1-nnode(ii,jj,opZ(kk)))*(1-nnode(opX(ii),jj,opZ(kk))) 
                  phi2=  phiTemp(opX(ii),jj,kk) + phiTemp(ii,jj,opZ(kk)) - phiTemp(opX(ii),jj,opZ(kk))
                  rho2=      rho(opX(ii),jj,kk) +     rho(ii,jj,opZ(kk)) -     rho(opX(ii),jj,opZ(kk))
                  !--- face3 ------------------------------------------------------------------------
                  key3= (1-nnode(ii,opY(jj),kk))*(1-nnode(ii,jj,opZ(kk)))*(1-nnode(ii,opY(jj),opZ(kk))) 
                  phi3=  phiTemp(ii,opY(jj),kk) + phiTemp(ii,jj,opZ(kk)) - phiTemp(ii,opY(jj),opZ(kk))
                  rho3=      rho(ii,opY(jj),kk) +     rho(ii,jj,opZ(kk)) -    rho(ii,opY(jj),opZ(kk))

                  IF ((key1+key2+key3) .GT. 0) THEN !--- At least one face of the cube connected to this node can be used for bilinear interpolation of drug concentration
                     phiTemp(ii,jj,kk) = (key1*phi1+key2*phi2+key3*phi3)/(key1+key2+key3) !--- Concentration is assigned to this solid node
                     rho(ii,jj,kk)     = (key1*rho1+key2*rho2+key3*rho3)/(key1+key2+key3) !--- Concentration is assigned to this solid node
                     nnode(ii,jj,kk) = 0 !--- This node is marked as non-solid in case any other node needs its value for interpolation       
                  ELSE !--- None of the 3 faces connected to this node can be used for bilinear interpolation of teh drug concentration  
                     Counter1=Counter1+1  ! Number of nodes in the cell that cannot be dealt with using linear interpolation
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF ((Counter1.GT.0).AND.(Counter1.LT.4)) THEN
         N_Case=7
         GOTO 100
      ELSEIF (Counter1.GE.4) THEN
         N_Case=8
200      Counter2=0
         DO ii=ix0,ix1
            DO jj=iy0,iy1
               DO kk=iz0,iz1
                  IF (nnode(ii,jj,kk).EQ.1) THEN
                     phi_N   = (1-nnode(opX(ii),jj,kk)) + (1-nnode(ii,opY(jj),kk)) +(1-nnode(ii,jj,opZ(kk)))
                     IF (phi_N .GT.0) THEN
                        phi_sum = (1-nnode(opX(ii),jj,kk)) * phiTemp(opX(ii),jj,kk) &
                                + (1-nnode(ii,opY(jj),kk)) * phiTemp(ii,opY(jj),kk) &
                                + (1-nnode(ii,jj,opZ(kk))) * phiTemp(ii,jj,opZ(kk))
                        rho_sum = (1-nnode(opX(ii),jj,kk)) *     rho(opX(ii),jj,kk) &
                                + (1-nnode(ii,opY(jj),kk)) *     rho(ii,opY(jj),kk) &
                                + (1-nnode(ii,jj,opZ(kk))) *     rho(ii,jj,opZ(kk))
                        phiTemp(ii,jj,kk) = phi_sum/phi_N
                            rho(ii,jj,kk) = rho_sum/phi_N
                        nnode(ii,jj,kk) = 0
                     ELSE
                        Counter2=Counter2+1  
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         IF (Counter2.GT.0) THEN
            GOTO 200
            N_Case=9
         ENDIF   
      ENDIF   
   ENDIF
   !--- 7 Solid nodes: Case V ------------------------------------------------------------------------------------------------------------------------------------
   !IF (P1_N_Solid_nodes .EQ. 7)  THEN !--- Only one fluid node exist in this cube --> concentration at P1 set equal to concentration at the fluid node 

   IF (P1_N_Solid_nodes .EQ. 7)  THEN !--- Only one fluid node exist in this cube --> concentration at P1 set equal to concentration at the fluid node 
      DO ii=ix0,ix1
         DO jj=iy0,iy1
            DO kk=iz0,iz1
               IF (node(ii,jj,kk).EQ.1) THEN 
                  phiWall= (1-node(opX(ii),jj,kk))           * phiTemp(opX(ii),jj,kk)      &
                         + (1-node(ii,opY(jj),kk))           * phiTemp(ii,opY(jj),kk)      &
                         + (1-node(ii,jj,opZ(kk)))           * phiTemp(ii,jj,opZ(kk))      &
                         + (1-node(opX(ii),opY(jj),kk))      * phiTemp(opX(ii),opY(jj),kk) &
                         + (1-node(opX(ii),jj,opZ(kk)))      * phiTemp(opX(ii),jj,opZ(kk)) &
                         + (1-node(ii,opY(jj),opZ(kk)))      * phiTemp(ii,opY(jj),opZ(kk)) &
                         + (1-node(opX(ii),opY(jj),opZ(kk))) * phiTemp(opX(ii),opY(jj),opZ(kk)) 
                  rhoWall= (1-node(opX(ii),jj,kk))           *     rho(opX(ii),jj,kk)      &
                         + (1-node(ii,opY(jj),kk))           *     rho(ii,opY(jj),kk)      &
                         + (1-node(ii,jj,opZ(kk)))           *     rho(ii,jj,opZ(kk))      &
                         + (1-node(opX(ii),opY(jj),kk))      *     rho(opX(ii),opY(jj),kk) &
                         + (1-node(opX(ii),jj,opZ(kk)))      *     rho(opX(ii),jj,opZ(kk)) &
                         + (1-node(ii,opY(jj),opZ(kk)))      *     rho(ii,opY(jj),opZ(kk)) &
                         + (1-node(opX(ii),opY(jj),opZ(kk))) *     rho(opX(ii),opY(jj),opZ(kk)) 
                   N_Case=10
                   GOTO 300
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDIF   
!===================================================================================================================================================
!=== Trilinear interpolation of drug concentration to the location of A* (phiWall)
!===================================================================================================================================================
IF (ix1 /= ix0) THEN
   xd= (P1_x-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
ELSE
   xd= 0.0_dbl
END IF
IF (iy1 /= iy0) THEN
   yd= (P1_y-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
ELSE
   yd= 0.0_dbl
END IF
IF (iz1 /= iz0) THEN
   zd= (P1_z-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
ELSE
   zd= 0.0_dbl
END IF
!--- Interpolation in x-direction
c00= phiTemp(ix0,iy0,iz0) * (1.0000_dbl-xd) + phiTemp(ix1,iy0,iz0) * xd
c01= phiTemp(ix0,iy0,iz1) * (1.0000_dbl-xd) + phiTemp(ix1,iy0,iz1) * xd
c10= phiTemp(ix0,iy1,iz0) * (1.0000_dbl-xd) + phiTemp(ix1,iy1,iz0) * xd
c11= phiTemp(ix0,iy1,iz1) * (1.0000_dbl-xd) + phiTemp(ix1,iy1,iz1) * xd
!--- Interpolation in y-direction
c0 = c00 * (1.0000_dbl-yd) + c10 * yd
c1 = c01 * (1.0000_dbl-yd) + c11 * yd
!--- Interpolation in z-direction
phiWall = c0 * (1.0000_dbl-zd) + c1 * zd
!--- Interpolation in x-direction
c00= rho(ix0,iy0,iz0) * (1.0000_dbl-xd) + rho(ix1,iy0,iz0) * xd
c01= rho(ix0,iy0,iz1) * (1.0000_dbl-xd) + rho(ix1,iy0,iz1) * xd
c10= rho(ix0,iy1,iz0) * (1.0000_dbl-xd) + rho(ix1,iy1,iz0) * xd
c11= rho(ix0,iy1,iz1) * (1.0000_dbl-xd) + rho(ix1,iy1,iz1) * xd
!--- Interpolation in y-direction
c0 = c00 * (1.0000_dbl-yd) + c10 * yd
c1 = c01 * (1.0000_dbl-yd) + c11 * yd
!--- Interpolation in z-direction
rhoWall = c0 * (1.0000_dbl-zd) + c1 * zd

!N_Case=11


!===================================================================================================================================================
!=== Treat the Neumann BC similar to the way we treated Dirichlet BC now that the phi_wall is defined
!===================================================================================================================================================
!----- neighboring node (fluid side) ---------------------------------------------------------------
300 ip1 = i + ex(m)
    jp1 = j + ey(m)
    kp1 = k + ez(m)
WRITE(*,*) 'phiWall,N_Case,i,j,k,im1,jm1,km1,q,iter',phiWall,N_Case,i,j,k,im1,jm1,km1,q,iter
!------ This rarely happens (both neighboring nodes over a line are solid)
!------  use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF

!----- Computing values at A* & the scalar streamed from A* (Chpter 3 paper) -----------------------
phiAstar= phiWall                                                  ! phi at solid surface
rhoAstar= rhoWall                                                  ! phi at solid surface
!rhoAstar= (rho(i,j,k)-rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)		 ! Extrapolate density
CALL Equilibrium_LOCAL(m,rhoAstar,ub,vb,wb,feq_Astar)    		       ! f_eq in mth direction
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar		             ! Contribution from A* to B*  

!----- Using only A and A* for interpolation (instead of A* and B*) 
PkA= (fplus(m,i,j,k)/rho(i,j,k) - wt(m)*Delta)*phiTemp(i,j,k)      ! contribution from current node to next in the mth direction

IF(q .LT. 0.25) THEN
  q = 0.25_dbl
END IF
phiBC	= ((PkAstar - PkA)/q) + PkA	
!===================================================================================================
END SUBROUTINE BC_Scalar_Permeability_1storder
!===================================================================================================





!===================================================================================================
SUBROUTINE BC_Scalar_Permeability(m,i,j,k,im1,jm1,km1,phiBC)				! implements the scalar BCs 
!===================================================================================================
IMPLICIT NONE

REAL(dbl),   INTENT(OUT):: phiBC           				! scalar contribution from the boundary condition
INTEGER(lng),INTENT(IN) :: m,i,j,k,im1,jm1,km1		! index variables
INTEGER(lng) :: ix0,ix1,iy0,iy1,iz0,iz00,iz1,iz11	! Trilinear interpolation parameters
INTEGER(lng) :: P1_N_Solid_nodes, P2_N_Solid_nodes
INTEGER(lng) :: ip1,jp1,kp1   	                  ! First neighboring node location
INTEGER(lng) :: xaxis,yaxis
REAL(dbl)    :: c00,c01,c10,c11,c0,c1,c,xd,yd,zd  ! Trilinear interpolation parameters
REAL(dbl)    :: xt,yt,zt,rt,vt,q                  ! Location of the boundary between i,j,k node and im1,jm1,km1 node
REAL(dbl)    :: Geom_nx,Geom_ny,Geom_nz,Geom_n_mag
REAL(dbl)    :: Ax,Ay,Az,A_mag
REAL(dbl)    :: Prod
REAL(dbl)    :: rhoAstar, phiAstar, feq_Astar,  PkAstar,PkA	! density at boundary and contribution of scalar from boundary
REAL(dbl)    :: rhoBstar, phiBstar, fPlusBstar, PkBstar 		! Values interpolated to Bstar location
REAL(dbl)    :: cosTheta, sinTheta					 
REAL(dbl)    :: ub, vb, wb
REAL(dbl)    :: P1_x,P1_y,P1_z,P2_x,P2_y,P2_z
REAL(dbl)    :: P1_phi,P2_phi,phiWall_new, Del_phiWall,DphiDn
REAL(dbl)    :: alpha,Diffusivity,Pww

CALL qCalc_iter(m,i,j,k,im1,jm1,km1,xt,yt,zt,rt,q)

IF (Flag_Couette) THEN
   ub= 0.0_dbl
   vb= 0.0_dbl
   xaxis= ANINT(0.5_dbl*(nx+1))
   IF ((iMin-1+im1) .GT. xaxis) THEN 
      wb = vel(k)
   ELSE 
      wb= -vel(k)
   END IF
ELSE
   cosTheta= xt/rt
   sinTheta= yt/rt
   IF (k.NE.km1) THEN
      vt = ((zt-z(k))*vel(km1)+(z(km1)-zt)*vel(k))/(z(km1)-z(k))
   ELSE
      vt = (vel(k)+vel(km1))*0.5_dbl
   ENDIF
   ub = vt* cosTheta						
   vb = vt* sinTheta					
   wb = -s_movingF/vcf 		
END IF

!--- Computing the normal to the geometry based on the equation for boundary 
Geom_nx= -xt/rt
Geom_ny= -yt/rt
Geom_nz= -amp1 *(2.0_dbl*PI/lambda1) *SIN(PI+(2.0_dbl *PI*zt/lambda1)) 
!---normalizing the geometry normal vector 
Geom_n_mag=sqrt(Geom_nx**2.0_dbl + Geom_ny**2.0_dbl + Geom_nz**2.0_dbl)
Geom_nx=Geom_nx/Geom_n_mag
Geom_ny=Geom_ny/Geom_n_mag
Geom_nz=Geom_nz/Geom_n_mag

!--------------------------------------------------------------------------------------------------
!--- Finding location of the point, P1, which is one mesh size away from (xt,yt,zt) at the boundary
alpha=1.5_dbl               ! The coefficient which defines the distance to walk away from the boundary
xaxis=ANINT(0.5_dbl*(nx+1))
yaxis=ANINT(0.5_dbl*(ny+1))
P1_x= ((xt + alpha*Geom_nx*xcf)/xcf) - iMin + xaxis + 1
P1_y= ((yt + alpha*Geom_ny*xcf)/xcf) - jMin + yaxis + 1
P1_z= ((zt + alpha*Geom_nz*xcf)/xcf) - kMin + 1

ix0= FLOOR(P1_x)
ix1= CEILING(P1_x)
iy0= FLOOR(P1_y)
iy1= CEILING(P1_y)
iz0= FLOOR(P1_z)
iz1= CEILING(P1_z)

P1_N_Solid_nodes =   node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1) 

!IF (P1_N_Solid_nodes .GT.0) THEN
!   WRITE(*,*) 'iter:',iter,'-----------------------'
!   WRITE(*,*) 'node: i,j,k,m:',i,j,k,m
!   WRITE(*,*) 'Geom_n:',Geom_nx,Geom_ny,Geom_nz
!   WRITE(*,*) 'P1: x,y,z,N',P1_x,P1_y,P1_z,P1_N_Solid_nodes 
!   WRITE(*,*) 'P1: ix,iy,iz:',ix0,ix1,iy0,iy1,iz0,iz1
!   WRITE(*,*) 'Status', node(ix0,iy0,iz0),node(ix1,iy0,iz0),node(ix0,iy1,iz0),node(ix0,iy0,iz1),node(ix1,iy1,iz0),node(ix1,iy0,iz1),node(ix0,iy1,iz1),node(ix1,iy1,iz1) 
!ENDIF

IF (ix1 /= ix0) THEN
   xd= (P1_x-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
ELSE
   xd= 0.0_dbl
END IF
IF (iy1 /= iy0) THEN
   yd= (P1_y-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
ELSE
   yd= 0.0_dbl
END IF
IF (iz1 /= iz0) THEN
   zd= (P1_z-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
ELSE
   zd= 0.0_dbl
END IF
!--- Interpolation in x-direction
c00= phiTemp(ix0,iy0,iz0) * (1.0000_dbl-xd) + phiTemp(ix1,iy0,iz0) * xd
c01= phiTemp(ix0,iy0,iz1) * (1.0000_dbl-xd) + phiTemp(ix1,iy0,iz1) * xd
c10= phiTemp(ix0,iy1,iz0) * (1.0000_dbl-xd) + phiTemp(ix1,iy1,iz0) * xd
c11= phiTemp(ix0,iy1,iz1) * (1.0000_dbl-xd) + phiTemp(ix1,iy1,iz1) * xd
!--- Interpolation in y-direction
c0 = c00 * (1.0000_dbl-yd) + c10 * yd
c1 = c01 * (1.0000_dbl-yd) + c11 * yd
!--- Interpolation in z-direction
P1_phi = c0 * (1.0000_dbl-zd) + c1 * zd

!--------------------------------------------------------------------------------------------------
!--- Finding location of the point, P2, which is two mesh size away from (xt,yt,zt) at the boundary
P2_x= ((xt + 2.0_dbl*alpha*Geom_nx*xcf)/xcf) - iMin + xaxis + 1
P2_y= ((yt + 2.0_dbl*alpha*Geom_ny*xcf)/xcf) - jMin + yaxis + 1
P2_z= ((zt + 2.0_dbl*alpha*Geom_nz*xcf)/xcf) - kMin + 1

ix0= FLOOR(P2_x)
ix1= CEILING(P2_x)
iy0= FLOOR(P2_y)
iy1= CEILING(P2_y)
iz0= FLOOR(P2_z)
iz1= CEILING(P2_z)
P2_N_Solid_nodes = node(ix0,iy0,iz0)+node(ix1,iy0,iz0)+node(ix0,iy1,iz0)+node(ix0,iy0,iz1)+node(ix1,iy1,iz0)+node(ix1,iy0,iz1)+node(ix0,iy1,iz1)+node(ix1,iy1,iz1) 

IF (ix1 /= ix0) THEN
   xd= (P2_x-REAL(ix0,dbl))/(REAL(ix1,dbl)-REAL(ix0,dbl))
ELSE
   xd= 0.0_dbl
END IF
IF (iy1 /= iy0) THEN
   yd= (P2_y-REAL(iy0,dbl))/(REAL(iy1,dbl)-REAL(iy0,dbl))
ELSE
   yd= 0.0_dbl
END IF
IF (iz1 /= iz0) THEN
   zd= (P2_z-REAL(iz0,dbl))/(REAL(iz1,dbl)-REAL(iz0,dbl))
ELSE
   zd= 0.0_dbl
END IF
!--- Interpolation in x-direction
c00= phiTemp(ix0,iy0,iz0) * (1.0_dbl-xd) + phiTemp(ix1,iy0,iz0) * xd
c01= phiTemp(ix0,iy0,iz1) * (1.0_dbl-xd) + phiTemp(ix1,iy0,iz1) * xd
c10= phiTemp(ix0,iy1,iz0) * (1.0_dbl-xd) + phiTemp(ix1,iy1,iz0) * xd
c11= phiTemp(ix0,iy1,iz1) * (1.0_dbl-xd) + phiTemp(ix1,iy1,iz1) * xd
!--- Interpolation in y-direction
c0 = c00 * (1.0_dbl-yd) + c10 * yd
c1 = c01 * (1.0_dbl-yd) + c11 * yd
!--- Interpolation in z-direction
P2_phi = c0 * (1.0_dbl-zd) + c1 * zd

IF (Flag_2step_Permeability) THEN
   Pww=1.0e-16  !No flux at this stage. Later in PassiveSacalr the preamibility is applied
ELSE
   Pww=Pw
ENDIF 
Diffusivity=((nuL/Sc)*(xcf**2.0_dbl)/tcf)*10000.0_dbl   ! Diffusivity used in LBM [cm2/s]
DphiDn=      0.0_dbl
Del_phiWall= 1.0_dbl
phiWall=     (4.0_dbl*P1_phi - 1.0_dbl*P2_phi -2.0_dbl*DphiDn) / 3.0_dbl 
DO WHILE (Del_phiWall .GE. 1.0e-10) 
   DphiDn= (Pww*phiWall/Diffusivity)*(100.0*xcf)
   phiWall_new= (4.0_dbl*P1_phi - 1.0_dbl*P2_phi -2.0_dbl*DphiDn) / 3.0_dbl 
   Del_phiWall =abs(phiWall_new-phiWall)
   phiWall=phiWall_new
ENDDO
!----- neighboring node (fluid side) ---------------------------------------------------------------
ip1 = i + ex(m)
jp1 = j + ey(m)
kp1 = k + ez(m)
!------ This rarely happens (both neighboring nodes over a line are solid)
!------  use values from the current node as an approximation
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF

!----- Computing values at A* & the scalar streamed from A* (Chpter 3 paper) -----------------------
rhoAstar= (rho(i,j,k)-rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)		 ! Extrapolate density
CALL Equilibrium_LOCAL(m,rhoAstar,ub,vb,wb,feq_Astar)    		       ! f_eq in mth direction
phiAstar= phiWall                                                  ! phi at solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar		             ! Contribution from A* to B*  

!----- Using only A and A* for interpolation (instead of A* and B*) 
PkA= (fplus(m,i,j,k)/rho(i,j,k) - wt(m)*Delta)*phiTemp(i,j,k)      ! contribution from current node to next in the mth direction

IF(q .LT. 0.25) THEN
  q = 0.25_dbl
END IF
phiBC	= ((PkAstar - PkA)/q) + PkAstar	
!===================================================================================================
END SUBROUTINE BC_Scalar_Permeability
!===================================================================================================









!===================================================================================================
SUBROUTINE AbsorbedScalarS(i,j,k,m,im1,jm1,km1,phiBC) 	  ! Monitoring the abosrption at boundaries	
!===================================================================================================
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: i,j,k,m,im1,jm1,km1			! index variables
REAL(dbl),    INTENT(IN) :: phiBC   		  		! scalar contribution from the boundary condition

INTEGER(lng) :: ip1,jp1,kp1
INTEGER(lng) :: xaxis
REAL(dbl)    :: phiOUT, phiIN					! scalar values exchanged with the wall
REAL(dbl)    :: phiAbsorbedSleft,phiAbsorbedSright
REAL(dbl)    :: phiINleft, phiINright, phiOUTright,phiOUTleft
REAL(dbl)    :: feq_AO_u0
REAL(dbl)    :: rhoAstar,phiAstar, PkAstar,feq_Astar,feq_Bstar
REAL(dbl)    :: rhoA, phiA, feq_A, PkA
REAL(dbl)    :: fPlusBstar, rhoBstar, phiBstar, PkBstar
REAL(dbl)    :: ub,vb,wb, ubb,vbb,wbb
REAL(dbl)    :: cosTheta, sinTheta                              ! COS(theta), SIN(theta)
REAL(dbl)    :: q
REAL(dbl)    :: xt,yt,zt,rt,vt			             	! boundary coordinate

CALL qCalc_iter(m,i,j,k,im1,jm1,km1,xt,yt,zt,rt,q)

IF (Flag_Couette) THEN
   ub= 0.0_dbl
   vb= 0.0_dbl
   xaxis= ANINT(0.5_dbl*(nx+1))
   IF ((iMin-1+im1) .GT. xaxis) THEN 
      wb = vel(k)
   ELSE 
      wb= -vel(k)
   END IF
ELSE
   cosTheta= xt/rt
   sinTheta= yt/rt
   IF (k.NE.km1) THEN
      vt = ((zt-z(k))*vel(km1)+(z(km1)-zt)*vel(k))/(z(km1)-z(k))
   ELSE
      vt = (vel(k)+vel(km1))*0.5_dbl
   ENDIF
   ub = vt* cosTheta						! x-component of the velocity at i,j,k
   vb = vt* sinTheta						! y-component of the velocity at i,j,k
   wb = -s_movingF/vcf  				! no z-component in this case)
END IF

ubb= u(i,j,k)- ub
vbb= v(i,j,k)- vb
wbb= w(i,j,k)- wb

!---------------------------------------------------------------------------------------------------
!----- Computing phiOUT ----------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
CALL Equilibrium_LOCAL(bb(m),rho(i,j,k),ubb,vbb,wbb,feq_AO_u0)
phiOUT= (feq_AO_u0/rho(i,j,k) - wt(bb(m))*Delta)*phiTemp(i,j,k)

!---------------------------------------------------------------------------------------------------
!---- Conmputing phiIN------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------
!----- neighboring node (fluid side)	
ip1 = i + ex(m) 			
jp1 = j + ey(m)			
kp1 = k + ez(m)		
IF(node(ip1,jp1,kp1) .NE. FLUID) THEN
  ip1 = i
  jp1 = j
  kp1 = k
END IF	

!----- Computing values at A* & scalar streamed from A* (Chpter 3 paper)
rhoAstar= (rho(i,j,k)- rho(ip1,jp1,kp1))*(1+q)+ rho(ip1,jp1,kp1)	! extrapolate the density
CALL Equilibrium_LOCAL(m,rhoAstar,0.0_dbl,0.0_dbl,0.0_dbl,feq_Astar)	! equibrium distribution function in mth direction, velocity relative to the boundary is zero
phiAstar= phiWall							! getting phi at the solid surface
PkAstar= (feq_Astar/rhoAstar- wt(m)*Delta)*phiAstar			! contribution from the wall in mth direction (0 if phiWall=0)

!------ Computing values at B* & scalar streamed from B* (Chpter 3 paper)
!rhoBstar=   (1-q)*rho(ip1,jp1,kp1)     + q*rho(i,j,k)
!CALL Equilibrium_LOCAL(m,rhoBstar,ubb,vbb,wbb,feq_Bstar)
!phiBstar=   (1-q)*phiTemp(ip1,jp1,kp1) + q*phiTemp(i,j,k)
!PkBstar=    (feq_Bstar/rhoBstar - wt(m)*Delta)*phiBstar

!phiIN= PkAstar+ (PkAstar- PkBstar)*(1-q)


!---- Modification for moving boundary in case of using only A and A* for BC
rhoA= rho(i,j,k)
phiA= phiTemp(i,j,k)
CALL Equilibrium_LOCAL(m,rhoA,ubb,vbb,wbb,feq_A) 			! Velocity relative to boundary is used
PkA= (feq_A/rhoA - wt(m)*Delta)* phiA 
IF(q .LT. 0.25) THEN
  q = 0.25_dbl
END IF 

phiIN   = ((PkAstar - PkA)/q) + PkAstar

!--- No Modifications in book-keeping for moving boundaries
!phiIN= phiBC                                                    	 ! contribution from wall to crrent node (in)
!phiOUT= (fplus(bb(m),i,j,k)/rho(i,j,k)-wt(bb(m))*Delta)*phiTemp(i,j,k)

phiAbsorbedS = phiAbsorbedS + (phiOUT-phiIN)				! scalar absorbed at current location in mth direction
!===================================================================================================
END SUBROUTINE AbsorbedScalarS
!===================================================================================================






!===================================================================================================
END MODULE BCscalar
!===================================================================================================
