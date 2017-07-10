!===================================================================================================
! This subroutine randomly creates the particle's locations  
!===================================================================================================
USE SetPrecision
USE Setup
USE Geometry

IMPLICIT NONE

REAL(dbl) :: D_Couette, L_Couette, Dx_Couette, Dy_Couette, Dz_Couette
REAL(dbl) :: x_center, y_center, z_center
REAL(dbl) :: rMax, teta1Max, teta2Max, rr, teta1, teta2
REAL(dbl) :: x_particle, y_particle, z_particle
REAL(dbl) :: R_Particle, R_Par_Max, D_Par_Max, R_Boundary
REAL(dbl) :: R_left, R_right, dz, Volume, Area
REAL(dbl) :: xp,yp,zp,up,vp,wp,U_slip,sh_conf,Sh_shear,Sh_slip,rp,bulk_conc,delNBbyCV,SSt,S,par_Conc 
REAL(dbl) :: eps,S_ratio,Sh,tmp,Cb,deltaR,zcf3,Drug_Released
INTEGER   :: n,num_Par, i, j, k, z_left, z_right, CPU, parid, cur_part
CHARACTER(7):: TEMP,iter_char                                       ! iteration stored as a character
CHARACTER(5):: sub_char,tmp_char
TYPE(ParRecord), POINTER :: CurPar
TYPE(ParRecord), POINTER :: currentt
TYPE(ParRecord), POINTER :: nextt

CALL ReadInput

nuL  = (2.0_dbl*tau -1.0_dbl)/6.0_dbl	! lattice kinematic viscosity
denL = 1.0_dbl		                	  ! arbitrary lattice density (1.0 for convenience)
kMin = 1
kMax = nz
nzSub= nz
numprocs=1
CPU=1
S_ratio =2.30196707 ! For no buffer case only
zcf=L/nz
zcf3=zcf*zcf*zcf

CALL Global_Setup
CALL AdvanceGeometry
!r(0:nzSub+1) = rDom(kMin-1:kMax+1)   ! Fill out the local radius array
!Volume=0.0
!DO i=1,nz
!   Volume= Volume+ xcf*PI*(r(i)**2.0)
!   Area  = Area  + xcf*PI* r(i) *2.0
!END DO

Write(*,*) 'Please enter the number of particles:'
read(*,*) np
Write(*,*) 'Please enter the ieration:'
read(*,*) iter
Write(*,*) 'Please enter the Released Drug:'
read(*,*) Drug_Released

WRITE(iter_char(1:7),'(I7.7)') iter
WRITE(sub_char(1:5),'(I5.5)') CPU
OPEN(160,FILE='pardat-'//iter_char//'-'//sub_char//'.csv')

CALL list_init(ParListHead)
CurPar => ParListHead

READ(160,*) tmp_char
DO i = 1, np
  write(*,*) i 
   READ(160,*) xp,yp,zp,up,vp,wp,parid,sh_conf,Sh_shear,Sh_slip,rp,bulk_conc,delNBbyCV,SSt,S,par_Conc,cur_part 
   CALL list_init(CurPar%next)		
   CurPar%next%prev => CurPar
   CurPar%next%next => null()      
   CurPar%next%pardata%parid = parid
   CurPar%next%pardata%xp = xp
   CurPar%next%pardata%yp = yp
   CurPar%next%pardata%zp = zp
   CurPar%next%pardata%up = up/vcf
   CurPar%next%pardata%vp = vp/vcf
   CurPar%next%pardata%wp = wp/vcf
   CurPar%next%pardata%rp = rp
   CurPar%next%pardata%xpold = 0.0_dbl 
   CurPar%next%pardata%ypold = 0.0_dbl
   CurPar%next%pardata%zpold = 0.0_dbl
   CurPar%next%pardata%upold = 0.0_dbl 
   CurPar%next%pardata%vpold = 0.0_dbl 
   CurPar%next%pardata%wpold = 0.0_dbl 
   CurPar%next%pardata%rpold = 0.0_dbl 
   CurPar%next%pardata%par_conc = S_ratio* S_intrinsic 
   CurPar%next%pardata%gamma_cont = 0.0000_dbl
   CurPar%next%pardata%sh_conf  = Sh_conf
   CurPar%next%pardata%sh_shear = Sh_shear
   CurPar%next%pardata%sh_slip  = Sh_slip 
   CurPar%next%pardata%S        = S
   CurPar%next%pardata%Sst      = Sst
   CurPar%next%pardata%Veff     = 0.0_dbl
   CurPar%next%pardata%Nbj      = 0.0_dbl
   CurPar%next%pardata%bulk_conc= bulk_conc* S_intrinsic
   CurPar%next%pardata%delNBbyCV= delNBbycv
   CurPar%next%pardata%cur_part = 1
   CurPar%next%pardata%new_part = 1
   CurPar => CurPar%next
END DO
CLOSE(160)

!--- Shear ------------------
S=0.0_dbl
currentt => ParListHead%next
DO WHILE (ASSOCIATED(currentt))
   nextt => currentt%next
   S= S + currentt%pardata%S
   currentt => nextt
ENDDO

!--- Cb --------------------
Cb=0.0_dbl
!---------------------------
open(5,file='Drug-Conservation-00001-Extrapolated.dat')
Do i=iter, 30000000
   np=0
   iter=iter+1
   currentt => ParListHead%next
   DO WHILE (ASSOCIATED(currentt))
      nextt => currentt%next
      IF (currentt%pardata%rp.GT.1e-16) THEN
         np= np+1
         currentt%pardata%rpold = currentt%pardata%rp
         S_ratio =2.30196707
         currentt%pardata%par_conc = S_ratio * S_intrinsic

         !---- Shear Effects -------------------------------------------------------------
         S= currentt%pardata%S
         Sst= S* (currentt%pardata%rp**2.0) / diffm
         currentt%pardata%Sst= Sst
         IF (Sst .LE. 5.0_dbl) THEN
            currentt%pardata%sh_shear = 1.0_dbl + 0.281_dbl*(Sst**0.5_dbl)          -1.0_dbl
         ELSE IF ((Sst .GT. 5.0_dbl).AND.(Sst .LE. 80.0)) THEN     
            currentt%pardata%sh_shear = 1.181_dbl*(Sst**0.2_dbl)                    -1.0_dbl
         ELSE IF (Sst.GT.80.0) THEN
            currentt%pardata%sh_shear = 4.5_dbl - (7.389/(Sst**(1.0_dbl/3.0_dbl)) ) -1.0_dbl 
         ENDIF
         !---- SLip Effects----------------------------------------------------------------

         !--- Release
         Sh= 1.0 + currentt%pardata%sh_shear  
         eps= (molarvol*diffm*tcf*Sh) * (currentt%pardata%par_conc-Cb)
         tmp= (currentt%pardata%rpold)**2 - 4*eps 
         If (tmp .LT. 0.0_dbl)THEN
            tmp=0.0_dbl
         ENDIF
         currentt%pardata%rp= 0.5_dbl * (currentt%pardata%rpold+ sqrt(tmp))
         deltaR= currentt%pardata%rpold-currentt%pardata%rp
         currentt%pardata%delNBbyCV = (4.0_dbl/3.0_dbl) * PI*(currentt%pardata%rpold**3.0_dbl - currentt%pardata%rp**3.0_dbl) /(molarvol*zcf3)
         Drug_Released= Drug_Released + currentt%pardata%delNBbyCV *1000000.0_dbl*zcf3
!        write(*,*) currentt%pardata%parid,Sh,currentt%pardata%par_conc,molarvol,diffm,tcf,Drug_Released 
!        write(*,*) currentt%pardata%parid,Sh,currentt%pardata%par_conc,molarvol,diffm,tcf,Drug_Released 
     ENDIF
     currentt => nextt
  ENDDO    
  write(5,*) iter,iter*tcf, np, Drug_Released 
  IF (np.EQ.0)THEN
     STOP
  ENDIF
ENDDO

stop
END
