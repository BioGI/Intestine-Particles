!===================================================================================================
! This subroutine randomly creates the particle's locations  
!===================================================================================================
USE SetPrecision
USE Setup
USE Geometry
USE ParticleDrug
IMPLICIT NONE

REAL(dbl) :: D_Couette, L_Couette, Dx_Couette, Dy_Couette, Dz_Couette
REAL(dbl) :: x_center, y_center, z_center
REAL(dbl) :: rMax, teta1Max, teta2Max, rr, teta1, teta2
REAL(dbl) :: x_particle, y_particle, z_particle
REAL(dbl) :: R_Particle, R_Par_Max, D_Par_Max, R_Boundary
REAL(dbl) :: R_left, R_right, dz, Volume, Area
REAL(dbl) :: xp,yp,zp,up,vp,wp,U_slip,sh_conf,Sh_shear,Sh_slip,rp,bulk_conc,delNBbyCV,SSt,S,par_pH,par_conc 
REAL(dbl) :: eps,S_ratio,Sh,tmp,Cb,deltaR,zcf3,Drug_Released
REAL(dbl) :: Drug_Released_del_diff, Drug_Released_del_shear, Drug_Released_del_slip
REAL(dbl) :: Pw,Gut_Surface,Gut_volume,Drug_Absorbed_Total
INTEGER   :: Continuation,Counter,n,num_Par, i, j, k, z_left, z_right, CPU, parid, cur_part
CHARACTER(7):: TEMP,iter_char                                       ! iteration stored as a character
CHARACTER(5):: sub_char,tmp_char
TYPE(ParRecord), POINTER :: CurPar
TYPE(ParRecord), POINTER :: currentt
TYPE(ParRecord), POINTER :: nextt

CALL ReadInput
Pw   =1.0e-4                          ! Permeability
Gut_surface= 7.0738                   ! cm2
Gut_volume = 1.9958                   ! cm3
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

Write(*,*) 'Please enter 0 if using particle.dat, 1 if using pardat files:'
read(*,*)  Continuation
Write(*,*) 'Please enter the ieration number:'
read(*,*) iter

IF (Continuation.EQ.1) THEN
   Write(*,*) 'Please enter the number of particles:'
   read(*,*) np
   Write(*,*) 'Please enter the Released Drug:'
   read(*,*) Drug_Released
ENDIF


CALL list_init(ParListHead)
CurPar => ParListHead

IF (Continuation.EQ.1) THEN
   WRITE(iter_char(1:7),'(I7.7)') iter
   WRITE(sub_char(1:5),'(I5.5)') CPU
   OPEN(160,FILE='pardat-'//iter_char//'-'//sub_char//'.csv')
   READ(160,*) tmp_char
   DO i = 1, np
      write(*,*) i 
      READ(160,*) xp,yp,zp,up,vp,wp,U_slip,parid,sh_conf,Sh_shear,Sh_slip,rp,bulk_conc,delNBbyCV,SSt,S,par_pH,par_conc,cur_part 
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
      CurPar%next%pardata%par_pH   =  par_pH
      CurPar%next%pardata%par_conc =  par_conc
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
ELSE
   OPEN(160,FILE='particle.dat')
   READ(160,*) np
   CALL list_init(ParListHead)
   CurPar => ParListHead
   DO i = 1, np
      READ(160,*) parid,xp,yp,zp,rp                                      ! read particle.dat file
      CALL list_init(CurPar%next)               
      CurPar%next%prev => CurPar
      CurPar%next%next => null()
      CurPar%next%pardata%parid = parid
      CurPar%next%pardata%xp = xp
      CurPar%next%pardata%yp = yp
      CurPar%next%pardata%zp = zp
      CurPar%next%pardata%up = 0.0_dbl
      CurPar%next%pardata%vp = 0.0_dbl
      CurPar%next%pardata%wp = 0.0_dbl
      CurPar%next%pardata%rp = rp
      CurPar%next%pardata%xpold = CurPar%next%pardata%xp
      CurPar%next%pardata%ypold = CurPar%next%pardata%yp
      CurPar%next%pardata%zpold = CurPar%next%pardata%zp
      CurPar%next%pardata%upold = CurPar%next%pardata%up
      CurPar%next%pardata%vpold = CurPar%next%pardata%vp
      CurPar%next%pardata%wpold = CurPar%next%pardata%wp
      CurPar%next%pardata%rpold = CurPar%next%pardata%rp
      CurPar%next%pardata%par_pH     = 0.0_dbl
      CurPar%next%pardata%par_conc   = 0.0_dbl
      CurPar%next%pardata%gamma_cont = 0.0000_dbl
      CurPar%next%pardata%sh_conf  = 0.0_dbl
      CurPar%next%pardata%sh_shear = 0.0_dbl
      CurPar%next%pardata%sh_slip  = 0.0_dbl
      CurPar%next%pardata%S        = 0.0_dbl
      CurPar%next%pardata%Sst      = 0.0_dbl
      CurPar%next%pardata%Veff     = 0.0_dbl
      CurPar%next%pardata%Nbj      = 0.0_dbl
      CurPar%next%pardata%bulk_conc= 0.0000_dbl
      CurPar%next%pardata%delNBbyCV= 0.00000_dbl
      CurPar%next%pardata%cur_part = 1
      CurPar%next%pardata%new_part = 1
      CurPar => CurPar%next
   END DO
   CLOSE(160)

ENDIF

!--- Shear ------------------
S=0.0_dbl
Counter = 0
IF (Continuation.EQ.1) THEN
   currentt => ParListHead%next
   DO WHILE (ASSOCIATED(currentt))
      nextt => currentt%next
      S= S + currentt%pardata%S
      Counter=Counter + 1
      currentt => nextt
   ENDDO
   S = S/Counter
ELSE
   WRITE(*,*) 'Please enter the average shear:'
   READ(*,*) S
ENDIF

!--- Cb --------------------
Cb=0.0_dbl
!---------------------------
Drug_Initial=0
Drug_Released_del_diff=0
Drug_Released_del_shear=0
Drug_Released_del_slip=0
Drug_Released_Total=1.0e-16
Drug_Absorbed=0
Drug_Absorbed_Total=0
Drug_Remained_in_Domain=0
Drug_Loss_Percent=0

OPEN (5,file='Drug-Conservation-00001-Extrapolated.dat')
WRITE(5,'(A145)') '#VARIABLES =iter,time, Initial, Released_del_diff, Released_del_shear, Released_del_slip, Released_Total, Absorbed, Remained_in_Domain, Loss_Percent'

Do i=iter, 2
   np=0
   iter=iter+1

   currentt => ParListHead%next
   DO WHILE (ASSOCIATED(currentt))
      nextt => currentt%next
      IF (currentt%pardata%rp.GT.1e-16) THEN
         !---- Shear Effects -------------------------------------------------------------
         currentt%pardata%S=S
         !S=currentt%pardata%S
         Sst= S* (currentt%pardata%rp**2.0) / diffm
         currentt%pardata%Sst= Sst
         IF (Sst .LE. 5.0_dbl) THEN
            currentt%pardata%sh_shear = 1.0_dbl + 0.281_dbl*(Sst**0.5_dbl)          -1.0_dbl
         ELSE IF ((Sst .GT. 5.0_dbl).AND.(Sst .LE. 80.0)) THEN     
            currentt%pardata%sh_shear = 1.181_dbl*(Sst**0.2_dbl)                    -1.0_dbl
         ELSE IF (Sst.GT.80.0) THEN
            currentt%pardata%sh_shear = 4.5_dbl - (7.389/(Sst**(1.0_dbl/3.0_dbl)) ) -1.0_dbl 
         ENDIF
      ENDIF
      currentt => nextt
   ENDDO

   CALL  Compute_C_surface

   currentt => ParListHead%next
   DO WHILE (ASSOCIATED(currentt))
      nextt => currentt%next
      IF (currentt%pardata%rp.GT.1e-16) THEN
         np= np+1
         currentt%pardata%bulk_conc =Cb
         currentt%pardata%rpold = currentt%pardata%rp
         Sh= 1.0 + currentt%pardata%sh_shear  
         eps= (molarvol*diffm*tcf*Sh) * (currentt%pardata%par_conc- currentt%pardata%bulk_conc)
         tmp= (currentt%pardata%rpold)**2.0_dbl - 4.0_dbl*eps 
         If (tmp .LT. 0.0_dbl)THEN
            tmp=0.0_dbl
         ENDIF
         currentt%pardata%rp= 0.5_dbl * (currentt%pardata%rpold+ sqrt(tmp))
         deltaR= currentt%pardata%rpold-currentt%pardata%rp
         currentt%pardata%delNBbyCV = (4.0_dbl/3.0_dbl) * PI*(currentt%pardata%rpold**3.0_dbl - currentt%pardata%rp**3.0_dbl) /(molarvol*zcf3)
         Drug_Released_Total= Drug_Released_Total + currentt%pardata%delNBbyCV *1000000.0_dbl*zcf3
     ENDIF
     currentt => nextt
  ENDDO    
  !Drug_Remained_in_Domain = Drug_Remained_in_Domain + Drug_Released_Total 
  Cb= Drug_Remained_in_Domain/Gut_volume 
  Drug_Absorbed= Cb * Pw *Gut_surface*tcf
  Drug_Absorbed_Total=Drug_Absorbed_Total + Drug_Absorbed
  Drug_Remained_in_Domain = Drug_Released_Total- Drug_Absorbed_Total 
  Drug_Loss_Percent=  100.0*(Drug_Released_Total-Drug_Remained_in_Domain-Drug_Absorbed_Total)/Drug_Released_Total 


  WRITE(5,'(I8, F13.4, 7E18.10, F11.6)') iter, iter*tcf, Drug_Initial, Drug_Released_del_diff, Drug_Released_del_shear, Drug_Released_del_slip, Drug_Released_Total, Drug_Absorbed_Total, Drug_Remained_in_Domain, Drug_Loss_Percent

  IF (np.EQ.0)THEN
     STOP
  ENDIF

  Output_Intervals=453
  IF ((MOD(iter, Output_Intervals) .EQ. 0))  THEN
     WRITE(iter_char(1:7),'(I7.7)') iter
     OPEN(170,FILE='pardat-'//iter_char//'-'//sub//'.csv')
     WRITE(170,*) '"x","y","z","u","v","w","U_slip", "ParID","Sh_conf","Sh_shear","Sh_slip","rp","Cb/Cs","delNBbyCV","Sst","S","C_surface","CPU"'
     currentt => ParListHead%next                                                       
     DO WHILE (ASSOCIATED(currentt))
        nextt => currentt%next   
        IF (currentt%pardata%rp .GT. Min_R_Acceptable) THEN                                           ! only write particle data when particle is not fully dissolved
            WRITE(170,1001) currentt%pardata%xp                   ,',', &
                            currentt%pardata%yp                   ,',', &
                            currentt%pardata%zp                   ,',', &
                            currentt%pardata%up*vcf               ,',', &
                            currentt%pardata%vp*vcf               ,',', &
                            currentt%pardata%wp*vcf               ,',', &
                            currentt%pardata%U_slip               ,',', &
                            currentt%pardata%parid                ,',', &
                            currentt%pardata%sh_conf              ,',', &
                            currentt%pardata%sh_shear             ,',', &
                            currentt%pardata%sh_slip              ,',', &
                            currentt%pardata%rp                   ,',', &
                            currentt%pardata%bulk_conc/S_intrinsic,',', &
                            currentt%pardata%delNBbyCV            ,',', &
                            currentt%pardata%Sst                  ,',', &
                            currentt%pardata%S                    ,',', &
                            currentt%pardata%par_pH               ,',', &
                            currentt%pardata%par_conc             ,',', &
                            currentt%pardata%cur_part
        END IF 
1001    format (6(F9.4,a2),(E18.9,a2),(I6,a2),3(F12.8,a2),3(F11.8,a2),3(F13.8,a2),I4)
        currentt => nextt
     ENDDO
     CLOSE(170)
  ENDIF
ENDDO

stop
END

