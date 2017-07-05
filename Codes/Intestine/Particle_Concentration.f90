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
REAL(dbl) :: R_left, R_right, dz, Volume, Area, Counter
REAL(dbl) :: xp,yp,zp,up,vp,wp,U_slip,sh_conf,Sh_shear,Sh_slip,rp,bulk_conc,delNBbyCV,SSt,S,par_Conc 
REAL(dbl) :: eps,S_ratio,Sh,tmp,Cb,deltaR,zcf3,Drug_Released
REAL(dbl) :: V_eff,L_eff,R_eff,R_Par,R_node,Vol_Percent
REAL(dbl) :: Md
REAL,DIMENSION(300,300,300) :: Node_N,Par_N,Par_Dose,Par_Con_N,Par_Con_Dose,phii
REAL(dbl) :: Par_Con_N_AVE,Par_Con_Dose_AVE,phii_AVE
REAL(dbl) :: Par_Con_N_VARIANCE,Par_Con_Dose_VARIANCE,phii_VARIANCE
REAL(dbl) :: Corr_ParN_ParDose,Corr_ParN_phii,Corr_ParDose_phii
REAL(dbl) :: Covar_ParN_ParDose,Covar_ParN_phii,Covar_ParDose_phii
INTEGER   :: NN,N_CPU, Lines_N,Fluid_N
INTEGER   :: n, num_Par
INTEGER   :: io,i, j, k,kk, i0,i1,j0,j1,k0,k1,iii, jjj, kkk 
INTEGER   :: q, z_left, z_right, CPU, parid, cur_part
CHARACTER(7):: TEMP,iter_char                                       ! iteration stored as a character
CHARACTER(5):: sub_char,tmp_char
TYPE(ParRecord), POINTER :: CurPar
TYPE(ParRecord), POINTER :: current_
TYPE(ParRecord), POINTER :: next_

Md=206.285

CALL ReadInput

nuL  = (2.0_dbl*tau -1.0_dbl)/6.0_dbl	! lattice kinematic viscosity
denL = 1.0_dbl		                	  ! arbitrary lattice density (1.0 for convenience)
kMin = 1
kMax = nz
nzSub= nz
numprocs=1
CPU=1
S_ratio =2.30196707    ! For no buffer case only
zcf=L/nz
zcf3=1.0e6*zcf*zcf*zcf !cm3
write(*,*) 'zcf3=',zcf3
CALL Global_Setup
CALL AdvanceGeometry
r(0:nzSub+1) = rDom(kMin-1:kMax+1)  


!--- Computing the total number  of fluid nodes ----------------------------------------------------
Fluid_N=0
DO i=1,nx
   DO j=1,ny
      Do k=1,nz
         IF (node(i,j,k) .EQ. FLUID) THEN
            Fluid_N= Fluid_N+1
         ENDIF  
      ENDDO
   ENDDO
ENDDO
WRITE(*,*) 'Total number of fluid nodes:',Fluid_N

!--- Inputs from the user --------------------------------------------------------------------------
WRITE(*,*) 'Please enter the iteration:'
read(*,*) iter
WRITE(*,*) 'Please enter the efective volume around the nodes (% of the total volume):'
read(*,*) Vol_Percent

!--- Dealing with effective volume around each node ------------------------------------------------
Vol_Percent=Vol_Percent/100
V_eff= Fluid_N*Vol_Percent
L_eff=V_eff**(1.0_dbl/3.0_dbl)
R_eff=((3.0_dbl/4.0_dbl)*V_eff/PI)**(1.0_dbl/3.0_dbl)
WRITE(*,*) 'V_eff,L_eff,R_eff',V_eff,L_eff,R_eff

!--- Reading the pardat file ------------------------------------------------------------------------
WRITE(iter_char(1:7),'(I7.7)') iter
WRITE(sub_char(1:5),'(I5.5)') CPU
OPEN(160,FILE='pardat-'//iter_char//'-'//sub_char//'.csv')
np=0
READ(160,*,iostat=io)
DO       
   READ(160,*,iostat=io)
   if (io/=0) exit
   np=np +1
ENDDO
WRITE(*,*) "number of particles calculated",np
REWIND(160)

CALL list_init(ParListHead)
CurPar => ParListHead
READ(160,*) 
DO i = 1, np
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

!--- Reading the output files to extract teh concentration data ------------------------------------
WRITE(*,*)'Please enter the number of output files to read (CPUs)'
READ(*,*) N_CPU
DO NN=1,N_CPU
   WRITE(sub_char(1:5),'(I5.5)') NN
   OPEN(60,FILE='out-'//iter_char//'-'//sub_char//'.dat')

   !--- Finding the number of lines in each subdomain output file
   READ(60,*)
   READ(60,*)
   Lines_N=0
   DO       
      READ(60,*,iostat=io)
      if (io/=0) exit
      Lines_N=Lines_N +1
   ENDDO
   WRITE(*,*) NN, Lines_N
   !--- Reading the output files
   REWIND(60)
   READ(60,*)
   READ(60,*)
   DO kk=1,Lines_N
      READ(60,*)i,j,k, phii(i,j,k), Node(i,j,k)
   ENDDO
   CLOSE(60)
ENDDO



DO i=1,nx
   DO j=1,ny
      Do k=1,nz
         IF (node(i,j,k) .EQ. FLUID) THEN
            !--- Loop over Particels ----------------------------------------------------
            Par_N(i,j,k) = 0
            Par_Dose(i,j,k) = 0.0_dbl
            current_ => ParListHead%next
            DO WHILE (ASSOCIATED(current_))
               next_ => current_%next
               IF (current_%pardata%rp .GT. 1e-16) THEN
                  R_Par=sqrt((current_%pardata%xp-i)**2.0_dbl+ (current_%pardata%yp-j)**2.0_dbl +(current_%pardata%zp-k)**2.0_dbl)
                  IF (R_Par.LT.R_eff)THEN
                     Par_Dose(i,j,k)=Par_Dose(i,j,k)+((4.0_dbl/3.0_dbl)*PI*(100.0_dbl*current_%pardata%rp)**3.0_dbl)*Md/molarvol
                     Par_N(i,j,k)=Par_N(i,j,k)+ 1
                  ENDIF
               ENDIF
               current_ => next_                  
            ENDDO  

            !--- Loop over surrounding nodes -------------------------------------------- 
            q=ceiling(R_eff) 
            i0=i-q
            i1=i+q
            j0=j-q
            j1=j+q
            k0=k-q
            k1=k+q
            iF (i0.LT.1)  i0=1
            iF (i1.GT.nx) i1=nx
            iF (j0.LT.1)  j0=1
            iF (j1.GT.ny) j1=ny
            iF (k0.LT.1)  k0=1
            iF (k1.GT.nz) k1=nz

            Node_N(i,j,k)= 0
            DO iii= i0,i1  
               DO jjj= j0,j1   
                  DO kkk= k0,k1  
                     R_Node=sqrt((iii-i)**2.0_dbl+ (jjj-j)**2.0_dbl +(kkk-k)**2.0_dbl)
                     IF (R_Node.LT.R_eff) THEN
                        IF (node(iii,jjj,kkk) .EQ. FLUID) THEN
                           Node_N(i,j,k)=Node_N(i,j,k) +1
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO   

            !--- Compute the particle concentrations ------------------------------------  
            Par_Con_N(i,j,k)   =Par_N(i,j,k)/   (Node_N(i,j,k) *zcf3)
            Par_Con_Dose(i,j,k)=Par_Dose(i,j,k)/(Node_N(i,j,k) *zcf3)
         ENDIF   
      ENDDO
   ENDDO    
ENDDO

!--- Outputing the particle concentration data for visualization, and averaging over the domain ---------------
OPEN(60,FILE='out-par-con-'//iter_char//'-'//sub//'.dat')
WRITE(60,*) 'VARIABLES = "x" "y" "z" "Node" "Node_N" "Par_N" "Par_Con_N" "Par_Con_Dose" "phi"'
WRITE(60,'(A10,E15.5,A5,I4,A5,I4,A5,I4,A8)') 'ZONE T="',iter/(nt/nPers),'" I=',nx,' J=',ny,' K=',nz,'F=POINT'

Par_Con_N_AVE =0
Par_Con_Dose_AVE=0
phii_AVE =0
DO k=1,nz
   DO j=1,ny
      DO i=1,nx
         IF (node(i,j,k) .EQ. FLUID) THEN
            Par_Con_N_AVE    = Par_Con_N_AVE    + Par_Con_N(i,j,k)
            Par_Con_Dose_AVE = Par_Con_Dose_AVE + Par_Con_Dose(i,j,k)
            phii_AVE         = phii_AVE         + phii(i,j,k)
            WRITE(60,'(4I4,2F12.2,3E16.6)') i,j,k,node(i,j,k),Node_N(i,j,k)/Fluid_N,Par_N(i,j,k), Par_Con_N(i,j,k),Par_Con_Dose(i,j,k),phii(i,j,k)
         ELSE
            WRITE(60,'(9I4)') i,j,k,node(i,j,k),0,0,0,0,0
         END IF
      ENDDO
   ENDDO
ENDDO
Par_Con_N_AVE    = Par_Con_N_AVE    / Fluid_N
Par_Con_Dose_AVE = Par_Con_Dose_AVE / Fluid_N
phii_AVE         = phii_AVE         / Fluid_N

!--- Statistical Analysis -----------------------------------------------------------------------------------------
!--- Correlation & Dependance: https://en.wikipedia.org/wiki/Correlation_and_dependence
!--- Standard Deviation:       https://en.wikipedia.org/wiki/Standard_deviation#Corrected_sample_standard_deviation
!------------------------------------------------------------------------------------------------------------------
Par_Con_N_VARIANCE    = 0.0_dbl
Par_Con_Dose_VARIANCE = 0.0_dbl
phii_VARIANCE         = 0.0_dbl
Covar_ParN_ParDose    = 0.0_dbl
Covar_ParN_phii       = 0.0_dbl
Covar_ParDose_phii    = 0.0_dbl

DO k=1,nz
   DO j=1,ny
      DO i=1,nx
         IF (node(i,j,k) .EQ. FLUID) THEN
             Par_Con_N_VARIANCE    = Par_Con_N_VARIANCE    + (Par_Con_N(i,j,k)-Par_Con_N_AVE)**2.0_dbl
             Par_Con_Dose_VARIANCE = Par_Con_Dose_VARIANCE + (Par_Con_Dose(i,j,k)-Par_Con_Dose_AVE)**2.0_dbl
             phii_VARIANCE         = phii_VARIANCE         + (phii(i,j,k)-phii_AVE)**2.0_dbl
             Covar_ParN_ParDose    = Covar_ParN_ParDose    + (Par_Con_N(i,j,k)-Par_Con_N_AVE)       * (Par_Con_Dose(i,j,k)-Par_Con_Dose_AVE)
             Covar_ParN_phii       = Covar_ParN_phii       + (Par_Con_N(i,j,k)-Par_Con_N_AVE)       * (phii(i,j,k)-phii_AVE)
             Covar_ParDose_phii    = Covar_ParDose_phii     +(Par_Con_Dose(i,j,k)-Par_Con_Dose_AVE) * (phii(i,j,k)-phii_AVE)
         END IF
      ENDDO
   ENDDO
ENDDO
Corr_ParN_ParDose= Covar_ParN_ParDose / (sqrt(Par_Con_N_VARIANCE   *Par_Con_Dose_VARIANCE))
Corr_ParN_phii   = Covar_ParN_phii    / (sqrt(Par_Con_N_VARIANCE   *phii_VARIANCE        ))
Corr_ParDose_phii= Covar_ParDose_phii / (sqrt(Par_Con_Dose_VARIANCE*phii_VARIANCE        ))

Par_Con_N_VARIANCE    = Par_Con_N_VARIANCE   /(Fluid_N-1)
Par_Con_Dose_VARIANCE = Par_Con_Dose_VARIANCE/(Fluid_N-1)
phii_VARIANCE         = phii_VARIANCE        /(Fluid_N-1)

WRITE(*,*)'Average:    Dose, N, phi',Par_Con_Dose_AVE, Par_Con_N_AVE,phii_AVE
WRITE(*,*)'Variance:   Dose, N, phi',Par_Con_N_VARIANCE,Par_Con_Dose_VARIANCE,phii_VARIANCE
WRITE(*,*)'Covariance: N_Dose, N_phi, Dose_phi',Covar_ParN_ParDose, Covar_ParN_phii,Covar_ParDose_phii 
WRITE(*,*)'Correlation:N_Dose, N_phi, Dose_phi',Corr_ParN_ParDose, Corr_ParN_phii,Corr_ParDose_phii 

END

