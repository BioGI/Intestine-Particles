!===================================================================================================
!This subroutine creates the particle distribution list 
!===================================================================================================
IMPLICIT NONE
INTEGER			:: Flag_Normal_Dist, Counter,i,j
REAL*8      :: Dose, nu_m, M_d, dmin,dmax,sigma,miu,pi
REAL*8      :: d_end,d_up,d1,d2,eps,V_integ,V_par_i,Vf,Vf2,Vp_tot,Vp_tot_Achieved, Rlist(10000000)

Write(*,*) 'Please enter the dosage (\mu g):'
read(*,*) Dose
Write(*,*) 'Please enter the maximum cut-off diameter, D_max (micron):'
read(*,*) dmax
Write(*,*) 'Please enter the average particle diameter, D_ave (micron):'
read(*,*) miu
Write(*,*) 'Please enter the minimum cut-off diameter, D_min (micron):'
read(*,*) dmin
Write(*,*) 'Please enter the standard deviation of the particle distribution, Sigma (micron):'
read(*,*) sigma
Write(*,*) 'Please enter 1 for Normal_Distribution and 0 for Modified_Normal_Distribution:'
read(*,*) Flag_Normal_Dist

pi              = 4.0 * atan(1.0)
nu_m            = 268.0e6  	           ! (\mu m)^3 / (\mu mol)  		
M_d             = 206.285              ! (\mu g)   / (\mu mol)
Vp_tot_Achieved = 0.0
Vp_tot          = (Dose/M_d) * nu_m    ! (\mu m) ^3
d_up            = dmax
d1              = dmax
Counter         = 0
eps=1e-6


open(49,file='Particle_Dist_info.dat')
DO WHILE (d1.GE.dmin) 
   Counter = Counter + 1
   V_integ = 0.0
   V_par_i = 1.0
   d2      = d_up
   DO WHILE (V_par_i .GE. V_integ)
!      eps = eps + 1e-11  
      d1 =  d2 - eps
      IF (Flag_Normal_Dist .EQ. 1) THEN
         V_integ = V_integ + (1.0/sqrt(2.*pi*sigma**2)) * exp(-(((d2+d1)/2.0)-miu)**2/(2.*sigma**2)) * (d2-d1) *Vp_tot
      ELSE IF (Flag_Normal_Dist .EQ. 0) THEN
         V_integ = V_integ + ((1.0/sqrt(2.*pi*sigma**2)) * exp(-(((d2+d1)/2.0)-miu)**2/(2.*sigma**2))) *( 0.5*(1.0+tanh(0.2*((d2+d1)/2.0) -3.0))) *((d2-d1) *Vp_tot)
      ELSE
         write(*,*) 'wrong choice for choosing between Normal_Distribution and Modified_Normal_Distribution'
      END IF    
      V_par_i = (4.0/3.0)* PI* ((d1+d2)/4.0)**3.0
      d2 = d1
   END DO
   Rlist(Counter)= (d1+d2)/2.0
   Vf = (V_par_i/Vp_tot)/(d_up-d1) 
   Vf2= V_par_i/Vp_tot 
   Vp_tot_Achieved= Vp_tot_Achieved + V_Par_i
   write(49,*) Counter, Rlist(Counter),  V_integ, V_par_i, Vf,vf2
   d_up= d1
END DO
Close(49)

open(50,file='Particle_Sizes.dat')
write(50,*) Counter 
Do i=1,Counter
   write(50,*) Rlist(i)
END DO
Close(50)

Write(*,*) 'Total dosage                 ', Dose
write(*,*) 'Total Particle volume desired', Vp_tot 
write(*,*) 'Minimum cut-off diameter     ', dmin 
write(*,*) 'Average particle diameter    ', miu
Write(*,*) 'Maximum cut-off diameter     ', dmax
write(*,*) 'Number of particle generated:', Counter
write(*,*) 'Total particle volume desired and achieved', Vp_tot, Vp_tot_Achieved
write(*,*) 'Achieved Dosage', (1.0- ((Vp_tot-Vp_tot_Achieved)/Vp_tot)) * Dose

stop

end
